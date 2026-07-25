"""Pipeline orchestrator: drives an analysis run through alignment -> variant calling -> annotation
-> AI prediction -> evidence aggregation -> disease ranking -> transcript discovery -> ASO design.

Each stage updates the AnalysisRun row so the UI can poll progress. A failure in any stage marks the
run failed with the error message rather than crashing the server. The AI-prediction stage is
explicitly non-fatal: the spec requires the rest of the pipeline to still deliver validated
(ClinVar/GWAS) results when the model is unavailable.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

from sqlalchemy import create_engine
from sqlalchemy.orm import Session, sessionmaker

from backend.annotation.annotate import annotate_variants
from backend.annotation.genes import assign_genes
from backend.aso.designer import design_asos
from backend.database.models import AnalysisRun, GenomeUpload
from backend.genomics.align import align_genome
from backend.genomics.load_vcf import load_vcf
from backend.genomics.variant_call import call_variants
from backend.services.disease_ranker import rank_diseases
from backend.services.evidence_engine import score_variants
from backend.settings import REPO_ROOT, settings
from backend.transcript.transcript_ranker import discover_transcripts
from backend.variant_lm.predictor import predict_variants

logger = logging.getLogger(__name__)

WORK_DIR = REPO_ROOT / "data" / "runs"


@dataclass
class StageResult:
    variants: int = 0
    annotations: int = 0
    predictions: int = 0
    scores: int = 0
    diseases: int = 0
    transcripts: int = 0
    asos: int = 0


def _sync_sessionmaker() -> sessionmaker[Session]:
    """The pipeline runs in a worker thread, so it uses a synchronous engine of its own rather than
    the app's async session.
    """
    sync_url = settings.database_url.replace("sqlite+aiosqlite", "sqlite")
    return sessionmaker(create_engine(sync_url))


def _set_stage(session: Session, run_id: str, stage: str) -> None:
    run = session.get(AnalysisRun, run_id)
    run.status = stage
    run.current_stage = stage
    session.commit()
    logger.info("Analysis %s: stage %s", run_id, stage)


def run_pipeline(analysis_run_id: str) -> StageResult:
    """Execute every stage for one analysis run. Safe to call from a worker thread."""
    Session_ = _sync_sessionmaker()
    result = StageResult()
    work_dir = WORK_DIR / analysis_run_id

    with Session_() as session:
        run = session.get(AnalysisRun, analysis_run_id)
        upload = session.get(GenomeUpload, run.genome_upload_id)
        genome_path = upload.file_path

    try:
        with Session_() as session:
            _set_stage(session, analysis_run_id, "aligning")
        bam_path = align_genome(genome_path, work_dir)

        with Session_() as session:
            _set_stage(session, analysis_run_id, "calling")
        vcf_path = call_variants(bam_path, work_dir)

        with Session_() as session:
            result.variants = load_vcf(session, analysis_run_id, vcf_path)

            _set_stage(session, analysis_run_id, "annotating")
            # Gene/transcript assignment first: it is coordinate-based, and every later stage keys
            # off Variant.gene, including for variants that carry no clinical evidence.
            assign_genes(session, analysis_run_id)
            counts = annotate_variants(session, analysis_run_id)
            result.annotations = counts.clinvar + counts.gwas

            _set_stage(session, analysis_run_id, "predicting")
            try:
                result.predictions = predict_variants(session, analysis_run_id)
            except Exception:
                # Non-fatal by design: validated evidence must still reach the user.
                logger.exception("AI prediction stage failed for %s; continuing", analysis_run_id)

            _set_stage(session, analysis_run_id, "aggregating")
            result.scores = score_variants(session, analysis_run_id)

            _set_stage(session, analysis_run_id, "ranking")
            result.diseases = rank_diseases(session, analysis_run_id)

            _set_stage(session, analysis_run_id, "designing_aso")
            result.transcripts = discover_transcripts(session, analysis_run_id)
            result.asos = design_asos(session, analysis_run_id)

            run = session.get(AnalysisRun, analysis_run_id)
            run.status = "done"
            run.current_stage = None
            session.commit()

    except Exception as exc:
        logger.exception("Pipeline failed for %s", analysis_run_id)
        with Session_() as session:
            run = session.get(AnalysisRun, analysis_run_id)
            run.status = "failed"
            run.error = str(exc)
            session.commit()
        raise

    return result
