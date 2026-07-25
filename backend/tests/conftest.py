"""Shared test fixtures. Every pipeline-stage test runs against a throwaway in-memory SQLite
database - never the real gene_doctor.db.
"""

from __future__ import annotations

import pytest
from sqlalchemy import create_engine
from sqlalchemy.orm import Session

from backend.database.models import (
    AiPrediction,
    AnalysisRun,
    Annotation,
    Base,
    EvidenceWeights,
    GenomeUpload,
    Variant,
)


@pytest.fixture
def session() -> Session:
    engine = create_engine("sqlite://")
    Base.metadata.create_all(engine)
    with Session(engine) as s:
        yield s


@pytest.fixture
def analysis_run(session: Session) -> AnalysisRun:
    upload = GenomeUpload(file_path="/tmp/fake.fa", status="valid", progress=1.0)
    session.add(upload)
    session.flush()
    run = AnalysisRun(genome_upload_id=upload.id, status="queued")
    session.add(run)
    session.add(EvidenceWeights(id=1))
    session.commit()
    return run


def make_variant(
    session: Session,
    run: AnalysisRun,
    chrom: str = "1",
    pos: int = 1000,
    ref: str = "A",
    alt: str = "G",
    gene: str | None = "BRCA1",
    transcript_id: str | None = None,
) -> Variant:
    variant = Variant(
        analysis_run_id=run.id,
        chrom=chrom,
        pos=pos,
        ref=ref,
        alt=alt,
        gene=gene,
        transcript_id=transcript_id,
    )
    session.add(variant)
    session.commit()
    return variant


def add_annotation(
    session: Session,
    variant: Variant,
    source: str = "clinvar",
    disease: str = "breast cancer",
    evidence_strength: float = 1.0,
    clinical_significance: str | None = "Pathogenic",
) -> Annotation:
    annotation = Annotation(
        variant_id=variant.id,
        source=source,
        disease=disease,
        clinical_significance=clinical_significance,
        evidence_strength=evidence_strength,
    )
    session.add(annotation)
    session.commit()
    return annotation


def add_ai_prediction(
    session: Session,
    variant: Variant,
    source: str = "variant_lm",
    functional_impact: float | None = 0.8,
    regulatory_impact: float | None = 0.5,
    pathogenicity_probability: float | None = 0.7,
    confidence: float | None = 0.6,
) -> AiPrediction:
    prediction = AiPrediction(
        variant_id=variant.id,
        source=source,
        functional_impact=functional_impact,
        regulatory_impact=regulatory_impact,
        pathogenicity_probability=pathogenicity_probability,
        confidence=confidence,
    )
    session.add(prediction)
    session.commit()
    return prediction
