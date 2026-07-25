"""Tests for coordinate-based gene assignment.

This is what makes a novel variant visible to the rest of the pipeline: gene symbols used to be a
side effect of ClinVar/GWAS annotation, so any variant without published evidence had gene=None and
was silently dropped before transcript discovery.
"""

from __future__ import annotations

from sqlalchemy.orm import Session

from backend.annotation.genes import assign_genes
from backend.database.models import AnalysisRun, RefEnsemblTranscript
from backend.tests.conftest import make_variant

EXONS = [
    {"start": 1000, "end": 1100, "exon_number": "1"},
    {"start": 2000, "end": 2100, "exon_number": "2"},
]


def _add_ref(
    session: Session,
    gene: str,
    transcript_id: str,
    chrom: str = "17",
    exons: list[dict] | None = None,
) -> RefEnsemblTranscript:
    exons = exons if exons is not None else EXONS
    ref = RefEnsemblTranscript(
        gene=gene,
        ensembl_gene_id=f"ENSG_{gene}",
        ensembl_transcript_id=transcript_id,
        chrom=chrom,
        start=exons[0]["start"],
        end=exons[-1]["end"],
        strand="+",
        exons=exons,
    )
    session.add(ref)
    session.commit()
    return ref


def test_variant_with_no_evidence_still_gets_a_gene(
    session: Session, analysis_run: AnalysisRun
) -> None:
    variant = make_variant(session, analysis_run, chrom="17", pos=2050, gene=None)
    _add_ref(session, "HTT", "ENST_HTT")

    assert assign_genes(session, analysis_run.id) == 1
    assert variant.gene == "HTT"
    assert variant.transcript_id == "ENST_HTT"


def test_intronic_position_inside_the_locus_still_maps(
    session: Session, analysis_run: AnalysisRun
) -> None:
    variant = make_variant(session, analysis_run, chrom="17", pos=1500, gene=None)
    _add_ref(session, "HTT", "ENST_HTT")

    assign_genes(session, analysis_run.id)
    assert variant.gene == "HTT"


def test_position_outside_every_transcript_is_left_alone(
    session: Session, analysis_run: AnalysisRun
) -> None:
    variant = make_variant(session, analysis_run, chrom="17", pos=99_000, gene=None)
    _add_ref(session, "HTT", "ENST_HTT")

    assert assign_genes(session, analysis_run.id) == 0
    assert variant.gene is None


def test_wrong_chromosome_does_not_match(session: Session, analysis_run: AnalysisRun) -> None:
    variant = make_variant(session, analysis_run, chrom="4", pos=2050, gene=None)
    _add_ref(session, "HTT", "ENST_HTT", chrom="17")

    assert assign_genes(session, analysis_run.id) == 0
    assert variant.gene is None


def test_exonic_transcript_beats_one_that_leaves_it_intronic(
    session: Session, analysis_run: AnalysisRun
) -> None:
    """The exonic placement is the one that matters for ASO design."""
    variant = make_variant(session, analysis_run, chrom="17", pos=2050, gene=None)
    # Spans the position but places it between exons.
    _add_ref(
        session,
        "INTRONIC",
        "ENST_INTRONIC",
        exons=[
            {"start": 1000, "end": 1100, "exon_number": "1"},
            {"start": 3000, "end": 3100, "exon_number": "2"},
            {"start": 4000, "end": 4100, "exon_number": "3"},
            {"start": 5000, "end": 5100, "exon_number": "4"},
        ],
    )
    _add_ref(session, "EXONIC", "ENST_EXONIC")

    assign_genes(session, analysis_run.id)
    assert variant.gene == "EXONIC"


def test_fuller_isoform_wins_a_tie(session: Session, analysis_run: AnalysisRun) -> None:
    variant = make_variant(session, analysis_run, chrom="17", pos=2050, gene=None)
    _add_ref(session, "FRAGMENT", "ENST_FRAGMENT", exons=[EXONS[1]])
    _add_ref(session, "FULL", "ENST_FULL", exons=[*EXONS, {"start": 2500, "end": 2600, "exon_number": "3"}])

    assign_genes(session, analysis_run.id)
    assert variant.gene == "FULL"
