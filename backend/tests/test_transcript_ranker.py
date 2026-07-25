"""Stage 9 tests: exon/variant geometry -> mechanisms, and transcript ranking."""

from __future__ import annotations

from sqlalchemy import select
from sqlalchemy.orm import Session

from backend.database.models import (
    AnalysisRun,
    Disease,
    DiseaseVariantLink,
    RefEnsemblTranscript,
    Transcript,
    Variant,
    VariantScore,
)
from backend.tests.conftest import add_ai_prediction, add_annotation, make_variant
from backend.transcript.transcript_ranker import discover_transcripts

# A 3-exon toy gene on chromosome 17. Exon 2 is internal, so it is skippable/re-includable.
EXONS = [
    {"start": 1000, "end": 1100, "exon_number": "1"},
    {"start": 2000, "end": 2200, "exon_number": "2"},
    {"start": 3000, "end": 3050, "exon_number": "3"},
]


def _add_ref(
    session: Session,
    gene: str = "BRCA1",
    transcript_id: str = "ENST00000000001",
    exons: list[dict] | None = None,
    strand: str = "+",
    chrom: str = "17",
) -> RefEnsemblTranscript:
    exons = exons if exons is not None else EXONS
    ref = RefEnsemblTranscript(
        gene=gene,
        ensembl_gene_id=f"ENSG_{gene}",
        ensembl_transcript_id=transcript_id,
        chrom=chrom,
        start=exons[0]["start"],
        end=exons[-1]["end"],
        strand=strand,
        exons=exons,
    )
    session.add(ref)
    session.commit()
    return ref


def _implicate(
    session: Session,
    run: AnalysisRun,
    variants: list[Variant],
    scores: list[float],
    disease_name: str = "breast cancer",
    rank_score: float = 0.9,
) -> Disease:
    disease = Disease(analysis_run_id=run.id, name=disease_name, rank_score=rank_score)
    session.add(disease)
    session.flush()
    for variant, score in zip(variants, scores, strict=True):
        # One variant can support many diseases, and VariantScore.variant_id is unique.
        existing = session.scalars(
            select(VariantScore).where(VariantScore.variant_id == variant.id)
        ).first()
        if existing is None:
            session.add(VariantScore(variant_id=variant.id, overall_score=score))
        session.add(
            DiseaseVariantLink(
                disease_id=disease.id, variant_id=variant.id, contribution_score=score
            )
        )
    session.commit()
    return disease


def _transcripts(session: Session, run: AnalysisRun) -> list[Transcript]:
    return list(
        session.scalars(
            select(Transcript)
            .where(Transcript.analysis_run_id == run.id)
            .order_by(Transcript.rank_score.desc())
        ).all()
    )


def test_no_diseases_creates_nothing(session: Session, analysis_run: AnalysisRun) -> None:
    _add_ref(session)
    assert discover_transcripts(session, analysis_run.id) == 0


def test_variant_in_exon_body_gives_exon_skip_and_knockdown(
    session: Session, analysis_run: AnalysisRun
) -> None:
    variant = make_variant(session, analysis_run, chrom="17", pos=2100, gene="BRCA1")
    _implicate(session, analysis_run, [variant], [0.9])
    _add_ref(session)

    assert discover_transcripts(session, analysis_run.id) == 1
    transcript = _transcripts(session, analysis_run)[0]

    assert [exon["exon_number"] for exon in transcript.candidate_exons] == ["2"]
    mechanisms = transcript.candidate_exons[0]["mechanisms"]
    assert "exon_skip" in mechanisms
    assert "knockdown" in mechanisms
    assert "splice_switch" not in mechanisms
    # A single-base substitution is allele-distinguishing.
    assert "allele_specific" in mechanisms


def test_variant_at_splice_boundary_gives_splice_switch(
    session: Session, analysis_run: AnalysisRun
) -> None:
    # 3bp outside the exon 2 acceptor site - inside SPLICE_WINDOW.
    variant = make_variant(session, analysis_run, chrom="17", pos=1997, gene="BRCA1")
    _implicate(session, analysis_run, [variant], [0.9])
    _add_ref(session)

    discover_transcripts(session, analysis_run.id)
    transcript = _transcripts(session, analysis_run)[0]

    mechanisms = transcript.candidate_exons[0]["mechanisms"]
    assert "splice_switch" in mechanisms
    assert "exon_skip" in mechanisms
    # Exon 2 is internal, so forcing inclusion is also on the table.
    assert "exon_include" in mechanisms


def test_terminal_exon_boundary_offers_no_exon_include(
    session: Session, analysis_run: AnalysisRun
) -> None:
    variant = make_variant(session, analysis_run, chrom="17", pos=1000, gene="BRCA1")
    _implicate(session, analysis_run, [variant], [0.9])
    _add_ref(session)

    discover_transcripts(session, analysis_run.id)
    transcript = _transcripts(session, analysis_run)[0]

    assert transcript.candidate_exons[0]["exon_number"] == "1"
    assert "exon_include" not in transcript.candidate_exons[0]["mechanisms"]


def test_knockdown_always_offered_even_without_candidate_exons(
    session: Session, analysis_run: AnalysisRun
) -> None:
    # Variant is in the gene's locus but nowhere near an exon of this isoform.
    variant = make_variant(session, analysis_run, chrom="17", pos=1500, gene="BRCA1")
    _implicate(session, analysis_run, [variant], [0.9])
    _add_ref(session)

    discover_transcripts(session, analysis_run.id)
    transcript = _transcripts(session, analysis_run)[0]

    assert transcript.candidate_exons == []
    assert transcript.opportunities == ["knockdown"]


def test_low_scoring_variant_does_not_make_a_candidate_exon(
    session: Session, analysis_run: AnalysisRun
) -> None:
    variant = make_variant(session, analysis_run, chrom="17", pos=2100, gene="BRCA1")
    _implicate(session, analysis_run, [variant], [0.05])
    _add_ref(session)

    discover_transcripts(session, analysis_run.id)
    assert _transcripts(session, analysis_run)[0].candidate_exons == []


def test_ranking_puts_strongly_implicated_transcript_first(
    session: Session, analysis_run: AnalysisRun
) -> None:
    strong_hits = [
        make_variant(session, analysis_run, chrom="17", pos=pos, gene="BRCA1")
        for pos in (1050, 2100, 3020)
    ]
    weak_hit = make_variant(session, analysis_run, chrom="17", pos=5100, gene="TP53")
    _implicate(
        session,
        analysis_run,
        [*strong_hits, weak_hit],
        [0.95, 0.9, 0.85, 0.3],
    )
    _add_ref(session, gene="BRCA1", transcript_id="ENST_STRONG")
    _add_ref(
        session,
        gene="TP53",
        transcript_id="ENST_WEAK",
        exons=[
            {"start": 5000, "end": 5200, "exon_number": "1"},
            {"start": 6000, "end": 6100, "exon_number": "2"},
        ],
    )

    assert discover_transcripts(session, analysis_run.id) == 2
    ranked = _transcripts(session, analysis_run)
    assert ranked[0].ensembl_transcript_id == "ENST_STRONG"
    assert ranked[0].rank_score > ranked[1].rank_score
    assert len(ranked[0].candidate_exons) == 3


def test_full_length_isoform_outranks_a_fragment_isoform(
    session: Session, analysis_run: AnalysisRun
) -> None:
    # Same variant, same candidate exon, so evidence/breadth/mechanisms tie exactly - only the
    # isoform-completeness tie-break separates them.
    variant = make_variant(session, analysis_run, chrom="17", pos=2100, gene="BRCA1")
    _implicate(session, analysis_run, [variant], [0.9])
    _add_ref(session, transcript_id="ENST_FRAGMENT", exons=EXONS)
    _add_ref(
        session,
        transcript_id="ENST_FULL",
        exons=[
            *EXONS,
            {"start": 4000, "end": 4100, "exon_number": "4"},
            {"start": 5000, "end": 5100, "exon_number": "5"},
        ],
    )

    discover_transcripts(session, analysis_run.id)
    ranked = _transcripts(session, analysis_run)
    assert [t.ensembl_transcript_id for t in ranked] == ["ENST_FULL", "ENST_FRAGMENT"]


def test_gene_budget_prefers_the_best_ranked_disease(
    session: Session, analysis_run: AnalysisRun
) -> None:
    top = make_variant(session, analysis_run, chrom="17", pos=2100, gene="BRCA1")
    ignored = make_variant(session, analysis_run, chrom="17", pos=5100, gene="TP53")
    _implicate(session, analysis_run, [top], [0.9], "breast cancer", rank_score=0.9)
    _implicate(session, analysis_run, [ignored], [0.9], "long qt syndrome", rank_score=0.1)
    _add_ref(session, gene="BRCA1", transcript_id="ENST_BRCA1")
    _add_ref(
        session,
        gene="TP53",
        transcript_id="ENST_TP53",
        exons=[{"start": 5000, "end": 5200, "exon_number": "1"}],
    )

    assert discover_transcripts(session, analysis_run.id, max_genes=1) == 1
    assert _transcripts(session, analysis_run)[0].gene == "BRCA1"


def test_gwas_dense_gene_does_not_starve_other_genes(
    session: Session, analysis_run: AnalysisRun
) -> None:
    """A locus carrying many trait associations must not consume the whole gene budget.

    HBB really does attract >100 separate GWAS traits, which under a disease-count cap pushed every
    other implicated gene out of transcript discovery entirely.
    """
    noisy = make_variant(session, analysis_run, chrom="11", pos=2100, gene="HBB")
    quiet = make_variant(session, analysis_run, chrom="17", pos=2100, gene="BRCA1")
    for i in range(40):
        _implicate(session, analysis_run, [noisy], [0.9], f"blood trait {i}", rank_score=0.9)
    _implicate(session, analysis_run, [quiet], [0.9], "breast cancer", rank_score=0.2)

    _add_ref(session, gene="HBB", transcript_id="ENST_HBB", chrom="11")
    _add_ref(session, gene="BRCA1", transcript_id="ENST_BRCA1", chrom="17")

    discover_transcripts(session, analysis_run.id)

    genes = {t.gene for t in _transcripts(session, analysis_run)}
    assert genes == {"HBB", "BRCA1"}


def test_ai_only_variant_still_reaches_a_design(
    session: Session, analysis_run: AnalysisRun
) -> None:
    """A gene with no validated evidence must still be designable off AI evidence alone.

    AI predictions never create a Disease, so before this these genes could not reach transcript
    discovery at all - the novel-variant arm of the pipeline could never produce a candidate.
    """
    novel = make_variant(session, analysis_run, chrom="17", pos=2100, gene="SMN2")
    add_ai_prediction(session, novel, functional_impact=0.9)
    session.add(VariantScore(variant_id=novel.id, overall_score=0.9))
    session.commit()
    _add_ref(session, gene="SMN2", transcript_id="ENST_SMN2")

    assert discover_transcripts(session, analysis_run.id) == 1
    transcript = _transcripts(session, analysis_run)[0]
    assert transcript.gene == "SMN2"
    assert transcript.evidence_basis == "ai_hypothesis"


def test_validated_genes_are_marked_validated(
    session: Session, analysis_run: AnalysisRun
) -> None:
    backed = make_variant(session, analysis_run, chrom="17", pos=2100, gene="BRCA1")
    add_annotation(session, backed, disease="breast cancer")
    _implicate(session, analysis_run, [backed], [0.9], "breast cancer", rank_score=0.9)
    _add_ref(session, gene="BRCA1", transcript_id="ENST_BRCA1")

    discover_transcripts(session, analysis_run.id)
    assert _transcripts(session, analysis_run)[0].evidence_basis == "validated"


def test_validated_genes_take_the_budget_before_ai_leads(
    session: Session, analysis_run: AnalysisRun
) -> None:
    backed = make_variant(session, analysis_run, chrom="17", pos=2100, gene="BRCA1")
    add_annotation(session, backed, disease="breast cancer")
    _implicate(session, analysis_run, [backed], [0.9], "breast cancer", rank_score=0.9)

    novel = make_variant(session, analysis_run, chrom="17", pos=2150, gene="SMN2")
    add_ai_prediction(session, novel, functional_impact=0.99)
    session.add(VariantScore(variant_id=novel.id, overall_score=0.99))
    session.commit()

    _add_ref(session, gene="BRCA1", transcript_id="ENST_BRCA1")
    _add_ref(session, gene="SMN2", transcript_id="ENST_SMN2")

    # Budget of one: validated evidence wins even though the AI variant scores higher.
    assert discover_transcripts(session, analysis_run.id, max_genes=1) == 1
    assert _transcripts(session, analysis_run)[0].gene == "BRCA1"


def test_weak_ai_only_variant_is_not_designed(
    session: Session, analysis_run: AnalysisRun
) -> None:
    novel = make_variant(session, analysis_run, chrom="17", pos=2100, gene="SMN2")
    add_ai_prediction(session, novel, functional_impact=0.05)
    session.add(VariantScore(variant_id=novel.id, overall_score=0.05))
    session.commit()
    _add_ref(session, gene="SMN2", transcript_id="ENST_SMN2")

    assert discover_transcripts(session, analysis_run.id) == 0


def test_discover_transcripts_is_idempotent(session: Session, analysis_run: AnalysisRun) -> None:
    variant = make_variant(session, analysis_run, chrom="17", pos=2100, gene="BRCA1")
    _implicate(session, analysis_run, [variant], [0.9])
    _add_ref(session)

    first = discover_transcripts(session, analysis_run.id)
    second = discover_transcripts(session, analysis_run.id)

    assert first == second == 1
    assert len(_transcripts(session, analysis_run)) == 1
