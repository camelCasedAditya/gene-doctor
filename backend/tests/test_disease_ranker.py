"""Disease ranking: support-weighted ordering, the validated/predicted provenance split and
idempotent re-ranking.
"""

from __future__ import annotations

import math

import pytest
from sqlalchemy import func, select
from sqlalchemy.orm import Session

from backend.database.models import AnalysisRun, Disease, DiseaseVariantLink
from backend.services.disease_ranker import rank_diseases
from backend.services.evidence_engine import score_variants
from backend.tests.conftest import add_ai_prediction, add_annotation, make_variant


def disease_by_name(session: Session, run: AnalysisRun, name: str) -> Disease:
    return session.scalars(
        select(Disease).where(Disease.analysis_run_id == run.id, Disease.name == name)
    ).one()


def test_well_supported_disease_outranks_single_weak_variant(
    session: Session, analysis_run: AnalysisRun
) -> None:
    for pos in (1000, 2000, 3000):
        variant = make_variant(session, analysis_run, pos=pos, gene="APOE")
        add_annotation(session, variant, disease="alzheimer disease", evidence_strength=1.0)
    weak = make_variant(session, analysis_run, pos=9000, gene="XYZ1")
    add_annotation(session, weak, disease="rare thing", evidence_strength=0.2)

    score_variants(session, analysis_run.id)
    assert rank_diseases(session, analysis_run.id) == 2

    supported = disease_by_name(session, analysis_run, "alzheimer disease")
    lonely = disease_by_name(session, analysis_run, "rare thing")
    # three variants at 1.0 -> 3 / sqrt(3); one variant at 0.2 -> 0.2
    assert supported.rank_score == pytest.approx(3.0 / math.sqrt(3))
    assert lonely.rank_score == pytest.approx(0.2)
    assert supported.rank_score > lonely.rank_score


def test_many_weak_variants_do_not_beat_strong_support(
    session: Session, analysis_run: AnalysisRun
) -> None:
    for pos in range(1000, 1010):
        variant = make_variant(session, analysis_run, pos=pos, gene="NOISE1")
        add_annotation(session, variant, disease="noisy trait", evidence_strength=0.05)
    for pos in range(5000, 5003):
        variant = make_variant(session, analysis_run, pos=pos, gene="REAL1")
        add_annotation(session, variant, disease="real disease", evidence_strength=0.95)

    score_variants(session, analysis_run.id)
    rank_diseases(session, analysis_run.id)

    assert (
        disease_by_name(session, analysis_run, "real disease").rank_score
        > disease_by_name(session, analysis_run, "noisy trait").rank_score
    )


def test_known_predicted_split(session: Session, analysis_run: AnalysisRun) -> None:
    validated = make_variant(session, analysis_run, pos=1000, gene="BRCA1")
    add_annotation(session, validated, disease="breast cancer", evidence_strength=0.9)
    ai_same_gene = make_variant(session, analysis_run, pos=1500, gene="BRCA1")
    add_ai_prediction(session, ai_same_gene)  # 0.8 / 0.5 / 0.7 -> ai signal 0.8
    ai_other_gene = make_variant(session, analysis_run, pos=8000, gene="TP53")
    add_ai_prediction(session, ai_other_gene)

    score_variants(session, analysis_run.id)
    assert rank_diseases(session, analysis_run.id) == 1

    disease = disease_by_name(session, analysis_run, "breast cancer")
    assert disease.known_variant_count == 1
    assert disease.predicted_variant_count == 1
    links = session.scalars(
        select(DiseaseVariantLink).where(DiseaseVariantLink.disease_id == disease.id)
    ).all()
    assert {link.variant_id for link in links} == {validated.id, ai_same_gene.id}
    # (0.9 + 0.8) / sqrt(2), confidence = mean 0.85 discounted by half-AI support
    assert disease.rank_score == pytest.approx(1.7 / math.sqrt(2))
    assert disease.confidence == pytest.approx(0.85 * 0.75)
    assert 0.0 <= disease.confidence <= 1.0


def test_variant_with_both_evidence_types_counts_as_known(
    session: Session, analysis_run: AnalysisRun
) -> None:
    variant = make_variant(session, analysis_run, gene="SMN1")
    add_annotation(session, variant, disease="spinal muscular atrophy", evidence_strength=1.0)
    add_ai_prediction(session, variant)

    score_variants(session, analysis_run.id)
    rank_diseases(session, analysis_run.id)

    disease = disease_by_name(session, analysis_run, "spinal muscular atrophy")
    assert (disease.known_variant_count, disease.predicted_variant_count) == (1, 0)


def test_fully_validated_disease_keeps_full_confidence(
    session: Session, analysis_run: AnalysisRun
) -> None:
    variant = make_variant(session, analysis_run, gene="CFTR")
    add_annotation(session, variant, disease="cystic fibrosis", evidence_strength=1.0)

    score_variants(session, analysis_run.id)
    rank_diseases(session, analysis_run.id)

    disease = disease_by_name(session, analysis_run, "cystic fibrosis")
    assert disease.confidence == pytest.approx(1.0)


def test_rank_diseases_is_idempotent(session: Session, analysis_run: AnalysisRun) -> None:
    variant = make_variant(session, analysis_run, gene="APOE")
    add_annotation(session, variant, disease="alzheimer disease", evidence_strength=1.0)
    score_variants(session, analysis_run.id)

    assert rank_diseases(session, analysis_run.id) == 1
    first_score = disease_by_name(session, analysis_run, "alzheimer disease").rank_score
    assert rank_diseases(session, analysis_run.id) == 1

    assert session.scalar(select(func.count()).select_from(Disease)) == 1
    assert session.scalar(select(func.count()).select_from(DiseaseVariantLink)) == 1
    assert disease_by_name(session, analysis_run, "alzheimer disease").rank_score == first_score


def test_unannotated_run_ranks_nothing(session: Session, analysis_run: AnalysisRun) -> None:
    variant = make_variant(session, analysis_run)
    add_ai_prediction(session, variant)
    score_variants(session, analysis_run.id)

    assert rank_diseases(session, analysis_run.id) == 0
    assert session.scalar(select(func.count()).select_from(Disease)) == 0
