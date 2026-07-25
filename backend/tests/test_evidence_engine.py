"""Weighted evidence engine: hand-computed scores, renormalization over missing sources,
weight configurability and idempotent re-scoring.
"""

from __future__ import annotations

import pytest
from sqlalchemy import func, select
from sqlalchemy.orm import Session

from backend.database.models import AnalysisRun, EvidenceWeights, VariantScore
from backend.services.evidence_engine import get_weights, score_variants
from backend.tests.conftest import add_ai_prediction, add_annotation, make_variant


def score_of(session: Session, variant_id: str) -> float:
    return session.scalars(
        select(VariantScore.overall_score).where(VariantScore.variant_id == variant_id)
    ).one()


def test_get_weights_creates_defaults(session: Session) -> None:
    weights = get_weights(session)
    assert (weights.clinvar_weight, weights.gwas_weight, weights.ai_weight) == (0.45, 0.30, 0.25)
    get_weights(session)
    assert session.scalar(select(func.count()).select_from(EvidenceWeights)) == 1


def test_all_three_sources(session: Session, analysis_run: AnalysisRun) -> None:
    variant = make_variant(session, analysis_run)
    add_annotation(session, variant, source="clinvar", evidence_strength=1.0)
    add_annotation(session, variant, source="gwas", evidence_strength=0.6)
    add_ai_prediction(
        session,
        variant,
        functional_impact=0.8,
        regulatory_impact=0.5,
        pathogenicity_probability=0.7,
    )

    assert score_variants(session, analysis_run.id) == 1
    # 0.45*1.0 + 0.30*0.6 + 0.25*0.8 (ai = max of the three signals), weights sum to 1.0
    assert score_of(session, variant.id) == pytest.approx(0.83)


def test_clinvar_only_renormalizes(session: Session, analysis_run: AnalysisRun) -> None:
    variant = make_variant(session, analysis_run)
    add_annotation(session, variant, source="clinvar", evidence_strength=0.8)

    score_variants(session, analysis_run.id)
    # 0.45*0.8 / 0.45 == 0.8, not 0.36
    assert score_of(session, variant.id) == pytest.approx(0.8)


def test_gwas_only_renormalizes(session: Session, analysis_run: AnalysisRun) -> None:
    variant = make_variant(session, analysis_run)
    add_annotation(session, variant, source="gwas", evidence_strength=0.8)

    score_variants(session, analysis_run.id)
    # 0.30*0.8 / 0.30 == 0.8, not 0.24
    assert score_of(session, variant.id) == pytest.approx(0.8)


def test_ai_only_renormalizes(session: Session, analysis_run: AnalysisRun) -> None:
    variant = make_variant(session, analysis_run)
    add_ai_prediction(
        session,
        variant,
        functional_impact=None,
        regulatory_impact=0.9,
        pathogenicity_probability=None,
    )

    score_variants(session, analysis_run.id)
    # None signals are absent, not zero: ai = 0.9, and 0.25*0.9 / 0.25 == 0.9
    assert score_of(session, variant.id) == pytest.approx(0.9)


def test_two_of_three_sources(session: Session, analysis_run: AnalysisRun) -> None:
    variant = make_variant(session, analysis_run)
    add_annotation(session, variant, source="clinvar", evidence_strength=1.0)
    add_annotation(session, variant, source="gwas", evidence_strength=0.5)

    score_variants(session, analysis_run.id)
    # (0.45*1.0 + 0.30*0.5) / 0.75 == 0.8
    assert score_of(session, variant.id) == pytest.approx(0.8)


def test_max_within_a_source(session: Session, analysis_run: AnalysisRun) -> None:
    variant = make_variant(session, analysis_run)
    add_annotation(session, variant, source="clinvar", evidence_strength=0.2)
    add_annotation(session, variant, source="clinvar", evidence_strength=0.95)
    add_ai_prediction(session, variant, source="variant_lm", functional_impact=0.3,
                      regulatory_impact=0.3, pathogenicity_probability=0.3)
    add_ai_prediction(session, variant, source="alphamissense", functional_impact=0.6,
                      regulatory_impact=None, pathogenicity_probability=None)

    score_variants(session, analysis_run.id)
    # (0.45*0.95 + 0.25*0.6) / 0.70
    assert score_of(session, variant.id) == pytest.approx((0.45 * 0.95 + 0.25 * 0.6) / 0.70)


def test_no_evidence_scores_zero(session: Session, analysis_run: AnalysisRun) -> None:
    variant = make_variant(session, analysis_run)
    assert score_variants(session, analysis_run.id) == 1
    assert score_of(session, variant.id) == 0.0


def test_null_evidence_strength_is_absent(session: Session, analysis_run: AnalysisRun) -> None:
    variant = make_variant(session, analysis_run)
    add_annotation(session, variant, source="clinvar", evidence_strength=None)
    add_annotation(session, variant, source="gwas", evidence_strength=0.4)

    score_variants(session, analysis_run.id)
    assert score_of(session, variant.id) == pytest.approx(0.4)


def test_weights_are_configurable(session: Session, analysis_run: AnalysisRun) -> None:
    variant = make_variant(session, analysis_run)
    add_annotation(session, variant, source="clinvar", evidence_strength=1.0)
    add_annotation(session, variant, source="gwas", evidence_strength=0.5)

    score_variants(session, analysis_run.id)
    assert score_of(session, variant.id) == pytest.approx(0.8)

    weights = get_weights(session)
    weights.clinvar_weight = 0.1
    weights.gwas_weight = 0.9
    session.commit()

    score_variants(session, analysis_run.id)
    # (0.1*1.0 + 0.9*0.5) / 1.0
    assert score_of(session, variant.id) == pytest.approx(0.55)


def test_score_variants_is_idempotent(session: Session, analysis_run: AnalysisRun) -> None:
    variant = make_variant(session, analysis_run)
    add_annotation(session, variant, source="clinvar", evidence_strength=0.8)

    assert score_variants(session, analysis_run.id) == 1
    assert score_variants(session, analysis_run.id) == 1
    assert session.scalar(select(func.count()).select_from(VariantScore)) == 1


def test_empty_run_writes_nothing(session: Session, analysis_run: AnalysisRun) -> None:
    assert score_variants(session, analysis_run.id) == 0
    assert session.scalar(select(func.count()).select_from(VariantScore)) == 0
