"""Weighted evidence aggregation: collapse a variant's evidence into one comparable score.

Three sources are scored independently - ClinVar and GWAS Catalog annotations (validated) and
AiPrediction rows (AI hypotheses) - then combined with the configurable weights held in the
single EvidenceWeights row, so the user can re-weight via the API and re-score.
"""

from __future__ import annotations

from sqlalchemy import select
from sqlalchemy.orm import Session

from backend.database.models import (
    AiPrediction,
    Annotation,
    EvidenceWeights,
    Variant,
    VariantScore,
)


def get_weights(session: Session) -> EvidenceWeights:
    """Return the single EvidenceWeights row, creating it with defaults if absent."""
    weights = session.get(EvidenceWeights, 1)
    if weights is None:
        weights = EvidenceWeights(id=1)
        session.add(weights)
        session.commit()
    return weights


def score_variants(session: Session, analysis_run_id: str) -> int:
    """Compute a VariantScore row for every variant in the run. Returns count of scores written."""
    weights = get_weights(session)
    variant_ids = session.scalars(
        select(Variant.id).where(Variant.analysis_run_id == analysis_run_id)
    ).all()
    if not variant_ids:
        return 0

    # ponytail: one in_() batch per query; chunk the ids if a run ever exceeds SQLite's
    # bound-parameter limit (~32k).
    clinvar, gwas = _annotation_scores(session, variant_ids)
    ai = _ai_scores(session, variant_ids)
    existing = {
        score.variant_id: score
        for score in session.scalars(
            select(VariantScore).where(VariantScore.variant_id.in_(variant_ids))
        )
    }

    for variant_id in variant_ids:
        overall = _combine(
            weights,
            clinvar.get(variant_id),
            gwas.get(variant_id),
            ai.get(variant_id),
        )
        score = existing.get(variant_id)
        if score is None:
            # variant_id is UNIQUE, so updating in place keeps re-scoring idempotent.
            session.add(VariantScore(variant_id=variant_id, overall_score=overall))
        else:
            score.overall_score = overall

    session.commit()
    return len(variant_ids)


def _combine(
    weights: EvidenceWeights,
    clinvar_score: float | None,
    gwas_score: float | None,
    ai_score: float | None,
) -> float:
    """Weighted mean over the sources that actually have data.

    Renormalization matters: a variant missing from ClinVar has not been called benign there, it
    simply was not measured. Zero-filling that source would punish it for a gap in the reference
    database, so the weighted sum is divided by the weights of the present sources only. A variant
    whose sole evidence is a GWAS signal of 0.8 therefore scores 0.8, not 0.30 * 0.8 = 0.24, and
    stays comparable to a ClinVar-backed variant. No evidence at all -> 0.0.
    """
    present = [
        (weight, score)
        for weight, score in (
            (weights.clinvar_weight, clinvar_score),
            (weights.gwas_weight, gwas_score),
            (weights.ai_weight, ai_score),
        )
        if score is not None
    ]
    total_weight = sum(weight for weight, _ in present)
    if not total_weight:
        return 0.0
    return sum(weight * score for weight, score in present) / total_weight


def _annotation_scores(
    session: Session, variant_ids: list[str]
) -> tuple[dict[str, float], dict[str, float]]:
    """Per variant, the strongest ClinVar and GWAS evidence_strength. A row with a NULL strength
    carries no signal and is treated as absent; a 0.0 strength (e.g. a ClinVar VUS) is real data
    and does count as that source being present.
    """
    clinvar: dict[str, float] = {}
    gwas: dict[str, float] = {}
    rows = session.execute(
        select(Annotation.variant_id, Annotation.source, Annotation.evidence_strength).where(
            Annotation.variant_id.in_(variant_ids)
        )
    )
    for variant_id, source, strength in rows:
        if strength is None:
            continue
        bucket = {"clinvar": clinvar, "gwas": gwas}.get(source)
        if bucket is None:
            continue
        bucket[variant_id] = max(bucket.get(variant_id, strength), strength)
    return clinvar, gwas


def _ai_scores(session: Session, variant_ids: list[str]) -> dict[str, float]:
    """Per variant, the max AI signal across every AI source and signal type.

    Max rather than a mean: functional impact, regulatory impact and pathogenicity probability
    describe different mechanisms, so a strong hit on any one of them is a real hypothesis and
    averaging it against the mechanisms that do not apply would wash it out. None fields are
    absent signals, not zeros.
    """
    ai: dict[str, float] = {}
    rows = session.execute(
        select(
            AiPrediction.variant_id,
            AiPrediction.functional_impact,
            AiPrediction.regulatory_impact,
            AiPrediction.pathogenicity_probability,
        ).where(AiPrediction.variant_id.in_(variant_ids))
    )
    for variant_id, *signals in rows:
        present = [signal for signal in signals if signal is not None]
        if present:
            best = max(present)
            ai[variant_id] = max(ai.get(variant_id, best), best)
    return ai
