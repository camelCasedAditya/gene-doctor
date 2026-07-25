"""Rank candidate diseases for an analysis run by aggregating its scored variants.

Validated and predicted contributors stay separated end to end: a disease exists because
Annotation rows (ClinVar/GWAS) point at it, and AI-only variants can only ever join as *predicted*
contributors, counted in their own column so the UI can show the provenance split.
"""

from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Iterable

from sqlalchemy import delete, select
from sqlalchemy.orm import Session

from backend.annotation.annotate import normalize_disease
from backend.database.models import (
    AiPrediction,
    Annotation,
    Disease,
    DiseaseVariantLink,
    Variant,
    VariantScore,
)


def rank_diseases(session: Session, analysis_run_id: str) -> int:
    """Aggregate scored variants by disease into Disease + DiseaseVariantLink rows. Returns count
    of diseases ranked.
    """
    _clear_run(session, analysis_run_id)

    variants = {
        variant.id: variant
        for variant in session.scalars(
            select(Variant).where(Variant.analysis_run_id == analysis_run_id)
        )
    }
    if not variants:
        session.commit()
        return 0

    variant_ids = list(variants)
    scores = dict(
        session.execute(
            select(VariantScore.variant_id, VariantScore.overall_score).where(
                VariantScore.variant_id.in_(variant_ids)
            )
        ).all()
    )

    # Validated evidence defines which diseases are on the table at all.
    known: dict[str, set[str]] = defaultdict(set)
    annotated: set[str] = set()
    rows = session.execute(
        select(Annotation.variant_id, Annotation.disease).where(
            Annotation.variant_id.in_(variant_ids)
        )
    )
    for variant_id, disease in rows:
        annotated.add(variant_id)
        if disease:
            known[normalize_disease(disease)].add(variant_id)

    ai_only_by_gene = _ai_only_by_gene(session, variants, annotated)

    for name, known_ids in known.items():
        genes = {variants[v].gene for v in known_ids if variants[v].gene}
        predicted_ids = {v for gene in genes for v in ai_only_by_gene.get(gene, ())}
        disease = Disease(
            analysis_run_id=analysis_run_id,
            name=name,
            known_variant_count=len(known_ids),
            predicted_variant_count=len(predicted_ids),
        )
        contributions = {v: scores.get(v, 0.0) for v in known_ids | predicted_ids}
        disease.rank_score = _rank_score(contributions.values())
        disease.confidence = _confidence(contributions.values(), len(known_ids))
        session.add(disease)
        session.flush()
        for variant_id, contribution in contributions.items():
            session.add(
                DiseaseVariantLink(
                    disease_id=disease.id,
                    variant_id=variant_id,
                    contribution_score=contribution,
                )
            )

    session.commit()
    return len(known)


def _rank_score(contributions: Iterable[float]) -> float:
    """mean(variant scores) * sqrt(number of supporting variants).

    Equivalently sum / sqrt(n): breadth of support raises the ranking with diminishing returns, so
    a well-supported disease cannot be outranked by a single weak variant (five variants at 0.4
    score 0.89, one variant at 0.2 scores 0.2), while a pile of near-zero variants still cannot
    manufacture a high rank. Gene disruption and regulatory impact already enter through the AI
    component of each VariantScore, so they are not re-applied here.
    """
    values = list(contributions)
    if not values:
        return 0.0
    return sum(values) / math.sqrt(len(values))


def _confidence(contributions: Iterable[float], known_count: int) -> float:
    """Mean variant score, halved as the evidence base shifts from validated to AI-predicted.

    A disease supported entirely by ClinVar/GWAS variants keeps its full mean; one carried mostly
    by AI hypotheses is discounted toward 0.5x. The known/predicted counts on the Disease row
    remain the authoritative split - this is only a summary for sorting.
    """
    values = list(contributions)
    if not values:
        return 0.0
    validated_fraction = known_count / len(values)
    return min(1.0, (sum(values) / len(values)) * (0.5 + 0.5 * validated_fraction))


def _ai_only_by_gene(
    session: Session, variants: dict[str, Variant], annotated: set[str]
) -> dict[str, set[str]]:
    """Variants whose ONLY evidence is an AiPrediction, indexed by gene.

    These have no disease of their own (AiPrediction carries no disease field), so they attach as
    predicted contributors to the diseases that validated variants in the same gene point at.

    ponytail: gene-level attachment only; upgrade to transcript/pathway-aware propagation if the
    hypotheses get too broad.
    """
    by_gene: dict[str, set[str]] = defaultdict(set)
    predicted_ids = session.scalars(
        select(AiPrediction.variant_id)
        .where(AiPrediction.variant_id.in_(list(variants)))
        .distinct()
    )
    for variant_id in predicted_ids:
        gene = variants[variant_id].gene
        if gene and variant_id not in annotated:
            by_gene[gene].add(variant_id)
    return by_gene


def _clear_run(session: Session, analysis_run_id: str) -> None:
    """Drop the run's previous ranking so re-ranking is idempotent."""
    prior = session.scalars(
        select(Disease.id).where(Disease.analysis_run_id == analysis_run_id)
    ).all()
    if not prior:
        return
    session.execute(delete(DiseaseVariantLink).where(DiseaseVariantLink.disease_id.in_(prior)))
    session.execute(delete(Disease).where(Disease.id.in_(prior)))
    session.flush()
