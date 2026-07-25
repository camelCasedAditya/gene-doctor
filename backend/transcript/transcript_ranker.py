"""Stage 9: discover and rank the transcripts worth targeting with an ASO.

Starting from the run's ranked diseases, walk back to the variants behind them, take those
variants' genes, and pull every annotated transcript of those genes out of the bulk-ingested
ref_ensembl_transcripts table. For each transcript we record which exons sit on top of
high-scoring variants (candidate_exons) and which ASO mechanisms that geometry makes plausible
(opportunities).

Everything here is a research heuristic: exon/variant geometry is a proxy for "this transcript
looks addressable", not evidence that an ASO against it would work.
"""

from __future__ import annotations

from sqlalchemy import delete, select
from sqlalchemy.orm import Session

from backend.database.models import (
    Aso,
    Disease,
    DiseaseVariantLink,
    RefEnsemblTranscript,
    Transcript,
    Variant,
    VariantScore,
)

# Canonical mechanism vocabulary, in the order we report it.
MECHANISM_ORDER = ("splice_switch", "exon_skip", "exon_include", "knockdown", "allele_specific")

# A variant this close to an exon edge is treated as splice-site adjacent. ~10bp covers the
# canonical donor/acceptor dinucleotides plus the immediately flanking consensus positions.
SPLICE_WINDOW = 10

# Variants below this overall_score don't make an exon a candidate - they are in the disease's
# evidence chain but too weak to steer a design.
MIN_VARIANT_SCORE = 0.25

# How many distinct genes to carry into transcript discovery.
MAX_TARGET_GENES = 25

# Saturation points for the ranking terms (see _rank_score).
EVIDENCE_SATURATION = 3.0
EXON_SATURATION = 3.0
EXON_COUNT_SATURATION = 30.0


def discover_transcripts(
    session: Session, analysis_run_id: str, max_genes: int = MAX_TARGET_GENES
) -> int:
    """For genes behind the top-ranked diseases, create ranked Transcript rows. Returns count created.

    Diseases are consumed in rank order until `max_genes` distinct genes are collected, rather than
    truncating to a fixed number of diseases. A gene in a GWAS-dense locus (HBB, say) can carry
    over a hundred separate trait associations, which under a disease-count cap would fill every
    slot and starve every other implicated gene of a design.
    """
    _clear_run(session, analysis_run_id)

    ranked_disease_ids = session.scalars(
        select(Disease.id)
        .where(Disease.analysis_run_id == analysis_run_id)
        .order_by(Disease.rank_score.desc())
    ).all()

    by_gene: dict[str, list[tuple[Variant, float]]] = {}
    for disease_id in ranked_disease_ids:
        for variant, score in _implicated_variants(session, analysis_run_id, [disease_id]):
            if not variant.gene:
                continue
            if variant.gene not in by_gene and len(by_gene) >= max_genes:
                continue
            by_gene.setdefault(variant.gene, []).append((variant, score))

    # Validated genes claim the budget first; whatever remains goes to AI-only leads.
    validated_genes = set(by_gene)
    for variant, score in _ai_only_variants(session, analysis_run_id):
        if variant.gene in validated_genes:
            continue  # already covered by validated evidence
        if variant.gene not in by_gene and len(by_gene) >= max_genes:
            continue
        by_gene.setdefault(variant.gene, []).append((variant, score))

    if not by_gene:
        return 0

    refs = session.scalars(
        select(RefEnsemblTranscript).where(RefEnsemblTranscript.gene.in_(by_gene))
    ).all()

    created = 0
    for ref in refs:
        gene_variants = [
            (variant, score)
            for variant, score in by_gene[ref.gene]
            if variant.chrom == ref.chrom
        ]
        candidate_exons = _candidate_exons(ref.exons, gene_variants)
        opportunities = _opportunities(candidate_exons)
        session.add(
            Transcript(
                analysis_run_id=analysis_run_id,
                gene=ref.gene,
                ensembl_transcript_id=ref.ensembl_transcript_id,
                chrom=ref.chrom,
                strand=ref.strand,
                exon_structure=ref.exons,
                candidate_exons=candidate_exons,
                opportunities=opportunities,
                rank_score=_rank_score(candidate_exons, opportunities, ref.exons),
                evidence_basis="validated" if ref.gene in validated_genes else "ai_hypothesis",
            )
        )
        created += 1

    session.commit()
    return created


def _clear_run(session: Session, analysis_run_id: str) -> None:
    """Make the stage idempotent: Aso rows FK to Transcript, so they go first."""
    run_transcripts = select(Transcript.id).where(Transcript.analysis_run_id == analysis_run_id)
    session.execute(delete(Aso).where(Aso.transcript_id.in_(run_transcripts)))
    session.execute(delete(Transcript).where(Transcript.analysis_run_id == analysis_run_id))
    session.commit()


def _ai_only_variants(session: Session, analysis_run_id: str) -> list[tuple[Variant, float]]:
    """High-scoring variants of this run whose only evidence is an AI prediction.

    These never form a Disease - AI output is not allowed to invent a disease - so without this
    they could never reach a design at all, and the whole novel-variant arm of the pipeline would
    be unable to produce a candidate. Genes reached this way are marked `ai_hypothesis` so the
    weaker basis stays visible everywhere downstream.
    """
    rows = session.execute(
        select(Variant, VariantScore.overall_score)
        .join(VariantScore, VariantScore.variant_id == Variant.id)
        .where(
            Variant.analysis_run_id == analysis_run_id,
            Variant.gene.isnot(None),
            VariantScore.overall_score >= MIN_VARIANT_SCORE,
            ~Variant.annotations.any(),
            Variant.ai_predictions.any(),
        )
        .order_by(VariantScore.overall_score.desc())
    ).all()
    return [(variant, score or 0.0) for variant, score in rows]


def _implicated_variants(
    session: Session, analysis_run_id: str, disease_ids: list[str]
) -> list[tuple[Variant, float]]:
    """This run's variants linked to the given diseases, paired with their overall_score."""
    rows = session.execute(
        select(Variant, VariantScore.overall_score)
        .join(DiseaseVariantLink, DiseaseVariantLink.variant_id == Variant.id)
        .outerjoin(VariantScore, VariantScore.variant_id == Variant.id)
        .where(
            DiseaseVariantLink.disease_id.in_(disease_ids),
            Variant.analysis_run_id == analysis_run_id,
        )
    ).all()
    # A variant can back several of the top diseases; keep one entry each.
    deduped = {variant.id: (variant, score or 0.0) for variant, score in rows}
    return list(deduped.values())


def _candidate_exons(
    exons: list[dict], gene_variants: list[tuple[Variant, float]]
) -> list[dict]:
    """Exons overlapping (or splice-adjacent to) a high-scoring variant, with their mechanisms."""
    candidates: list[dict] = []
    for index, exon in enumerate(exons):
        start, end = int(exon["start"]), int(exon["end"])
        hits = [
            (variant, score)
            for variant, score in gene_variants
            if score >= MIN_VARIANT_SCORE
            and start - SPLICE_WINDOW <= variant.pos <= end + SPLICE_WINDOW
        ]
        if not hits:
            continue

        internal = 0 < index < len(exons) - 1
        mechanisms: set[str] = set()
        for variant, _ in hits:
            mechanisms |= _variant_mechanisms(variant, start, end, internal)

        candidates.append({
            "exon_number": str(exon.get("exon_number", index + 1)),
            "start": start,
            "end": end,
            "variant_score": round(sum(score for _, score in hits), 4),
            "variant_positions": [variant.pos for variant, _ in hits],
            "variants": [
                f"{variant.chrom}:{variant.pos}{variant.ref}>{variant.alt}" for variant, _ in hits
            ],
            "mechanisms": _ordered(mechanisms),
        })
    return candidates


def _variant_mechanisms(variant: Variant, start: int, end: int, internal: bool) -> set[str]:
    """Heuristic mechanism rules for one variant against one exon.

    - within SPLICE_WINDOW of an exon edge -> the splice decision itself is addressable:
      splice_switch, and exon_skip because redirecting that decision drops the exon. On an
      internal exon the reverse is also plausible (blocking a silencer to force inclusion),
      so exon_include is offered too - first/last exons cannot be skipped or re-included.
    - inside the exon body -> exon_skip (block an exonic splicing enhancer) and knockdown
      (RNase-H cleavage of the mature transcript).
    - a single-nucleotide substitution gives a distinguishing allele, so an allele-selective
      ASO is conceivable. The schema carries no genotype, so this is a variant-class proxy
      for "heterozygous-style", not a real zygosity call.
    """
    mechanisms: set[str] = set()
    if min(abs(variant.pos - start), abs(variant.pos - end)) <= SPLICE_WINDOW:
        mechanisms |= {"splice_switch", "exon_skip"}
        if internal:
            mechanisms.add("exon_include")
    else:
        mechanisms |= {"exon_skip", "knockdown"}
    if len(variant.ref) == 1 and len(variant.alt) == 1 and variant.ref != variant.alt:
        mechanisms.add("allele_specific")
    return mechanisms


def _opportunities(candidate_exons: list[dict]) -> list[str]:
    """Union of the candidate exons' mechanisms. knockdown is always on the table: any expressed
    transcript can in principle be cleaved, regardless of where the variants sit.
    """
    mechanisms = {"knockdown"}
    for exon in candidate_exons:
        mechanisms |= set(exon["mechanisms"])
    return _ordered(mechanisms)


def _rank_score(
    candidate_exons: list[dict], opportunities: list[str], exon_structure: list[dict]
) -> float:
    """How strongly this transcript is implicated, in [0, 1]:

        0.60 * evidence + 0.25 * breadth + 0.13 * mechanisms + 0.02 * completeness

        evidence     = min(sum of overlapping variant scores / 3.0, 1.0)  how pathogenic the hits are
        breadth      = min(number of candidate exons / 3.0, 1.0)          how much of it is implicated
        mechanisms   = len(opportunities) / 5                             how many ways in there are
        completeness = min(annotated exon count / 30.0, 1.0)              isoform tie-break only

    Evidence dominates. The completeness term carries almost no weight on purpose: the same
    variants land on every isoform of a gene, which ties their scores exactly, and without a
    tie-break the arbitrary winner is often a two-exon fragment isoform rather than the
    full-length transcript. Weights are judgement calls, not fitted to any outcome data.
    """
    evidence = min(sum(exon["variant_score"] for exon in candidate_exons) / EVIDENCE_SATURATION, 1.0)
    breadth = min(len(candidate_exons) / EXON_SATURATION, 1.0)
    mechanisms = len(opportunities) / len(MECHANISM_ORDER)
    completeness = min(len(exon_structure) / EXON_COUNT_SATURATION, 1.0)
    return round(0.60 * evidence + 0.25 * breadth + 0.13 * mechanisms + 0.02 * completeness, 4)


def _ordered(mechanisms: set[str]) -> list[str]:
    return [mechanism for mechanism in MECHANISM_ORDER if mechanism in mechanisms]
