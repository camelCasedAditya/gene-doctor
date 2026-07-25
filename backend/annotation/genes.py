"""Assign a gene and transcript to variants from genomic coordinates.

This runs before annotation, because gene membership is a property of where a variant sits, not of
whether anyone has published evidence about it. Deriving the gene from ClinVar/GWAS records alone
leaves every novel variant with no gene at all, which in turn makes it invisible to transcript
discovery and ASO design - the entire novel-variant arm of the pipeline.
"""

from __future__ import annotations

from sqlalchemy import select
from sqlalchemy.orm import Session

from backend.database.models import RefEnsemblTranscript, Variant


def assign_genes(session: Session, analysis_run_id: str) -> int:
    """Fill Variant.gene / Variant.transcript_id by overlap with Ensembl transcripts.

    Returns the number of variants that got a gene.
    """
    variants = session.scalars(
        select(Variant)
        .where(Variant.analysis_run_id == analysis_run_id)
        .order_by(Variant.chrom, Variant.pos)
    ).all()
    if not variants:
        return 0

    # ponytail: per-chromosome cache plus a linear scan. Fine for the thousands of variants a
    # locus-scale run produces; swap in an interval tree if whole-genome variant counts land here.
    cache: dict[str, list[RefEnsemblTranscript]] = {}
    assigned = 0

    for variant in variants:
        if variant.chrom not in cache:
            cache[variant.chrom] = list(
                session.scalars(
                    select(RefEnsemblTranscript).where(
                        RefEnsemblTranscript.chrom == variant.chrom
                    )
                ).all()
            )

        overlapping = [
            transcript
            for transcript in cache[variant.chrom]
            if transcript.start <= variant.pos <= transcript.end
        ]
        if not overlapping:
            continue

        best = _best_transcript(overlapping, variant.pos)
        variant.gene = best.gene
        variant.transcript_id = best.ensembl_transcript_id
        assigned += 1

    session.commit()
    return assigned


def _best_transcript(candidates: list[RefEnsemblTranscript], pos: int) -> RefEnsemblTranscript:
    """Prefer a transcript whose exon actually contains the position, then the most complete one.

    A variant inside an exon is the interesting case for ASO design, so a transcript that places it
    exonically beats one that leaves it intronic; exon count breaks the remaining ties so a full
    isoform wins over a two-exon fragment.
    """
    return max(
        candidates,
        key=lambda t: (
            any(exon["start"] <= pos <= exon["end"] for exon in t.exons),
            len(t.exons),
        ),
    )
