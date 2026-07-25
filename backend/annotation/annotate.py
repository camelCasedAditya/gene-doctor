"""Annotate called variants with validated clinical evidence from ClinVar and the GWAS Catalog.

ClinVar is queried directly from its tabix-indexed VCF; GWAS associations come from the
bulk-ingested ref_gwas_associations table. Both produce Annotation rows - validated evidence only,
never mixed with AiPrediction rows.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import pysam
from sqlalchemy import select
from sqlalchemy.orm import Session

from backend.database.models import Annotation, RefGwasAssociation, Variant
from backend.settings import settings

# ClinVar CLNSIG values that count as meaningful pathogenic/benign evidence, mapped to an
# evidence_strength in [0, 1]. Anything not listed (VUS, conflicting, not_provided) scores 0.0
# and therefore leaves the variant to be picked up by the AI predictor instead.
CLNSIG_STRENGTH = {
    "Pathogenic": 1.0,
    "Pathogenic/Likely_pathogenic": 0.95,
    "Likely_pathogenic": 0.8,
    "Benign": 0.1,
    "Benign/Likely_benign": 0.1,
    "Likely_benign": 0.15,
}

GWAS_WINDOW = 100  # bp; GWAS tag SNPs rarely sit exactly on a called variant


@dataclass
class AnnotationCounts:
    clinvar: int = 0
    gwas: int = 0
    diseases: set[str] = field(default_factory=set)


def normalize_disease(name: str) -> str:
    """Normalized disease key used to aggregate the same disease across ClinVar and GWAS.

    ponytail: string normalization only, upgrade to EFO/MONDO ontology mapping if cross-source
    aggregation gets noisy.
    """
    return " ".join(name.replace("_", " ").lower().split())


def _clinvar_strength(clnsig: str) -> float:
    # Longest prefix wins: "Pathogenic/Likely_pathogenic" also starts with "Pathogenic", and must
    # not be scored as if it were an unqualified Pathogenic call.
    for key in sorted(CLNSIG_STRENGTH, key=len, reverse=True):
        if clnsig.startswith(key):
            return CLNSIG_STRENGTH[key]
    return 0.0


def annotate_variants(session: Session, analysis_run_id: str) -> AnnotationCounts:
    """Annotate every variant of a run against ClinVar + GWAS. Returns counts of evidence found."""
    counts = AnnotationCounts()
    variants = session.scalars(
        select(Variant).where(Variant.analysis_run_id == analysis_run_id)
    ).all()
    if not variants:
        return counts

    clinvar = pysam.VariantFile(str(settings.clinvar_vcf))

    for variant in variants:
        for annotation in _clinvar_annotations(clinvar, variant):
            session.add(annotation)
            counts.clinvar += 1
            if annotation.disease:
                counts.diseases.add(annotation.disease)

        for annotation in _gwas_annotations(session, variant):
            session.add(annotation)
            counts.gwas += 1
            if annotation.disease:
                counts.diseases.add(annotation.disease)

    clinvar.close()
    session.commit()
    return counts


def _clinvar_annotations(clinvar: pysam.VariantFile, variant: Variant) -> list[Annotation]:
    """Exact position+allele matches in ClinVar for one variant."""
    annotations: list[Annotation] = []
    try:
        records = clinvar.fetch(variant.chrom, variant.pos - 1, variant.pos)
    except (ValueError, KeyError):
        return annotations  # contig absent from ClinVar (e.g. an unplaced scaffold)

    for record in records:
        if record.pos != variant.pos or variant.alt not in (record.alts or ()):
            continue

        clnsig = ",".join(record.info.get("CLNSIG", ()))
        gene_info = record.info.get("GENEINFO", "")
        publications = [citation for citation in record.info.get("CLNVI", ()) if citation]

        if gene_info and not variant.gene:
            variant.gene = gene_info.split(":")[0]

        for disease_name in _clinvar_diseases(record):
            annotations.append(
                Annotation(
                    variant_id=variant.id,
                    source="clinvar",
                    disease=disease_name,
                    clinical_significance=clnsig or None,
                    publications=publications or None,
                    evidence_strength=_clinvar_strength(clnsig),
                )
            )
    return annotations


def _clinvar_diseases(record: pysam.VariantRecord) -> list[str]:
    """Extract normalized disease names from a ClinVar record's CLNDN field.

    CLNDN separates multiple diseases with '|' and encodes spaces as '_'. Individual disease names
    legitimately contain commas ("Breast-ovarian_cancer,_familial,_susceptibility_to,_1"), but
    pysam splits every INFO list field on commas - so the tuple must be rejoined before splitting
    on '|', or one disease fragments into several junk entries ("familial", "1", ...).
    """
    raw = ",".join(record.info.get("CLNDN", ()))
    diseases = []
    for name in raw.split("|"):
        if not name or name.lower() in {"not_provided", "not_specified"}:
            continue
        diseases.append(normalize_disease(name))
    return diseases


def _gwas_annotations(session: Session, variant: Variant) -> list[Annotation]:
    """GWAS Catalog associations within GWAS_WINDOW bp of one variant."""
    rows = session.scalars(
        select(RefGwasAssociation).where(
            RefGwasAssociation.chrom == variant.chrom,
            RefGwasAssociation.pos >= variant.pos - GWAS_WINDOW,
            RefGwasAssociation.pos <= variant.pos + GWAS_WINDOW,
        )
    ).all()

    annotations: list[Annotation] = []
    for row in rows:
        if row.gene and not variant.gene:
            variant.gene = row.gene.split(",")[0].strip()
        annotations.append(
            Annotation(
                variant_id=variant.id,
                source="gwas",
                disease=normalize_disease(row.disease_trait),
                clinical_significance=None,
                publications=[f"PMID:{row.pubmed_id}"] if row.pubmed_id else None,
                evidence_strength=_gwas_strength(row.p_value),
            )
        )
    return annotations


def _gwas_strength(p_value: float | None) -> float:
    """Map a GWAS p-value onto [0, 1]. Genome-wide significance (5e-8) lands at ~0.5, and
    stronger associations approach 1.0.
    """
    if p_value is None or p_value <= 0:
        return 0.0
    import math

    neg_log_p = -math.log10(p_value)
    return min(neg_log_p / 15.0, 1.0)
