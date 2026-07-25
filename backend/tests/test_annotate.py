"""Tests for ClinVar + GWAS annotation.

The disease-name tests matter more than they look: ClinVar disease names contain commas, and pysam
splits INFO list fields on commas, so a naive read fragments one disease into several junk entries.
"""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from backend.annotation.annotate import (
    _clinvar_diseases,
    _clinvar_strength,
    _gwas_annotations,
    _gwas_strength,
    normalize_disease,
)
from backend.database.models import RefGwasAssociation
from backend.tests.conftest import make_variant


def _record(clndn: tuple[str, ...]) -> SimpleNamespace:
    return SimpleNamespace(info={"CLNDN": clndn})


def test_comma_containing_disease_stays_one_disease():
    """The exact shape pysam hands back for 'Breast-ovarian_cancer,_familial,_susceptibility_to,_1'."""
    record = _record(("Breast-ovarian_cancer", "_familial", "_susceptibility_to", "_1"))
    assert _clinvar_diseases(record) == ["breast-ovarian cancer, familial, susceptibility to, 1"]


def test_pipe_separates_multiple_diseases():
    record = _record(("Lung_cancer|Breast_cancer",))
    assert _clinvar_diseases(record) == ["lung cancer", "breast cancer"]


def test_pipe_and_comma_together():
    record = _record(("Cancer_A,_type_1|Cancer_B",))
    assert _clinvar_diseases(record) == ["cancer a, type 1", "cancer b"]


def test_placeholder_diseases_dropped():
    assert _clinvar_diseases(_record(("not_provided",))) == []
    assert _clinvar_diseases(_record(("not_specified|Real_disease",))) == ["real disease"]
    assert _clinvar_diseases(_record(())) == []


def test_normalize_disease_collapses_whitespace_and_case():
    assert normalize_disease("  Breast   Cancer  ") == "breast cancer"
    assert normalize_disease("Long_QT_syndrome") == "long qt syndrome"
    # Idempotent, since the ranking stage may re-apply it.
    assert normalize_disease(normalize_disease("Long_QT_syndrome")) == "long qt syndrome"


@pytest.mark.parametrize(
    ("clnsig", "expected"),
    [
        ("Pathogenic", 1.0),
        ("Pathogenic/Likely_pathogenic", 0.95),
        ("Likely_pathogenic", 0.8),
        ("Benign", 0.1),
        ("Uncertain_significance", 0.0),
        ("Conflicting_classifications_of_pathogenicity", 0.0),
        ("", 0.0),
    ],
)
def test_clinvar_strength(clnsig, expected):
    assert _clinvar_strength(clnsig) == expected


def test_gwas_strength_scales_with_significance():
    assert _gwas_strength(None) == 0.0
    assert _gwas_strength(0.0) == 0.0
    # Genome-wide significance should register as real but not maximal evidence.
    assert 0.3 < _gwas_strength(5e-8) < 0.7
    # A far stronger association saturates at 1.0.
    assert _gwas_strength(1e-30) == 1.0
    assert _gwas_strength(1e-9) > _gwas_strength(1e-6)


def test_gwas_annotation_matches_within_window(session, analysis_run):
    variant = make_variant(session, analysis_run, chrom="7", pos=100_000, gene=None)
    session.add_all([
        RefGwasAssociation(
            chrom="7", pos=100_050, gene="EGFR", disease_trait="Lung carcinoma", p_value=1e-12
        ),
        RefGwasAssociation(  # too far away to count
            chrom="7", pos=200_000, gene="OTHER", disease_trait="Unrelated trait", p_value=1e-12
        ),
        RefGwasAssociation(  # right position, wrong chromosome
            chrom="8", pos=100_050, gene="OTHER", disease_trait="Unrelated trait", p_value=1e-12
        ),
    ])
    session.commit()

    annotations = _gwas_annotations(session, variant)
    assert [a.disease for a in annotations] == ["lung carcinoma"]
    assert annotations[0].source == "gwas"
    # A gene-less variant should pick up the gene symbol from the matching association.
    assert variant.gene == "EGFR"
