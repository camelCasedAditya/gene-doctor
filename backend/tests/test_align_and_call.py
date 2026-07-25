"""Integration test for the alignment -> variant calling -> VCF loading chain.

Uses a small synthetic reference and a copy of it carrying known SNVs, so the whole shell-out chain
(minimap2, samtools, bcftools) is exercised without touching the real 3GB GRCh38 reference.
"""

from __future__ import annotations

import random
import subprocess

import pysam
import pytest

from backend.database.models import Variant
from backend.genomics.align import align_genome
from backend.genomics.load_vcf import load_vcf
from backend.genomics.variant_call import call_variants
from backend.settings import settings

REF_LEN = 60_000
SNV_OFFSETS = (20_000, 40_000)
# Substitution cycle, so the alt base is always different from whatever the reference holds.
_SUBSTITUTE = {"A": "C", "C": "G", "G": "T", "T": "A"}


def _write_fasta(path, name: str, sequence: str) -> None:
    with open(path, "w") as f:
        f.write(f">{name}\n")
        for i in range(0, len(sequence), 60):
            f.write(sequence[i : i + 60] + "\n")


@pytest.fixture
def synthetic_reference(tmp_path, monkeypatch):
    """A random-but-alignable reference plus a mutated 'patient' copy carrying known SNVs."""
    rng = random.Random(1234)
    ref_seq = "".join(rng.choice("ACGT") for _ in range(REF_LEN))

    ref_path = tmp_path / "ref.fa"
    _write_fasta(ref_path, "21", ref_seq)
    subprocess.run(["samtools", "faidx", str(ref_path)], check=True)
    mmi_path = tmp_path / "ref.mmi"
    subprocess.run(
        ["minimap2", "-x", "asm5", "-d", str(mmi_path), str(ref_path)],
        check=True,
        capture_output=True,
    )

    patient = list(ref_seq)
    expected: dict[int, str] = {}
    for offset in SNV_OFFSETS:
        alt = _SUBSTITUTE[patient[offset]]
        patient[offset] = alt
        expected[offset + 1] = alt  # VCF positions are 1-based
    patient_path = tmp_path / "patient.fa"
    _write_fasta(patient_path, "21", "".join(patient))

    monkeypatch.setattr(settings, "reference_fasta", ref_path)
    monkeypatch.setattr(settings, "reference_mmi", mmi_path)
    return patient_path, expected


def test_align_produces_indexed_bam(synthetic_reference, tmp_path):
    patient_path, _ = synthetic_reference
    bam = align_genome(str(patient_path), tmp_path / "work")
    assert bam.exists()
    assert bam.with_suffix(".bam.bai").exists()
    with pysam.AlignmentFile(str(bam)) as af:
        assert af.mapped > 0


def test_call_variants_recovers_known_snvs(synthetic_reference, tmp_path):
    patient_path, expected = synthetic_reference
    bam = align_genome(str(patient_path), tmp_path / "work")
    vcf = call_variants(bam, tmp_path / "work")
    assert vcf.exists()

    called = {}
    with pysam.VariantFile(str(vcf)) as vf:
        for record in vf:
            called[record.pos] = record.alts[0]

    for pos, alt in expected.items():
        assert pos in called, f"expected a call at {pos}, got {sorted(called)}"
        assert called[pos] == alt


def test_load_vcf_creates_variant_rows(synthetic_reference, tmp_path, session, analysis_run):
    patient_path, expected = synthetic_reference
    bam = align_genome(str(patient_path), tmp_path / "work")
    vcf = call_variants(bam, tmp_path / "work")

    count = load_vcf(session, analysis_run.id, vcf)
    assert count == len(expected)

    variants = session.query(Variant).filter_by(analysis_run_id=analysis_run.id).all()
    assert {v.pos for v in variants} == set(expected)
    assert all(v.chrom == "21" for v in variants)
