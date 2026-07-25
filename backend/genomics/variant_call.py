"""Call variants from a sorted/indexed BAM against the GRCh38 reference using bcftools, then
normalize (left-align indels, split multiallelics) into a clean VCF.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

from backend.settings import settings


class VariantCallError(RuntimeError):
    pass


def _run(cmd: list[str], output_path: Path) -> None:
    with open(output_path, "wb") as out_file:
        result = subprocess.run(cmd, stdout=out_file, stderr=subprocess.PIPE)
    if result.returncode != 0:
        raise VariantCallError(f"Command failed: {' '.join(cmd)}\n{result.stderr.decode()}")


def call_variants(bam_path: Path, output_dir: Path) -> Path:
    """Run bcftools mpileup | call | norm on bam_path. Returns the path to the normalized VCF."""
    output_dir.mkdir(parents=True, exist_ok=True)
    pileup_vcf = output_dir / "pileup.vcf"
    called_vcf = output_dir / "called.vcf"
    normalized_vcf = output_dir / "variants.vcf"

    _run(["bcftools", "mpileup", "-f", str(settings.reference_fasta), str(bam_path)], pileup_vcf)
    _run(["bcftools", "call", "-mv", "-Ov", str(pileup_vcf)], called_vcf)
    _run(["bcftools", "norm", "-f", str(settings.reference_fasta), str(called_vcf)], normalized_vcf)

    pileup_vcf.unlink(missing_ok=True)
    called_vcf.unlink(missing_ok=True)
    return normalized_vcf
