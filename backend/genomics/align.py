"""Align an assembled genome FASTA against the GRCh38 reference with minimap2, producing a
sorted, indexed BAM. Uses the asm5 preset (assembly-to-reference), since the uploaded FASTA is
an assembled per-chromosome genome, not raw sequencing reads.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

from backend.settings import settings


class AlignmentError(RuntimeError):
    pass


def _run(cmd: list[str]) -> None:
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise AlignmentError(f"Command failed: {' '.join(cmd)}\n{result.stderr}")


def align_genome(query_fasta: str, output_dir: Path) -> Path:
    """Align query_fasta against the GRCh38 reference. Returns the path to the sorted, indexed BAM."""
    output_dir.mkdir(parents=True, exist_ok=True)
    sam_path = output_dir / "aligned.sam"
    bam_path = output_dir / "aligned.bam"
    sorted_bam_path = output_dir / "aligned.sorted.bam"

    _run([
        "minimap2",
        "-a",
        "-t",
        "4",
        str(settings.reference_mmi),
        query_fasta,
        "-o",
        str(sam_path),
    ])
    _run(["samtools", "view", "-b", "-o", str(bam_path), str(sam_path)])
    _run(["samtools", "sort", "-o", str(sorted_bam_path), str(bam_path)])
    _run(["samtools", "index", str(sorted_bam_path)])

    sam_path.unlink(missing_ok=True)
    bam_path.unlink(missing_ok=True)

    return sorted_bam_path
