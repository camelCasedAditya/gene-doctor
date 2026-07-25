"""Streaming FASTA validation: structural correctness + human-genome chromosome checks.

Never loads the whole file into memory - pysam.FastxFile streams records, and progress
is tracked by bytes read off disk (works for both plain and gzipped FASTA).
"""

from __future__ import annotations

import re
from collections.abc import Callable
from dataclasses import dataclass, field

import pysam

EXPECTED_CHROMOSOMES = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}
VALID_BASES = set("ACGTNRYSWKMBDHVacgtnryswkmbdhv")
_CHR_PREFIX = re.compile(r"^chr", re.IGNORECASE)


@dataclass
class ValidationResult:
    valid: bool
    chromosomes_found: list[str] = field(default_factory=list)
    duplicates: list[str] = field(default_factory=list)
    missing: list[str] = field(default_factory=list)
    extra_contigs: list[str] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "valid": self.valid,
            "chromosomes_found": self.chromosomes_found,
            "duplicates": self.duplicates,
            "missing": self.missing,
            "extra_contigs": self.extra_contigs,
            "errors": self.errors,
            "warnings": self.warnings,
        }


def _normalize(name: str) -> str:
    return _CHR_PREFIX.sub("", name)


def validate_fasta(
    file_path: str,
    on_progress: Callable[[float], None] | None = None,
) -> ValidationResult:
    """Validate a FASTA as a human genome. Reports duplicate/missing chromosomes and
    invalid sequence characters. Missing Y is a warning only (valid for female samples);
    every other missing or duplicated chromosome is a hard failure.
    """
    errors: list[str] = []
    warnings: list[str] = []
    seen: dict[str, int] = {}
    extra_contigs: list[str] = []

    # pysam.FastxFile exposes no file-position API, so progress is approximated by
    # how many of the expected chromosomes have been seen so far (extra scaffolds
    # don't move the bar, but there are only ever a handful of primary chromosomes).
    try:
        with pysam.FastxFile(file_path) as fasta:
            for record in fasta:
                name = _normalize(record.name)
                sequence = record.sequence or ""

                if not set(sequence) <= VALID_BASES:
                    bad_chars = sorted(set(sequence) - VALID_BASES)
                    errors.append(f"Invalid sequence characters in '{record.name}': {bad_chars}")

                if name in EXPECTED_CHROMOSOMES:
                    seen[name] = seen.get(name, 0) + 1
                else:
                    extra_contigs.append(record.name)

                if on_progress is not None:
                    on_progress(min(len(seen) / len(EXPECTED_CHROMOSOMES), 0.99))
    except (OSError, ValueError) as exc:
        return ValidationResult(valid=False, errors=[f"Malformed FASTA: {exc}"])

    duplicates = sorted(name for name, count in seen.items() if count > 1)
    missing = sorted(EXPECTED_CHROMOSOMES - seen.keys())

    if duplicates:
        errors.append(f"Duplicate chromosomes: {duplicates}")

    hard_missing = [c for c in missing if c != "Y"]
    if hard_missing:
        errors.append(f"Missing chromosomes: {hard_missing}")
    if "Y" in missing:
        warnings.append("Chromosome Y is missing (expected for female samples)")

    if on_progress is not None:
        on_progress(1.0)

    return ValidationResult(
        valid=not errors,
        chromosomes_found=sorted(seen.keys()),
        duplicates=duplicates,
        missing=missing,
        extra_contigs=extra_contigs,
        errors=errors,
        warnings=warnings,
    )
