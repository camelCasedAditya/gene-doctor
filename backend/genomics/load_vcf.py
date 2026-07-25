"""Load a called/normalized VCF into Variant rows for an analysis run."""

from __future__ import annotations

import uuid
from pathlib import Path

import pysam
from sqlalchemy.orm import Session

from backend.database.models import Variant

BATCH_SIZE = 5000


def load_vcf(session: Session, analysis_run_id: str, vcf_path: Path) -> int:
    """Insert one Variant row per ALT allele in the VCF. Returns the number of variants loaded."""
    batch: list[dict] = []
    count = 0

    with pysam.VariantFile(str(vcf_path)) as vcf:
        for record in vcf:
            for alt in record.alts or ():
                batch.append({
                    "id": uuid.uuid4().hex,
                    "analysis_run_id": analysis_run_id,
                    "chrom": record.chrom,
                    "pos": record.pos,
                    "ref": record.ref,
                    "alt": alt,
                    "gene": None,
                    "transcript_id": None,
                })
                count += 1
                if len(batch) >= BATCH_SIZE:
                    session.execute(Variant.__table__.insert(), batch)
                    session.commit()
                    batch.clear()

    if batch:
        session.execute(Variant.__table__.insert(), batch)
        session.commit()

    return count
