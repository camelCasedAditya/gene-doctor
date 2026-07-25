"""One-time ingestion of the GWAS Catalog associations TSV into ref_gwas_associations.

Run via: python -m backend.annotation.ingest_gwas
"""

from __future__ import annotations

import csv
import sys
import uuid

from sqlalchemy import create_engine, delete
from sqlalchemy.orm import Session

from backend.database.models import Base, RefGwasAssociation
from backend.settings import settings

BATCH_SIZE = 5000


def _to_float(value: str) -> float | None:
    try:
        return float(value)
    except ValueError:
        return None


def _to_int(value: str) -> int | None:
    try:
        return int(value)
    except ValueError:
        return None


def ingest_gwas() -> int:
    csv.field_size_limit(sys.maxsize)
    sync_url = settings.database_url.replace("sqlite+aiosqlite", "sqlite")
    engine = create_engine(sync_url)
    Base.metadata.create_all(engine)

    count = 0
    batch: list[dict] = []

    with Session(engine) as session:
        session.execute(delete(RefGwasAssociation))

        with open(settings.gwas_associations_tsv, newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                chrom = row.get("CHR_ID", "").strip()
                pos = _to_int(row.get("CHR_POS", "").strip())
                disease_trait = row.get("DISEASE/TRAIT", "").strip()
                if not chrom or pos is None or not disease_trait or "x" in chrom.lower():
                    continue  # skip multi-chromosome / unmapped associations

                batch.append({
                    "id": uuid.uuid4().hex,
                    "chrom": chrom,
                    "pos": pos,
                    "gene": row.get("MAPPED_GENE", "").strip() or None,
                    "disease_trait": disease_trait,
                    "p_value": _to_float(row.get("P-VALUE", "").strip()),
                    "pubmed_id": row.get("PUBMEDID", "").strip() or None,
                    "study_accession": row.get("STUDY ACCESSION", "").strip() or None,
                })
                count += 1

                if len(batch) >= BATCH_SIZE:
                    session.execute(RefGwasAssociation.__table__.insert(), batch)
                    session.commit()
                    batch.clear()

            if batch:
                session.execute(RefGwasAssociation.__table__.insert(), batch)
                session.commit()

    return count


if __name__ == "__main__":
    n = ingest_gwas()
    print(f"Ingested {n} GWAS associations")
