"""One-time ingestion of the Ensembl GRCh38 GTF into ref_ensembl_transcripts, one row per
transcript with its exon structure as JSON.

Run via: python -m backend.annotation.ingest_ensembl_gtf
"""

from __future__ import annotations

import gzip
import re

from sqlalchemy import create_engine, delete
from sqlalchemy.orm import Session

from backend.database.models import Base, RefEnsemblTranscript
from backend.settings import settings

BATCH_SIZE = 2000
_ATTR_RE = re.compile(r'(\w+) "([^"]*)"')


def _parse_attributes(field: str) -> dict[str, str]:
    return dict(_ATTR_RE.findall(field))


def ingest_ensembl_gtf() -> int:
    transcripts: dict[str, dict] = {}

    with gzip.open(settings.ensembl_gtf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            chrom, _source, feature, start, end, _score, strand, _frame, attr_field = line.rstrip(
                "\n"
            ).split("\t", 8)

            if feature == "transcript":
                attrs = _parse_attributes(attr_field)
                transcript_id = attrs.get("transcript_id")
                if not transcript_id:
                    continue
                transcripts[transcript_id] = {
                    "gene": attrs.get("gene_name", attrs.get("gene_id", "")),
                    "ensembl_gene_id": attrs.get("gene_id", ""),
                    "ensembl_transcript_id": transcript_id,
                    "chrom": chrom,
                    "start": int(start),
                    "end": int(end),
                    "strand": strand,
                    "exons": [],
                }
            elif feature == "exon":
                attrs = _parse_attributes(attr_field)
                transcript_id = attrs.get("transcript_id")
                if transcript_id not in transcripts:
                    continue
                transcripts[transcript_id]["exons"].append({
                    "start": int(start),
                    "end": int(end),
                    "exon_number": attrs.get("exon_number", ""),
                })

    sync_url = settings.database_url.replace("sqlite+aiosqlite", "sqlite")
    engine = create_engine(sync_url)
    Base.metadata.create_all(engine)

    count = 0
    batch: list[dict] = []
    with Session(engine) as session:
        session.execute(delete(RefEnsemblTranscript))
        for row in transcripts.values():
            row["exons"].sort(key=lambda e: e["start"])
            batch.append(row)
            count += 1
            if len(batch) >= BATCH_SIZE:
                session.execute(RefEnsemblTranscript.__table__.insert(), batch)
                session.commit()
                batch.clear()
        if batch:
            session.execute(RefEnsemblTranscript.__table__.insert(), batch)
            session.commit()

    return count


if __name__ == "__main__":
    n = ingest_ensembl_gtf()
    print(f"Ingested {n} Ensembl transcripts")
