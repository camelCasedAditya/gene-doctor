"""Drop and recreate the per-analysis result tables, preserving the bulk-ingested reference tables.

The ref_* tables take minutes to rebuild (974k GWAS associations, 388k transcripts), so a schema
change to a result table shouldn't force a full re-ingest. SQLAlchemy's create_all only creates
missing tables - it never ALTERs an existing one - so a model change needs the stale table dropped.

Run via: python -m backend.database.reset_results
"""

from __future__ import annotations

from sqlalchemy import create_engine

from backend.database.models import (
    AiPrediction,
    AnalysisRun,
    Annotation,
    Aso,
    Base,
    Disease,
    DiseaseVariantLink,
    EvidenceWeights,
    GenomeUpload,
    Transcript,
    Variant,
    VariantScore,
)
from backend.settings import settings

# Dropped in FK-safe order (children before parents).
RESULT_TABLES = [
    Aso,
    DiseaseVariantLink,
    Disease,
    Transcript,
    VariantScore,
    AiPrediction,
    Annotation,
    Variant,
    AnalysisRun,
    GenomeUpload,
    EvidenceWeights,
]


def reset_results() -> None:
    engine = create_engine(settings.database_url.replace("sqlite+aiosqlite", "sqlite"))
    tables = [model.__table__ for model in RESULT_TABLES]
    Base.metadata.drop_all(engine, tables=tables)
    Base.metadata.create_all(engine, tables=tables)


if __name__ == "__main__":
    reset_results()
    print("Result tables reset; reference tables left intact.")
