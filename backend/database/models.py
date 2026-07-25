import uuid
from datetime import datetime, timezone

from sqlalchemy import JSON, Float, ForeignKey, Integer, String, UniqueConstraint
from sqlalchemy.orm import DeclarativeBase, Mapped, mapped_column, relationship


def _uuid() -> str:
    return str(uuid.uuid4())


def _now() -> datetime:
    return datetime.now(timezone.utc)


class Base(DeclarativeBase):
    pass


class GenomeUpload(Base):
    __tablename__ = "genome_uploads"

    id: Mapped[str] = mapped_column(String, primary_key=True, default=_uuid)
    file_path: Mapped[str] = mapped_column(String)
    status: Mapped[str] = mapped_column(String, default="pending")  # pending|validating|valid|invalid|error
    progress: Mapped[float] = mapped_column(Float, default=0.0)
    result: Mapped[dict | None] = mapped_column(JSON, nullable=True)
    created_at: Mapped[datetime] = mapped_column(default=_now)
    updated_at: Mapped[datetime] = mapped_column(default=_now, onupdate=_now)


class AnalysisRun(Base):
    __tablename__ = "analysis_runs"

    id: Mapped[str] = mapped_column(String, primary_key=True, default=_uuid)
    genome_upload_id: Mapped[str] = mapped_column(ForeignKey("genome_uploads.id"))
    status: Mapped[str] = mapped_column(String, default="queued")
    # queued|aligning|calling|annotating|predicting|aggregating|ranking|designing_aso|done|failed
    current_stage: Mapped[str | None] = mapped_column(String, nullable=True)
    error: Mapped[str | None] = mapped_column(String, nullable=True)
    created_at: Mapped[datetime] = mapped_column(default=_now)
    updated_at: Mapped[datetime] = mapped_column(default=_now, onupdate=_now)

    variants: Mapped[list["Variant"]] = relationship(back_populates="analysis_run")


class Variant(Base):
    __tablename__ = "variants"

    id: Mapped[str] = mapped_column(String, primary_key=True, default=_uuid)
    analysis_run_id: Mapped[str] = mapped_column(ForeignKey("analysis_runs.id"))
    chrom: Mapped[str] = mapped_column(String, index=True)
    pos: Mapped[int] = mapped_column(Integer, index=True)
    ref: Mapped[str] = mapped_column(String)
    alt: Mapped[str] = mapped_column(String)
    gene: Mapped[str | None] = mapped_column(String, nullable=True, index=True)
    transcript_id: Mapped[str | None] = mapped_column(String, nullable=True)

    analysis_run: Mapped["AnalysisRun"] = relationship(back_populates="variants")
    annotations: Mapped[list["Annotation"]] = relationship(back_populates="variant")
    ai_predictions: Mapped[list["AiPrediction"]] = relationship(back_populates="variant")
    score: Mapped["VariantScore | None"] = relationship(back_populates="variant", uselist=False)


class Annotation(Base):
    """Validated clinical evidence only (ClinVar, GWAS Catalog). Never mixed with AiPrediction rows."""

    __tablename__ = "annotations"

    id: Mapped[str] = mapped_column(String, primary_key=True, default=_uuid)
    variant_id: Mapped[str] = mapped_column(ForeignKey("variants.id"))
    source: Mapped[str] = mapped_column(String)  # "clinvar" | "gwas"
    disease: Mapped[str | None] = mapped_column(String, nullable=True, index=True)
    clinical_significance: Mapped[str | None] = mapped_column(String, nullable=True)
    publications: Mapped[list | None] = mapped_column(JSON, nullable=True)
    evidence_strength: Mapped[float | None] = mapped_column(Float, nullable=True)

    variant: Mapped["Variant"] = relationship(back_populates="annotations")


class AiPrediction(Base):
    """AI-generated research hypotheses only (variant_lm, alphamissense). Never mixed with Annotation rows."""

    __tablename__ = "ai_predictions"

    id: Mapped[str] = mapped_column(String, primary_key=True, default=_uuid)
    variant_id: Mapped[str] = mapped_column(ForeignKey("variants.id"))
    source: Mapped[str] = mapped_column(String)  # "variant_lm" | "alphamissense"
    functional_impact: Mapped[float | None] = mapped_column(Float, nullable=True)
    regulatory_impact: Mapped[float | None] = mapped_column(Float, nullable=True)
    pathogenicity_probability: Mapped[float | None] = mapped_column(Float, nullable=True)
    confidence: Mapped[float | None] = mapped_column(Float, nullable=True)

    variant: Mapped["Variant"] = relationship(back_populates="ai_predictions")


class EvidenceWeights(Base):
    """Single-row config table for the weighted evidence engine."""

    __tablename__ = "evidence_weights"

    id: Mapped[int] = mapped_column(Integer, primary_key=True, default=1)
    clinvar_weight: Mapped[float] = mapped_column(Float, default=0.45)
    gwas_weight: Mapped[float] = mapped_column(Float, default=0.30)
    ai_weight: Mapped[float] = mapped_column(Float, default=0.25)


class VariantScore(Base):
    __tablename__ = "variant_scores"

    id: Mapped[str] = mapped_column(String, primary_key=True, default=_uuid)
    variant_id: Mapped[str] = mapped_column(ForeignKey("variants.id"), unique=True)
    overall_score: Mapped[float] = mapped_column(Float)

    variant: Mapped["Variant"] = relationship(back_populates="score")


class Disease(Base):
    __tablename__ = "diseases"

    id: Mapped[str] = mapped_column(String, primary_key=True, default=_uuid)
    analysis_run_id: Mapped[str] = mapped_column(ForeignKey("analysis_runs.id"), index=True)
    name: Mapped[str] = mapped_column(String, index=True)  # normalized (lowercase/trimmed)
    description: Mapped[str | None] = mapped_column(String, nullable=True)
    rank_score: Mapped[float | None] = mapped_column(Float, nullable=True)
    confidence: Mapped[float | None] = mapped_column(Float, nullable=True)
    # Counts of each evidence type behind this ranking, so the UI can show validated vs AI
    # provenance without re-querying every linked variant.
    known_variant_count: Mapped[int] = mapped_column(Integer, default=0)
    predicted_variant_count: Mapped[int] = mapped_column(Integer, default=0)

    variant_links: Mapped[list["DiseaseVariantLink"]] = relationship(back_populates="disease")

    __table_args__ = (UniqueConstraint("analysis_run_id", "name"),)


class DiseaseVariantLink(Base):
    __tablename__ = "disease_variant_links"

    id: Mapped[str] = mapped_column(String, primary_key=True, default=_uuid)
    disease_id: Mapped[str] = mapped_column(ForeignKey("diseases.id"))
    variant_id: Mapped[str] = mapped_column(ForeignKey("variants.id"))
    contribution_score: Mapped[float] = mapped_column(Float)

    disease: Mapped["Disease"] = relationship(back_populates="variant_links")

    __table_args__ = (UniqueConstraint("disease_id", "variant_id"),)


class Transcript(Base):
    __tablename__ = "transcripts"

    id: Mapped[str] = mapped_column(String, primary_key=True, default=_uuid)
    analysis_run_id: Mapped[str] = mapped_column(ForeignKey("analysis_runs.id"), index=True)
    gene: Mapped[str] = mapped_column(String, index=True)
    ensembl_transcript_id: Mapped[str] = mapped_column(String, index=True)
    chrom: Mapped[str] = mapped_column(String)
    strand: Mapped[str] = mapped_column(String)
    exon_structure: Mapped[list] = mapped_column(JSON)
    # Exons overlapping high-scoring variants, and the ASO mechanisms they enable.
    candidate_exons: Mapped[list] = mapped_column(JSON, default=list)
    opportunities: Mapped[list] = mapped_column(JSON, default=list)
    rank_score: Mapped[float | None] = mapped_column(Float, nullable=True)
    # "validated" - reached via a disease backed by ClinVar/GWAS evidence.
    # "ai_hypothesis" - reached only via AI-scored variants with no validated evidence at all.
    # Carried through to every ASO designed here, so an AI-only lead is never presented as though
    # it had clinical backing.
    evidence_basis: Mapped[str] = mapped_column(String, default="validated")

    asos: Mapped[list["Aso"]] = relationship(back_populates="transcript")

    __table_args__ = (UniqueConstraint("analysis_run_id", "ensembl_transcript_id"),)


class Aso(Base):
    __tablename__ = "asos"

    id: Mapped[str] = mapped_column(String, primary_key=True, default=_uuid)
    transcript_id: Mapped[str] = mapped_column(ForeignKey("transcripts.id"))
    sequence: Mapped[str] = mapped_column(String)
    genomic_position: Mapped[str] = mapped_column(String)
    target_exon: Mapped[str | None] = mapped_column(String, nullable=True)
    mechanism: Mapped[str] = mapped_column(String)
    # splice_switch | exon_skip | exon_include | knockdown | allele_specific
    gc_pct: Mapped[float] = mapped_column(Float)
    tm: Mapped[float] = mapped_column(Float)
    predicted_efficacy: Mapped[float] = mapped_column(Float)
    predicted_specificity: Mapped[float] = mapped_column(Float)
    off_target_score: Mapped[float] = mapped_column(Float)
    confidence: Mapped[float] = mapped_column(Float)

    transcript: Mapped["Transcript"] = relationship(back_populates="asos")


class RefGwasAssociation(Base):
    """Bulk-ingested GWAS Catalog reference table, queried by position during annotation."""

    __tablename__ = "ref_gwas_associations"

    id: Mapped[str] = mapped_column(String, primary_key=True, default=_uuid)
    chrom: Mapped[str] = mapped_column(String, index=True)
    pos: Mapped[int] = mapped_column(Integer, index=True)
    gene: Mapped[str | None] = mapped_column(String, nullable=True, index=True)
    disease_trait: Mapped[str] = mapped_column(String, index=True)
    p_value: Mapped[float | None] = mapped_column(Float, nullable=True)
    pubmed_id: Mapped[str | None] = mapped_column(String, nullable=True)
    study_accession: Mapped[str | None] = mapped_column(String, nullable=True)


class RefEnsemblTranscript(Base):
    """Bulk-ingested Ensembl GTF reference table: one row per transcript, exons as JSON."""

    __tablename__ = "ref_ensembl_transcripts"

    id: Mapped[str] = mapped_column(String, primary_key=True, default=_uuid)
    gene: Mapped[str] = mapped_column(String, index=True)
    ensembl_gene_id: Mapped[str] = mapped_column(String, index=True)
    ensembl_transcript_id: Mapped[str] = mapped_column(String, unique=True, index=True)
    chrom: Mapped[str] = mapped_column(String, index=True)
    start: Mapped[int] = mapped_column(Integer)
    end: Mapped[int] = mapped_column(Integer)
    strand: Mapped[str] = mapped_column(String)
    exons: Mapped[list] = mapped_column(JSON)  # [{"start": int, "end": int, "exon_number": str}, ...]
