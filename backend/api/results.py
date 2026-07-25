"""Read endpoints for pipeline results: variants, diseases, transcripts, ASOs, and the configurable
evidence weights.

Every response keeps validated evidence (`annotations`) and AI hypotheses (`ai_predictions`) in
separate, explicitly-sourced fields - they are never merged into one list.
"""

from fastapi import APIRouter, Depends, HTTPException, Query
from pydantic import BaseModel
from sqlalchemy import func, select
from sqlalchemy.ext.asyncio import AsyncSession
from sqlalchemy.orm import selectinload

from backend.database.models import (
    Annotation,
    Aso,
    Disease,
    DiseaseVariantLink,
    EvidenceWeights,
    Transcript,
    Variant,
    VariantScore,
)
from backend.database.session import get_session

router = APIRouter(tags=["results"])


class AnnotationOut(BaseModel):
    source: str
    disease: str | None
    clinical_significance: str | None
    publications: list | None
    evidence_strength: float | None


class AiPredictionOut(BaseModel):
    source: str
    functional_impact: float | None
    regulatory_impact: float | None
    pathogenicity_probability: float | None
    confidence: float | None


class VariantOut(BaseModel):
    id: str
    chrom: str
    pos: int
    ref: str
    alt: str
    gene: str | None
    transcript_id: str | None
    overall_score: float | None
    annotations: list[AnnotationOut]
    ai_predictions: list[AiPredictionOut]


class VariantsResponse(BaseModel):
    total: int
    items: list[VariantOut]


def _variant_out(variant: Variant, score: float | None) -> VariantOut:
    return VariantOut(
        id=variant.id,
        chrom=variant.chrom,
        pos=variant.pos,
        ref=variant.ref,
        alt=variant.alt,
        gene=variant.gene,
        transcript_id=variant.transcript_id,
        overall_score=score,
        annotations=[
            AnnotationOut(
                source=a.source,
                disease=a.disease,
                clinical_significance=a.clinical_significance,
                publications=a.publications,
                evidence_strength=a.evidence_strength,
            )
            for a in variant.annotations
        ],
        ai_predictions=[
            AiPredictionOut(
                source=p.source,
                functional_impact=p.functional_impact,
                regulatory_impact=p.regulatory_impact,
                pathogenicity_probability=p.pathogenicity_probability,
                confidence=p.confidence,
            )
            for p in variant.ai_predictions
        ],
    )


@router.get("/variants", response_model=VariantsResponse)
async def list_variants(
    analysis_id: str,
    gene: str | None = None,
    disease: str | None = None,
    chrom: str | None = None,
    evidence: str | None = Query(None, pattern="^(known|predicted)$"),
    min_score: float | None = None,
    limit: int = Query(100, le=1000),
    offset: int = 0,
    session: AsyncSession = Depends(get_session),
) -> VariantsResponse:
    stmt = (
        select(Variant, VariantScore.overall_score)
        .outerjoin(VariantScore, VariantScore.variant_id == Variant.id)
        .options(selectinload(Variant.annotations), selectinload(Variant.ai_predictions))
        .where(Variant.analysis_run_id == analysis_id)
    )

    if gene:
        stmt = stmt.where(Variant.gene.ilike(f"%{gene}%"))
    if chrom:
        stmt = stmt.where(Variant.chrom == chrom)
    if min_score is not None:
        stmt = stmt.where(VariantScore.overall_score >= min_score)
    if disease:
        stmt = stmt.where(Variant.annotations.any(Annotation.disease.ilike(f"%{disease.lower()}%")))
    if evidence == "known":
        stmt = stmt.where(Variant.annotations.any())
    elif evidence == "predicted":
        stmt = stmt.where(~Variant.annotations.any(), Variant.ai_predictions.any())

    total = await session.scalar(
        select(func.count()).select_from(stmt.order_by(None).subquery())
    )
    rows = (
        await session.execute(
            stmt.order_by(VariantScore.overall_score.desc().nullslast())
            .limit(limit)
            .offset(offset)
        )
    ).all()

    return VariantsResponse(
        total=total or 0, items=[_variant_out(variant, score) for variant, score in rows]
    )


class DiseaseOut(BaseModel):
    id: str
    name: str
    description: str | None
    rank_score: float | None
    confidence: float | None
    known_variant_count: int
    predicted_variant_count: int
    genes: list[str]


class DiseasesResponse(BaseModel):
    items: list[DiseaseOut]


class DiseaseDetail(DiseaseOut):
    supporting_variants: list[VariantOut]


async def _disease_genes(session: AsyncSession, disease_id: str) -> list[str]:
    genes = (
        await session.scalars(
            select(Variant.gene)
            .join(DiseaseVariantLink, DiseaseVariantLink.variant_id == Variant.id)
            .where(DiseaseVariantLink.disease_id == disease_id, Variant.gene.isnot(None))
            .distinct()
        )
    ).all()
    return [g for g in genes if g]


@router.get("/diseases", response_model=DiseasesResponse)
async def list_diseases(
    analysis_id: str, session: AsyncSession = Depends(get_session)
) -> DiseasesResponse:
    diseases = (
        await session.scalars(
            select(Disease)
            .where(Disease.analysis_run_id == analysis_id)
            .order_by(Disease.rank_score.desc().nullslast())
        )
    ).all()

    items = []
    for disease in diseases:
        items.append(
            DiseaseOut(
                id=disease.id,
                name=disease.name,
                description=disease.description,
                rank_score=disease.rank_score,
                confidence=disease.confidence,
                known_variant_count=disease.known_variant_count,
                predicted_variant_count=disease.predicted_variant_count,
                genes=await _disease_genes(session, disease.id),
            )
        )
    return DiseasesResponse(items=items)


@router.get("/diseases/{disease_id}", response_model=DiseaseDetail)
async def get_disease(
    disease_id: str, session: AsyncSession = Depends(get_session)
) -> DiseaseDetail:
    disease = await session.get(Disease, disease_id)
    if disease is None:
        raise HTTPException(status_code=404, detail="Disease not found")

    rows = (
        await session.execute(
            select(Variant, VariantScore.overall_score)
            .join(DiseaseVariantLink, DiseaseVariantLink.variant_id == Variant.id)
            .outerjoin(VariantScore, VariantScore.variant_id == Variant.id)
            .options(selectinload(Variant.annotations), selectinload(Variant.ai_predictions))
            .where(DiseaseVariantLink.disease_id == disease_id)
            .order_by(VariantScore.overall_score.desc().nullslast())
        )
    ).all()

    return DiseaseDetail(
        id=disease.id,
        name=disease.name,
        description=disease.description,
        rank_score=disease.rank_score,
        confidence=disease.confidence,
        known_variant_count=disease.known_variant_count,
        predicted_variant_count=disease.predicted_variant_count,
        genes=await _disease_genes(session, disease.id),
        supporting_variants=[_variant_out(v, s) for v, s in rows],
    )


class TranscriptOut(BaseModel):
    id: str
    gene: str
    ensembl_transcript_id: str
    chrom: str
    strand: str
    exon_structure: list
    candidate_exons: list
    opportunities: list
    rank_score: float | None
    evidence_basis: str


class TranscriptsResponse(BaseModel):
    items: list[TranscriptOut]


@router.get("/transcripts", response_model=TranscriptsResponse)
async def list_transcripts(
    analysis_id: str, session: AsyncSession = Depends(get_session)
) -> TranscriptsResponse:
    transcripts = (
        await session.scalars(
            select(Transcript)
            .where(Transcript.analysis_run_id == analysis_id)
            .order_by(Transcript.rank_score.desc().nullslast())
        )
    ).all()
    return TranscriptsResponse(
        items=[
            TranscriptOut(
                id=t.id,
                gene=t.gene,
                ensembl_transcript_id=t.ensembl_transcript_id,
                chrom=t.chrom,
                strand=t.strand,
                exon_structure=t.exon_structure,
                candidate_exons=t.candidate_exons,
                opportunities=t.opportunities,
                rank_score=t.rank_score,
                evidence_basis=t.evidence_basis,
            )
            for t in transcripts
        ]
    )


class AsoOut(BaseModel):
    id: str
    transcript_id: str
    ensembl_transcript_id: str
    gene: str
    sequence: str
    genomic_position: str
    target_exon: str | None
    mechanism: str
    gc_pct: float
    tm: float
    predicted_efficacy: float
    predicted_specificity: float
    off_target_score: float
    confidence: float
    evidence_basis: str


class AsosResponse(BaseModel):
    items: list[AsoOut]


@router.get("/asos", response_model=AsosResponse)
async def list_asos(
    analysis_id: str,
    mechanism: str | None = None,
    transcript_id: str | None = None,
    session: AsyncSession = Depends(get_session),
) -> AsosResponse:
    stmt = (
        select(Aso, Transcript.gene, Transcript.ensembl_transcript_id, Transcript.evidence_basis)
        .join(Transcript, Aso.transcript_id == Transcript.id)
        .where(Transcript.analysis_run_id == analysis_id)
    )
    if mechanism:
        stmt = stmt.where(Aso.mechanism == mechanism)
    if transcript_id:
        stmt = stmt.where(Aso.transcript_id == transcript_id)

    rows = (await session.execute(stmt.order_by(Aso.predicted_efficacy.desc()))).all()
    return AsosResponse(
        items=[
            AsoOut(
                id=aso.id,
                transcript_id=aso.transcript_id,
                ensembl_transcript_id=ensembl_id,
                gene=gene,
                sequence=aso.sequence,
                genomic_position=aso.genomic_position,
                target_exon=aso.target_exon,
                mechanism=aso.mechanism,
                gc_pct=aso.gc_pct,
                tm=aso.tm,
                predicted_efficacy=aso.predicted_efficacy,
                predicted_specificity=aso.predicted_specificity,
                off_target_score=aso.off_target_score,
                confidence=aso.confidence,
                evidence_basis=basis,
            )
            for aso, gene, ensembl_id, basis in rows
        ]
    )


class WeightsOut(BaseModel):
    clinvar_weight: float
    gwas_weight: float
    ai_weight: float


@router.get("/config/weights", response_model=WeightsOut)
async def get_config_weights(session: AsyncSession = Depends(get_session)) -> WeightsOut:
    weights = await session.get(EvidenceWeights, 1)
    if weights is None:
        weights = EvidenceWeights(id=1)
        session.add(weights)
        await session.commit()
        await session.refresh(weights)
    return WeightsOut(
        clinvar_weight=weights.clinvar_weight,
        gwas_weight=weights.gwas_weight,
        ai_weight=weights.ai_weight,
    )


@router.put("/config/weights", response_model=WeightsOut)
async def put_config_weights(
    body: WeightsOut, session: AsyncSession = Depends(get_session)
) -> WeightsOut:
    weights = await session.get(EvidenceWeights, 1)
    if weights is None:
        weights = EvidenceWeights(id=1)
        session.add(weights)
    weights.clinvar_weight = body.clinvar_weight
    weights.gwas_weight = body.gwas_weight
    weights.ai_weight = body.ai_weight
    await session.commit()
    return body
