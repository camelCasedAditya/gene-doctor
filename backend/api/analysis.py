import asyncio

from fastapi import APIRouter, Depends, HTTPException
from pydantic import BaseModel
from sqlalchemy import func, select
from sqlalchemy.ext.asyncio import AsyncSession

from backend.database.models import (
    AiPrediction,
    AnalysisRun,
    Annotation,
    Aso,
    Disease,
    GenomeUpload,
    Transcript,
    Variant,
)
from backend.database.session import get_session
from backend.services.pipeline import run_pipeline

router = APIRouter(tags=["analysis"])


class AnalyzeRequest(BaseModel):
    genome_upload_id: str


class AnalysisCounts(BaseModel):
    variants: int = 0
    known_variants: int = 0
    unknown_variants: int = 0
    diseases: int = 0
    transcripts: int = 0
    asos: int = 0


class AnalysisStatus(BaseModel):
    id: str
    genome_upload_id: str
    status: str
    current_stage: str | None = None
    error: str | None = None
    counts: AnalysisCounts = AnalysisCounts()


@router.post("/analyze", response_model=AnalysisStatus)
async def start_analysis(
    req: AnalyzeRequest, session: AsyncSession = Depends(get_session)
) -> AnalysisStatus:
    upload = await session.get(GenomeUpload, req.genome_upload_id)
    if upload is None:
        raise HTTPException(status_code=404, detail="Genome upload not found")
    if upload.status != "valid":
        raise HTTPException(
            status_code=400,
            detail=f"Genome upload must be validated first (status: {upload.status})",
        )

    run = AnalysisRun(genome_upload_id=upload.id, status="queued")
    session.add(run)
    await session.commit()
    await session.refresh(run)

    # The pipeline is CPU/IO-heavy and synchronous (minimap2, bcftools, torch), so it runs in a
    # worker thread to keep the event loop responsive for status polling.
    asyncio.create_task(asyncio.to_thread(run_pipeline, run.id))

    return AnalysisStatus(id=run.id, genome_upload_id=run.genome_upload_id, status=run.status)


@router.get("/analysis/{analysis_id}", response_model=AnalysisStatus)
async def get_analysis(
    analysis_id: str, session: AsyncSession = Depends(get_session)
) -> AnalysisStatus:
    run = await session.get(AnalysisRun, analysis_id)
    if run is None:
        raise HTTPException(status_code=404, detail="Analysis not found")

    total = await session.scalar(
        select(func.count(Variant.id)).where(Variant.analysis_run_id == analysis_id)
    )
    known = await session.scalar(
        select(func.count(func.distinct(Variant.id)))
        .join(Annotation, Annotation.variant_id == Variant.id)
        .where(Variant.analysis_run_id == analysis_id)
    )
    predicted_only = await session.scalar(
        select(func.count(func.distinct(Variant.id)))
        .join(AiPrediction, AiPrediction.variant_id == Variant.id)
        .where(
            Variant.analysis_run_id == analysis_id,
            ~Variant.id.in_(select(Annotation.variant_id)),
        )
    )
    diseases = await session.scalar(
        select(func.count(Disease.id)).where(Disease.analysis_run_id == analysis_id)
    )
    transcripts = await session.scalar(
        select(func.count(Transcript.id)).where(Transcript.analysis_run_id == analysis_id)
    )
    asos = await session.scalar(
        select(func.count(Aso.id))
        .join(Transcript, Aso.transcript_id == Transcript.id)
        .where(Transcript.analysis_run_id == analysis_id)
    )

    return AnalysisStatus(
        id=run.id,
        genome_upload_id=run.genome_upload_id,
        status=run.status,
        current_stage=run.current_stage,
        error=run.error,
        counts=AnalysisCounts(
            variants=total or 0,
            known_variants=known or 0,
            unknown_variants=predicted_only or 0,
            diseases=diseases or 0,
            transcripts=transcripts or 0,
            asos=asos or 0,
        ),
    )


@router.get("/analyses", response_model=list[AnalysisStatus])
async def list_analyses(session: AsyncSession = Depends(get_session)) -> list[AnalysisStatus]:
    """Analysis runs, newest first - lets the UI offer an analysis selector."""
    runs = (
        await session.scalars(select(AnalysisRun).order_by(AnalysisRun.created_at.desc()))
    ).all()
    return [
        AnalysisStatus(
            id=run.id,
            genome_upload_id=run.genome_upload_id,
            status=run.status,
            current_stage=run.current_stage,
            error=run.error,
        )
        for run in runs
    ]
