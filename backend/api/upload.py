import asyncio
import os

from fastapi import APIRouter, Depends, HTTPException
from pydantic import BaseModel
from sqlalchemy import update
from sqlalchemy.ext.asyncio import AsyncSession

from backend.database.models import GenomeUpload
from backend.database.session import SessionLocal, get_session
from backend.genomics.validate import validate_fasta

router = APIRouter(prefix="/upload", tags=["upload"])


class UploadRequest(BaseModel):
    file_path: str


class UploadStatus(BaseModel):
    id: str
    status: str
    progress: float
    result: dict | None = None


async def _run_validation(upload_id: str) -> None:
    async with SessionLocal() as session:
        upload = await session.get(GenomeUpload, upload_id)
        upload.status = "validating"
        await session.commit()
        file_path = upload.file_path

    loop = asyncio.get_running_loop()

    async def report(progress: float) -> None:
        async with SessionLocal() as session:
            # One atomic conditional UPDATE, not a read-modify-write: these writes are scheduled
            # from the worker thread and complete in arbitrary order, so a load-then-commit would
            # let a stale value land on top of the final 1.0 and leave the bar stuck partway.
            # Guarding inside the statement makes progress monotonic regardless of arrival order.
            await session.execute(
                update(GenomeUpload)
                .where(GenomeUpload.id == upload_id, GenomeUpload.progress < progress)
                .values(progress=progress)
            )
            await session.commit()

    def on_progress(progress: float) -> None:
        asyncio.run_coroutine_threadsafe(report(progress), loop)

    result = await loop.run_in_executor(None, validate_fasta, file_path, on_progress)

    async with SessionLocal() as session:
        upload = await session.get(GenomeUpload, upload_id)
        upload.status = "valid" if result.valid else "invalid"
        upload.progress = 1.0
        upload.result = result.to_dict()
        await session.commit()


@router.post("", response_model=UploadStatus)
async def create_upload(req: UploadRequest, session: AsyncSession = Depends(get_session)) -> UploadStatus:
    if not os.path.isfile(req.file_path):
        raise HTTPException(status_code=400, detail=f"File not found: {req.file_path}")

    upload = GenomeUpload(file_path=req.file_path, status="pending", progress=0.0)
    session.add(upload)
    await session.commit()
    await session.refresh(upload)

    asyncio.create_task(_run_validation(upload.id))

    return UploadStatus(id=upload.id, status=upload.status, progress=upload.progress)


@router.get("/{upload_id}", response_model=UploadStatus)
async def get_upload(upload_id: str, session: AsyncSession = Depends(get_session)) -> UploadStatus:
    upload = await session.get(GenomeUpload, upload_id)
    if upload is None:
        raise HTTPException(status_code=404, detail="Upload not found")
    return UploadStatus(id=upload.id, status=upload.status, progress=upload.progress, result=upload.result)
