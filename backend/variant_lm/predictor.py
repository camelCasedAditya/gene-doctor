"""AI prediction for variants that lack validated clinical evidence.

Two independent signals, each written as its own AiPrediction row:

  - "variant_lm": Nucleotide Transformer ref-vs-alt log-likelihood delta over a genomic
    window, used as a proxy for functional/regulatory disruption. Works anywhere in the
    genome (coding or not).
  - "alphamissense": precomputed missense pathogenicity, looked up by tabix. Only coding
    missense variants have an entry.

Everything here is an AI-generated research hypothesis and lands exclusively in
AiPrediction - never in the Annotation table, which holds validated clinical evidence only.

Model loading and inference are best-effort: if weights can't be fetched or inference dies,
the stage logs a warning and returns the predictions it did manage to make.
"""

from __future__ import annotations

import logging
import math
from functools import lru_cache
from pathlib import Path
from typing import Any

import pysam
from sqlalchemy import select
from sqlalchemy.orm import Session

from backend.database.models import AiPrediction, Annotation, Variant
from backend.settings import settings

logger = logging.getLogger(__name__)

# A variant is "already explained" once any single annotation reaches this evidence strength;
# below that it goes to the AI predictor.
WEAK_EVIDENCE_THRESHOLD: float = 0.5

# Nucleotide Transformer v1, trained on the human reference - exactly the distribution we are
# measuring deviation from. The smaller v2 checkpoints would be cheaper, but they ship their
# architecture as remote code that imports symbols transformers 5.x removed, so they cannot be
# loaded here at all; this one is a plain EsmForMaskedLM. Swap freely for any masked DNA LM.
NT_MODEL_ID = "InstaDeepAI/nucleotide-transformer-500m-human-ref"

# Genomic context handed to the model, centred on the variant. NT v2 handles up to 12 kb;
# 1 kb is enough to cover local regulatory context while keeping inference cheap.
WINDOW_BP = 1000

# Calibration knob for the [0, 1] normalisation below: a per-token mean log-likelihood drop of
# this size maps to an impact of 1 - 1/e ≈ 0.63. Measured on real GRCh38 SNVs with the model
# below, one changed 6-mer token moves the window mean by 4e-4 to 3e-2, so 1e-2 puts a typical
# disruptive substitution in the upper half of the range. Retune if the checkpoint changes.
DELTA_SCALE = 0.01

# AlphaMissense is a published, benchmarked model, so its scores get a higher confidence than
# our own log-likelihood proxy.
ALPHAMISSENSE_CONFIDENCE = 0.9


@lru_cache(maxsize=1)
def _load_lm() -> tuple[Any, Any, str] | None:
    """Lazily load tokenizer + model, cached for the process. None if unavailable.

    ponytail: the failure is cached too, so one bad load (no network) disables the LM for the
    whole process rather than re-downloading per variant. Restart to retry.
    """
    try:
        import torch
        from transformers import AutoModelForMaskedLM, AutoTokenizer

        device = "mps" if torch.backends.mps.is_available() else "cpu"
        tokenizer = AutoTokenizer.from_pretrained(NT_MODEL_ID)
        model = AutoModelForMaskedLM.from_pretrained(NT_MODEL_ID).to(device).eval()
        logger.info("Loaded %s on %s", NT_MODEL_ID, device)
        return tokenizer, model, device
    except Exception:
        logger.warning("Could not load %s; skipping LM predictions", NT_MODEL_ID, exc_info=True)
        return None


def _open(opener: Any, path: Path) -> Any | None:
    """Open a pysam reference file, or None if it is missing/unindexed."""
    try:
        return opener(str(path))
    except (OSError, ValueError):
        logger.warning("Reference file unavailable: %s", path, exc_info=True)
        return None


def predict_variants(session: Session, analysis_run_id: str) -> int:
    """Score variants lacking strong validated evidence. Creates AiPrediction rows.
    Returns count of predictions created.
    """
    strong_evidence = select(Annotation.variant_id).where(
        Annotation.evidence_strength >= WEAK_EVIDENCE_THRESHOLD
    )
    variants = session.scalars(
        select(Variant).where(
            Variant.analysis_run_id == analysis_run_id,
            Variant.id.not_in(strong_evidence),
        )
    ).all()
    if not variants:
        return 0

    fasta = _open(pysam.FastaFile, settings.reference_fasta)
    alphamissense = _open(pysam.TabixFile, settings.alphamissense_tsv)
    lm = _load_lm() if fasta is not None else None

    created = 0
    for variant in variants:
        if lm is not None:
            try:
                prediction = _lm_prediction(variant, fasta, lm)
            except Exception:
                logger.warning(
                    "LM inference failed on %s:%s; skipping remaining LM predictions",
                    variant.chrom,
                    variant.pos,
                    exc_info=True,
                )
                lm = None  # one failure (OOM, bad tensor) will repeat - stop trying
                prediction = None
            if prediction is not None:
                session.add(prediction)
                created += 1

        if alphamissense is not None:
            prediction = _alphamissense_prediction(variant, alphamissense)
            if prediction is not None:
                session.add(prediction)
                created += 1

    for handle in (fasta, alphamissense):
        if handle is not None:
            handle.close()

    session.commit()
    return created


def _lm_prediction(variant: Variant, fasta: pysam.FastaFile, lm: tuple) -> AiPrediction | None:
    """Nucleotide Transformer ref-vs-alt log-likelihood delta for one variant."""
    window = _reference_window(variant, fasta)
    if window is None:
        return None
    ref_seq, offset = window

    if ref_seq[offset : offset + len(variant.ref)] != variant.ref.upper():
        logger.warning(
            "Reference mismatch at %s:%s (expected %s); skipping LM prediction",
            variant.chrom,
            variant.pos,
            variant.ref,
        )
        return None
    alt_seq = ref_seq[:offset] + variant.alt.upper() + ref_seq[offset + len(variant.ref) :]

    delta = _sequence_log_likelihood(alt_seq, lm) - _sequence_log_likelihood(ref_seq, lm)
    impact = _normalize_delta(delta)

    # The window score can't tell a coding disruption from a regulatory one, so the same
    # normalised delta is reported under both names rather than inventing a split. Confidence
    # is the share of unambiguous bases in the window - an N-heavy window means poor context.
    known_bases = sum(ref_seq.count(base) for base in "ACGT") / len(ref_seq)
    return AiPrediction(
        variant_id=variant.id,
        source="variant_lm",
        functional_impact=impact,
        regulatory_impact=impact,
        confidence=round(known_bases, 4),
    )


def _reference_window(variant: Variant, fasta: pysam.FastaFile) -> tuple[str, int] | None:
    """WINDOW_BP of reference centred on the variant, plus the variant's offset within it."""
    start = max(variant.pos - 1 - WINDOW_BP // 2, 0)
    try:
        sequence = fasta.fetch(variant.chrom, start, start + WINDOW_BP).upper()
    except (KeyError, ValueError):
        logger.warning("Contig %s absent from reference; skipping", variant.chrom)
        return None
    offset = variant.pos - 1 - start
    if offset + len(variant.ref) > len(sequence):
        return None  # variant runs off the end of the contig
    return sequence, offset


def _sequence_log_likelihood(sequence: str, lm: tuple) -> float:
    """Mean per-token log-likelihood of a sequence under the masked LM.

    ponytail: a single unmasked forward pass, not a per-token masked pseudo-likelihood - one
    pass instead of ~170, at the cost of letting the model see the base it is scoring. Upgrade
    to masked-marginal scoring of the variant token if the deltas prove too flat. Per-token
    mean (not sum) keeps ref and alt comparable when an indel changes the token count.
    """
    import torch

    tokenizer, model, device = lm
    encoded = tokenizer(sequence, return_tensors="pt").to(device)
    with torch.no_grad():
        logits = model(**encoded).logits
    log_probs = torch.log_softmax(logits.float(), dim=-1)
    token_ids = encoded["input_ids"]
    return log_probs.gather(-1, token_ids.unsqueeze(-1)).squeeze(-1).mean().item()


def _normalize_delta(delta: float) -> float:
    """Map a ref→alt log-likelihood delta onto [0, 1].

    Only drops count: alt sequence less likely than ref (delta < 0) means the model finds the
    alt allele unexpected in its genomic context, which is the disruption signal. An alt that
    looks *more* like normal genome scores 0.
    """
    return 1.0 - math.exp(-max(0.0, -delta) / DELTA_SCALE)


def _alphamissense_prediction(
    variant: Variant, alphamissense: pysam.TabixFile
) -> AiPrediction | None:
    """Precomputed AlphaMissense pathogenicity for one variant, if it has an entry.

    Most variants legitimately have none (the file covers missense substitutions only). Where
    several transcripts score the same substitution, the worst case is kept.
    """
    try:
        rows = alphamissense.fetch(f"chr{variant.chrom}", variant.pos - 1, variant.pos)
    except (KeyError, ValueError):
        return None

    scores = [
        float(fields[8])
        for fields in (row.split("\t") for row in rows)
        if len(fields) > 8
        and int(fields[1]) == variant.pos
        and fields[2] == variant.ref
        and fields[3] == variant.alt
    ]
    if not scores:
        return None

    return AiPrediction(
        variant_id=variant.id,
        source="alphamissense",
        pathogenicity_probability=max(scores),
        confidence=ALPHAMISSENSE_CONFIDENCE,
    )
