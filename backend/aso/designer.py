"""Stage 10: design candidate ASO sequences against the run's top-ranked transcripts.

An ASO hybridises to the mature transcript, so its sequence is the reverse complement of the
mRNA. The reference FASTA only stores the plus strand, which means the direction depends on the
gene's strand - see aso_from_genomic, that flip is the single most important detail in this file.

Target windows are chosen per mechanism from the transcript's candidate_exons, then scored with
published rules of thumb (GC content, nearest-neighbour Tm, sequence-motif liabilities). Those
scores are research heuristics for prioritising which candidates to look at first - they are not
validated efficacy predictions, and nothing here has been tested in a cell.
"""

from __future__ import annotations

import re
from collections.abc import Iterable
from dataclasses import dataclass

import pysam
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp, gc_fraction
from sqlalchemy import delete, select
from sqlalchemy.orm import Session

from backend.aso.offtarget import score_offtarget
from backend.database.models import Aso, Transcript
from backend.settings import settings

# Therapeutic ASOs run 18-25nt; 20 is the usual starting length.
ASO_LENGTH = 20
ASO_LENGTH_BAND = (18, 25)
WINDOW_STEP = 5
MAX_ASOS_PER_MECHANISM = 3

# Design windows, in bp.
SPLICE_FLANK = 15  # either side of an exon edge, for splice-modulating ASOs
EXON_BODY_SPAN = 40  # stretch of exon body sampled for RNase-H / enhancer targeting
INTRON_OFFSET = 6  # skip the donor consensus before targeting an intronic silencer
ALLELE_SPAN = 25  # centred on the variant, so every window covers the discriminating base

# Favourable ranges from standard ASO design practice.
GC_BAND = (40.0, 60.0)
TM_BAND = (45.0, 65.0)
GC_TOLERANCE = 20.0
TM_TOLERANCE = 15.0
G_QUADRUPLEX_PENALTY = 0.5
HOMOPOLYMER_PENALTY = 0.2


@dataclass
class _Candidate:
    transcript: Transcript
    mechanism: str
    target_exon: str | None
    start: int  # 1-based inclusive genomic start of the target window
    end: int
    sequence: str  # the ASO itself, already oriented against the transcript
    variant_backed: bool


def design_asos(session: Session, analysis_run_id: str, top_n_transcripts: int = 20) -> int:
    """Generate candidate Aso rows for the run's ranked transcripts. Returns count created."""
    session.execute(
        delete(Aso).where(
            Aso.transcript_id.in_(
                select(Transcript.id).where(Transcript.analysis_run_id == analysis_run_id)
            )
        )
    )
    session.commit()

    transcripts = session.scalars(
        select(Transcript)
        .where(Transcript.analysis_run_id == analysis_run_id)
        .order_by(Transcript.rank_score.desc())
        .limit(top_n_transcripts)
    ).all()
    if not transcripts:
        return 0

    fasta = pysam.FastaFile(str(settings.reference_fasta))
    try:
        candidates = _dedupe(
            candidate
            for transcript in transcripts
            for candidate in _transcript_candidates(fasta, transcript)
        )
    finally:
        fasta.close()
    if not candidates:
        return 0

    off_target = score_offtarget([candidate.sequence for candidate in candidates])
    for candidate in candidates:
        session.add(_build_aso(candidate, off_target[candidate.sequence]))
    session.commit()
    return len(candidates)


def _dedupe(candidates: Iterable[_Candidate]) -> list[_Candidate]:
    """Collapse candidates that repeat the same sequence for the same mechanism.

    A gene's isoforms share most of their exons, so designing across 40 BRCA1 isoforms yields the
    same 20mer many times over under a different exon label - noise, not 40 distinct reagents.
    Deduping on (sequence, mechanism) keeps a sequence that genuinely supports two different
    mechanisms, since those are two different therapeutic hypotheses for one reagent. Candidates
    arrive in transcript-rank order, so the retained copy is attributed to the best-ranked
    transcript.
    """
    seen: set[tuple[str, str]] = set()
    unique: list[_Candidate] = []
    for candidate in candidates:
        key = (candidate.sequence, candidate.mechanism)
        if key in seen:
            continue
        seen.add(key)
        unique.append(candidate)
    return unique


def aso_from_genomic(genomic_sequence: str, strand: str) -> str:
    """Convert plus-strand genomic sequence into the ASO that hybridises to that transcript.

    A '+' strand transcript's mRNA *is* the plus-strand sequence, so the ASO is its reverse
    complement. A '-' strand transcript's mRNA is already the reverse complement of the plus
    strand, so the ASO comes back out as the plus-strand sequence itself.
    """
    sequence = genomic_sequence.upper()
    if strand == "+":
        return str(Seq(sequence).reverse_complement())
    return sequence


def predicted_efficacy(sequence: str, gc_pct: float, tm: float) -> float:
    """Heuristic prioritisation score in [0, 1] - NOT a validated efficacy prediction.

        0.60 * gc_term + 0.40 * tm_term - motif penalties

    - gc_term: 1.0 inside 40-60% GC (the workable duplex-stability window), decaying to 0 over
      20 percentage points outside it. Low GC binds too weakly, high GC too promiscuously.
    - tm_term: 1.0 inside 45-65 C by nearest-neighbour Tm, decaying over 15 C outside it.
    - runs of >=4 G: heavy penalty, they form G-quadruplexes and aggregate rather than hybridise.
    - any other run of >=4 identical bases: smaller penalty for low-complexity stretches, which
      hybridise out of register.
    """
    score = 0.60 * _band_term(gc_pct, GC_BAND, GC_TOLERANCE) + 0.40 * _band_term(
        tm, TM_BAND, TM_TOLERANCE
    )
    if re.search(r"G{4,}", sequence):
        score -= G_QUADRUPLEX_PENALTY
    if re.search(r"([ACT])\1{3,}", sequence):
        score -= HOMOPOLYMER_PENALTY
    return round(_clamp(score), 4)


def _band_term(value: float, band: tuple[float, float], tolerance: float) -> float:
    low, high = band
    if low <= value <= high:
        return 1.0
    distance = low - value if value < low else value - high
    return max(0.0, 1.0 - distance / tolerance)


def _clamp(value: float) -> float:
    return min(max(value, 0.0), 1.0)


def _transcript_candidates(fasta: pysam.FastaFile, transcript: Transcript) -> list[_Candidate]:
    candidates: list[_Candidate] = []
    for mechanism in transcript.opportunities:
        made = 0
        for target_exon, region_start, region_end, variant_backed in _target_regions(
            transcript, mechanism
        ):
            for start in range(region_start, region_end - ASO_LENGTH + 2, WINDOW_STEP):
                if made >= MAX_ASOS_PER_MECHANISM:
                    break
                end = start + ASO_LENGTH - 1
                genomic = _fetch(fasta, transcript.chrom, start, end)
                if genomic is None:
                    continue
                candidates.append(
                    _Candidate(
                        transcript=transcript,
                        mechanism=mechanism,
                        target_exon=target_exon,
                        start=start,
                        end=end,
                        sequence=aso_from_genomic(genomic, transcript.strand),
                        variant_backed=variant_backed,
                    )
                )
                made += 1
            if made >= MAX_ASOS_PER_MECHANISM:
                break
    return candidates


def _target_regions(
    transcript: Transcript, mechanism: str
) -> list[tuple[str | None, int, int, bool]]:
    """Genomic windows (1-based inclusive) worth tiling for one mechanism.

    Boundary mechanisms use the highest-scoring candidate exon; knockdown falls back to the
    longest exon when no variant landed in this isoform, since RNase-H cleavage only needs
    accessible exon body. Returns (exon label, start, end, backed-by-a-variant).
    """
    exons: list[dict] = transcript.candidate_exons or []
    best = max(exons, key=lambda exon: exon["variant_score"], default=None)

    if mechanism == "knockdown":
        exon = best or _longest_exon(transcript.exon_structure)
        if exon is None:
            return []
        centre = (int(exon["start"]) + int(exon["end"])) // 2
        return [(
            str(exon.get("exon_number", "")) or None,
            centre - EXON_BODY_SPAN // 2,
            centre + EXON_BODY_SPAN // 2,
            best is not None,
        )]

    if best is None:
        return []
    label = str(best.get("exon_number", "")) or None
    start, end = int(best["start"]), int(best["end"])

    if mechanism == "splice_switch":
        # Straddle both exon edges - the acceptor and donor decisions are both addressable.
        return [
            (label, edge - SPLICE_FLANK, edge + SPLICE_FLANK, True) for edge in (start, end)
        ]
    if mechanism == "exon_skip":
        # Exonic splicing enhancers cluster near the exon edge; stay inside the exon.
        return [(label, start, min(start + EXON_BODY_SPAN - 1, end), True)]
    if mechanism == "exon_include":
        # Intronic splicing silencer just past the donor site (the ISS-N1/nusinersen pattern).
        return [(label, end + INTRON_OFFSET, end + INTRON_OFFSET + EXON_BODY_SPAN - 1, True)]
    if mechanism == "allele_specific":
        # Every window must cover the discriminating base, so centre tightly on the variant.
        position = (best.get("variant_positions") or [(start + end) // 2])[0]
        return [(label, position - ALLELE_SPAN // 2, position + ALLELE_SPAN // 2, True)]
    return []


def _longest_exon(exon_structure: list[dict]) -> dict | None:
    if not exon_structure:
        return None
    return max(exon_structure, key=lambda exon: int(exon["end"]) - int(exon["start"]))


def _fetch(fasta: pysam.FastaFile, chrom: str, start: int, end: int) -> str | None:
    """Plus-strand reference sequence for a 1-based inclusive span, or None if unusable."""
    try:
        sequence = fasta.fetch(chrom, start - 1, end).upper()
    except (KeyError, ValueError):
        return None
    if len(sequence) != end - start + 1 or set(sequence) - set("ACGT"):
        return None  # off the end of the contig, or an assembly gap
    return sequence


def _build_aso(candidate: _Candidate, off_target_score: float) -> Aso:
    gc_pct = gc_fraction(candidate.sequence) * 100
    tm = MeltingTemp.Tm_NN(candidate.sequence)
    # Confidence tracks the evidence chain, not the chemistry: how strongly the transcript was
    # implicated, plus whether this specific window was picked out by a variant or is a generic
    # knockdown fallback. Research confidence - it says nothing about clinical validity.
    evidence = 1.0 if candidate.variant_backed else 0.4
    return Aso(
        transcript_id=candidate.transcript.id,
        sequence=candidate.sequence,
        genomic_position=f"chr{candidate.transcript.chrom}:{candidate.start}-{candidate.end}",
        target_exon=candidate.target_exon,
        mechanism=candidate.mechanism,
        gc_pct=round(gc_pct, 2),
        tm=round(tm, 2),
        predicted_efficacy=predicted_efficacy(candidate.sequence, gc_pct, tm),
        predicted_specificity=round(1.0 - off_target_score, 4),
        off_target_score=off_target_score,
        confidence=round(_clamp(0.7 * (candidate.transcript.rank_score or 0.0) + 0.3 * evidence), 4),
    )
