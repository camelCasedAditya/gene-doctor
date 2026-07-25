"""Off-target screen for candidate ASO sequences, via minimap2 against the cDNA transcriptome.

An ASO hybridizes to RNA, so the relevant question is how many *other genes* a candidate could
bind - not how many places it occurs in the genome. Two consequences shape the scoring here:

  * Multiple hits to different isoforms of the same gene are on-target, so hits are collapsed to
    distinct gene symbols before counting.
  * A 20nt query cannot be seeded by the genome index's k=19 (it returns no alignments at all),
    so this runs against a dedicated k=12/w=1 transcriptome index.

Every candidate goes into one FASTA and one minimap2 invocation - the index takes seconds to load,
so per-sequence calls are not an option.
"""

from __future__ import annotations

import logging
import subprocess
import tempfile
from functools import lru_cache
from pathlib import Path

from backend.settings import settings

logger = logging.getLogger(__name__)

# Returned when minimap2 can't tell us anything - unresolved is not the same as proven unique.
NEUTRAL_OFFTARGET = 0.5

# Loose seeding surfaces thousands of 12-13bp seed hits that are far too short to hybridize
# stably. Only near-full-length matches are plausible off-target binders.
MIN_MATCH_BASES = 18

# Distinct off-target genes at which risk saturates at 1.0.
OFFTARGET_GENE_CAP = 10
MAX_HITS_REPORTED = 100
MINIMAP2_TIMEOUT = 900


@lru_cache(maxsize=1)
def _transcript_to_gene() -> dict[str, str]:
    """transcript_id -> gene_symbol, as written by the transcriptome ingest step."""
    try:
        with open(settings.transcript_gene_map) as f:
            return dict(
                line.rstrip("\n").split("\t", 1) for line in f if "\t" in line
            )
    except OSError:
        logger.warning("transcript->gene map unavailable; off-target genes cannot be resolved")
        return {}


def score_offtarget(sequences: list[str]) -> dict[str, float]:
    """Map each ASO sequence to an off-target score in [0,1] (higher = more off-target risk)."""
    if not sequences:
        return {}

    unique = sorted(set(sequences))
    with tempfile.TemporaryDirectory() as tmp_dir:
        query = Path(tmp_dir) / "aso_candidates.fa"
        query.write_text("".join(f">s{i}\n{seq}\n" for i, seq in enumerate(unique)))
        paf = _run_minimap2(query)

    if paf is None:
        return {seq: NEUTRAL_OFFTARGET for seq in unique}

    genes_per_query = _genes_per_query(paf)
    return {seq: _risk(genes_per_query.get(f"s{i}")) for i, seq in enumerate(unique)}


def _genes_per_query(paf: str) -> dict[str, set[str]]:
    """Distinct gene symbols hit by each query, counting only near-full-length matches."""
    mapping = _transcript_to_gene()
    hits: dict[str, set[str]] = {}

    for line in paf.splitlines():
        fields = line.split("\t")
        if len(fields) < 10:
            continue
        try:
            matching_bases = int(fields[9])
        except ValueError:
            continue
        if matching_bases < MIN_MATCH_BASES:
            continue
        # Fall back to the transcript id when the gene map is missing, so an unresolvable map
        # degrades to per-transcript counting rather than silently reporting zero risk.
        target = fields[5]
        hits.setdefault(fields[0], set()).add(mapping.get(target, target))

    return hits


def _risk(genes: set[str] | None) -> float:
    """One gene is the intended target (risk 0); every additional gene is an off-target risk."""
    if not genes:
        return NEUTRAL_OFFTARGET
    return round(min((len(genes) - 1) / OFFTARGET_GENE_CAP, 1.0), 4)


def _run_minimap2(query: Path) -> str | None:
    """Run the batch alignment. Returns PAF text, or None if the screen couldn't be performed.

    The sr preset suits short queries, but its seeding thresholds (-n2 -m20 -s40) reject 20nt reads
    outright, so they are loosened; -p0/-N/--secondary=yes keep every repeat copy in the output
    instead of collapsing them, since repeat copies are exactly what we are counting.
    """
    cmd = [
        "minimap2",
        "-x", "sr",
        "-n", "1",
        "-m", "10",
        "-s", "20",
        "-p", "0",
        "-N", str(MAX_HITS_REPORTED),
        "--secondary=yes",
        "-t", "4",
        str(settings.transcriptome_mmi),
        str(query),
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=MINIMAP2_TIMEOUT)
    except (OSError, subprocess.TimeoutExpired) as exc:
        logger.warning("off-target screen unavailable, using neutral scores: %s", exc)
        return None

    if result.returncode != 0:
        logger.warning(
            "minimap2 exited %d, using neutral off-target scores: %s",
            result.returncode,
            result.stderr[-500:],
        )
        return None
    return result.stdout
