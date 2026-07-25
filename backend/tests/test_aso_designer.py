"""Stage 10 tests: ASO orientation, sequence metrics, design heuristics, off-target degradation.

These never touch the real 3GB reference or shell out to minimap2 - a toy contig stands in for the
genome and the off-target screen is stubbed.
"""

from __future__ import annotations

import random
import subprocess
from pathlib import Path

import pysam
import pytest
from sqlalchemy import select
from sqlalchemy.orm import Session

from backend.aso import designer, offtarget
from backend.aso.designer import (
    ASO_LENGTH,
    ASO_LENGTH_BAND,
    aso_from_genomic,
    design_asos,
    predicted_efficacy,
)
from backend.aso.offtarget import NEUTRAL_OFFTARGET, score_offtarget
from backend.database.models import AnalysisRun, Aso, Transcript

# Toy contig is a plain ATGC repeat, so every window is exactly 50% GC and any 20mer starting at a
# position 1 (mod 4) is "ATGC" * 5 - which makes the metric assertions below exact.
CONTIG = "17"
CONTIG_LENGTH = 3400
KNOWN_TARGET = "ATGC" * 5
KNOWN_ASO = "GCATGCATGCATGCATGCAT"  # reverse complement of KNOWN_TARGET

EXONS = [
    {"start": 1001, "end": 1100, "exon_number": "1"},
    {"start": 2001, "end": 2200, "exon_number": "2"},
    {"start": 3001, "end": 3050, "exon_number": "3"},
]
CANDIDATE_EXON = {
    "exon_number": "2",
    "start": 2001,
    "end": 2200,
    "variant_score": 0.9,
    "variant_positions": [2101],
    "variants": ["17:2101A>G"],
    "mechanisms": ["splice_switch", "exon_skip", "exon_include", "knockdown", "allele_specific"],
}


@pytest.fixture
def toy_reference(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """A non-repeating contig with KNOWN_TARGET planted at the start of exon 2.

    The sequence must not repeat: with a uniform contig every window yields the same 20mer, so
    distinct loci would be indistinguishable and duplicate-collapsing behaviour untestable.
    """
    rng = random.Random(20240724)
    bases = [rng.choice("ACGT") for _ in range(CONTIG_LENGTH)]
    start = 2000  # 0-based; exon 2 begins at 1-based 2001
    bases[start : start + len(KNOWN_TARGET)] = list(KNOWN_TARGET)
    sequence = "".join(bases)

    fasta = tmp_path / "toy.fa"
    fasta.write_text(f">{CONTIG}\n" + "\n".join(
        sequence[i : i + 60] for i in range(0, len(sequence), 60)
    ) + "\n")
    pysam.faidx(str(fasta))
    monkeypatch.setattr(designer.settings, "reference_fasta", fasta)
    return fasta


@pytest.fixture
def stub_offtarget(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(designer, "score_offtarget", lambda sequences: {s: 0.25 for s in sequences})


def _add_transcript(
    session: Session,
    run: AnalysisRun,
    transcript_id: str = "ENST_TOY",
    strand: str = "+",
    opportunities: list[str] | None = None,
    candidate_exons: list[dict] | None = None,
    rank_score: float = 0.8,
    exon_structure: list[dict] | None = None,
) -> Transcript:
    transcript = Transcript(
        analysis_run_id=run.id,
        gene="BRCA1",
        ensembl_transcript_id=transcript_id,
        chrom=CONTIG,
        strand=strand,
        exon_structure=exon_structure if exon_structure is not None else EXONS,
        candidate_exons=candidate_exons if candidate_exons is not None else [CANDIDATE_EXON],
        opportunities=opportunities if opportunities is not None else ["exon_skip"],
        rank_score=rank_score,
    )
    session.add(transcript)
    session.commit()
    return transcript


def _asos(session: Session, run: AnalysisRun) -> list[Aso]:
    return list(
        session.scalars(
            select(Aso)
            .join(Transcript, Transcript.id == Aso.transcript_id)
            .where(Transcript.analysis_run_id == run.id)
        ).all()
    )


# --- orientation: the one detail everything else depends on ---------------------------------


def test_plus_strand_aso_is_the_reverse_complement_of_the_target() -> None:
    assert aso_from_genomic("ATGGCATTAC", "+") == "GTAATGCCAT"


def test_minus_strand_aso_is_the_plus_strand_sequence_itself() -> None:
    # The mRNA of a '-' strand gene is already the reverse complement of the reference plus
    # strand, so reverse-complementing again returns the plus strand.
    assert aso_from_genomic("atggcattac", "-") == "ATGGCATTAC"


def test_designed_plus_strand_aso_is_reverse_complemented(
    session: Session, analysis_run: AnalysisRun, toy_reference: Path, stub_offtarget: None
) -> None:
    _add_transcript(session, analysis_run, opportunities=["exon_skip"])
    design_asos(session, analysis_run.id)

    first = next(aso for aso in _asos(session, analysis_run) if aso.genomic_position.endswith("2001-2020"))
    assert first.sequence == KNOWN_ASO


def test_designed_minus_strand_aso_matches_the_plus_strand(
    session: Session, analysis_run: AnalysisRun, toy_reference: Path, stub_offtarget: None
) -> None:
    _add_transcript(session, analysis_run, strand="-", opportunities=["exon_skip"])
    design_asos(session, analysis_run.id)

    first = next(aso for aso in _asos(session, analysis_run) if aso.genomic_position.endswith("2001-2020"))
    assert first.sequence == KNOWN_TARGET


def test_isoforms_do_not_produce_duplicate_candidates(
    session: Session, analysis_run: AnalysisRun, toy_reference: Path, stub_offtarget: None
) -> None:
    """Isoforms of one gene share exons, so they must not multiply the same reagent."""
    _add_transcript(session, analysis_run, transcript_id="ENST_ISO_A", rank_score=0.9)
    single_isoform_count = design_asos(session, analysis_run.id)

    _add_transcript(session, analysis_run, transcript_id="ENST_ISO_B", rank_score=0.5)
    two_isoform_count = design_asos(session, analysis_run.id)

    assert single_isoform_count > 0
    assert two_isoform_count == single_isoform_count

    asos = _asos(session, analysis_run)
    keys = [(aso.sequence, aso.mechanism) for aso in asos]
    assert len(keys) == len(set(keys))


def test_dedupe_keeps_a_sequence_that_serves_two_mechanisms() -> None:
    """The dedupe key is (sequence, mechanism): one reagent supporting two mechanisms is two
    hypotheses, so only exact sequence+mechanism repeats collapse.
    """
    shared = "ACGTACGTACGTACGTACGT"

    def candidate(mechanism: str, exon: str, variant_backed: bool = True) -> designer._Candidate:
        return designer._Candidate(
            transcript=None,
            mechanism=mechanism,
            target_exon=exon,
            start=1,
            end=20,
            sequence=shared,
            variant_backed=variant_backed,
        )

    kept = designer._dedupe([
        candidate("exon_skip", "2"),
        candidate("knockdown", "2"),
        candidate("exon_skip", "7", variant_backed=False),  # exact seq+mechanism repeat
    ])

    assert [c.mechanism for c in kept] == ["exon_skip", "knockdown"]
    assert kept[0].target_exon == "2"  # first (best-ranked) attribution is the one retained


# --- sequence metrics -----------------------------------------------------------------------


def test_gc_and_tm_for_a_known_window(
    session: Session, analysis_run: AnalysisRun, toy_reference: Path, stub_offtarget: None
) -> None:
    _add_transcript(session, analysis_run, opportunities=["exon_skip"])
    design_asos(session, analysis_run.id)

    known = next(aso for aso in _asos(session, analysis_run) if aso.sequence == KNOWN_ASO)
    assert known.gc_pct == 50.0
    # Nearest-neighbour Tm. A GC-content approximation would give ~50.4 here, so this pins Tm_NN.
    assert known.tm == pytest.approx(55.5, abs=0.05)


def test_designed_asos_sit_in_the_therapeutic_length_band(
    session: Session, analysis_run: AnalysisRun, toy_reference: Path, stub_offtarget: None
) -> None:
    _add_transcript(
        session,
        analysis_run,
        opportunities=["splice_switch", "exon_skip", "exon_include", "knockdown", "allele_specific"],
    )
    design_asos(session, analysis_run.id)

    low, high = ASO_LENGTH_BAND
    assert low <= ASO_LENGTH <= high
    asos = _asos(session, analysis_run)
    assert asos
    assert all(low <= len(aso.sequence) <= high for aso in asos)


# --- design heuristics ----------------------------------------------------------------------


def test_poly_g_scores_worse_than_a_balanced_sequence() -> None:
    balanced = "ACTGACTGCATCAGTCAGTA"
    quadruplex = "ACTGGGGGCATCAGTCAGTA"
    assert predicted_efficacy(quadruplex, 55.0, 55.0) < predicted_efficacy(balanced, 55.0, 55.0)


def test_homopolymer_run_scores_worse_than_a_balanced_sequence() -> None:
    balanced = "ACTGACTGCATCAGTCAGTA"
    low_complexity = "ACTGAAAAACATCAGTCAGTA"
    assert predicted_efficacy(low_complexity, 55.0, 55.0) < predicted_efficacy(balanced, 55.0, 55.0)


def test_extreme_gc_and_tm_score_worse_than_the_favourable_band() -> None:
    sequence = "ACTGACTGCATCAGTCAGTA"
    assert predicted_efficacy(sequence, 50.0, 55.0) == 1.0
    assert predicted_efficacy(sequence, 90.0, 80.0) < predicted_efficacy(sequence, 50.0, 55.0)


# --- row contents ---------------------------------------------------------------------------


def test_rows_are_traceable_and_use_the_transcript_opportunities(
    session: Session, analysis_run: AnalysisRun, toy_reference: Path, stub_offtarget: None
) -> None:
    _add_transcript(session, analysis_run, opportunities=["splice_switch", "knockdown"])
    design_asos(session, analysis_run.id)

    asos = _asos(session, analysis_run)
    assert {aso.mechanism for aso in asos} == {"splice_switch", "knockdown"}
    for aso in asos:
        assert aso.target_exon == "2"
        assert aso.genomic_position.startswith("chr17:")
        # predicted_specificity is the inverse of the off-target risk.
        assert aso.predicted_specificity == pytest.approx(1.0 - aso.off_target_score)
        assert 0.0 <= aso.confidence <= 1.0


def test_generic_knockdown_fallback_has_lower_confidence(
    session: Session, analysis_run: AnalysisRun, toy_reference: Path, stub_offtarget: None
) -> None:
    """The two transcripts must sit at different loci: sharing an exon would (correctly) collapse
    their identical sequences to one candidate, leaving nothing to compare.
    """
    backed = _add_transcript(session, analysis_run, "ENST_BACKED", opportunities=["knockdown"])
    _add_transcript(
        session,
        analysis_run,
        "ENST_GENERIC",
        opportunities=["knockdown"],
        candidate_exons=[],
        exon_structure=[{"start": 1001, "end": 1100, "exon_number": "1"}],
    )
    design_asos(session, analysis_run.id)

    by_transcript = {aso.transcript_id: aso for aso in _asos(session, analysis_run)}
    generic = next(aso for tid, aso in by_transcript.items() if tid != backed.id)
    assert generic.confidence < by_transcript[backed.id].confidence


def test_no_transcripts_designs_nothing(
    session: Session, analysis_run: AnalysisRun, toy_reference: Path, stub_offtarget: None
) -> None:
    assert design_asos(session, analysis_run.id) == 0


def test_design_asos_is_idempotent(
    session: Session, analysis_run: AnalysisRun, toy_reference: Path, stub_offtarget: None
) -> None:
    _add_transcript(session, analysis_run, opportunities=["exon_skip", "knockdown"])

    first = design_asos(session, analysis_run.id)
    second = design_asos(session, analysis_run.id)

    assert first == second > 0
    assert len(_asos(session, analysis_run)) == first


# --- off-target screen ----------------------------------------------------------------------


def test_score_offtarget_is_neutral_when_minimap2_is_unavailable(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def missing_binary(*args: object, **kwargs: object) -> None:
        raise FileNotFoundError("minimap2")

    monkeypatch.setattr(offtarget.subprocess, "run", missing_binary)
    assert score_offtarget(["ACGT", "TTTT"]) == {
        "ACGT": NEUTRAL_OFFTARGET,
        "TTTT": NEUTRAL_OFFTARGET,
    }


def test_score_offtarget_is_neutral_when_minimap2_fails(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        offtarget.subprocess,
        "run",
        lambda *a, **k: subprocess.CompletedProcess(a, 1, "", "index not found"),
    )
    assert score_offtarget(["ACGT"]) == {"ACGT": NEUTRAL_OFFTARGET}


def _paf_line(query: str, target: str, matching_bases: int = 20) -> str:
    """A minimal PAF row: field 6 is the target name, field 10 the matching base count."""
    return f"{query}\t20\t0\t20\t+\t{target}\t100\t1\t20\t{matching_bases}\t20\t60"


def _stub_paf(monkeypatch: pytest.MonkeyPatch, paf: str, gene_map: dict[str, str]) -> None:
    monkeypatch.setattr(
        offtarget.subprocess, "run", lambda *a, **k: subprocess.CompletedProcess(a, 0, paf, "")
    )
    monkeypatch.setattr(offtarget, "_transcript_to_gene", lambda: gene_map)


def test_score_offtarget_counts_distinct_genes(monkeypatch: pytest.MonkeyPatch) -> None:
    # Sequences are written to the query FASTA sorted, so these become s0, s1, s2.
    unique, moderate, promiscuous = "AAAC", "AAAG", "AAAT"
    gene_map = {f"T{i}": f"GENE{i}" for i in range(12)}

    paf = "\n".join(
        [_paf_line("s0", "T0")]
        + [_paf_line("s1", f"T{i}") for i in range(6)]  # 5 off-target genes
        + [_paf_line("s2", f"T{i}") for i in range(12)]  # saturates the cap
    )
    _stub_paf(monkeypatch, paf, gene_map)

    scores = score_offtarget([promiscuous, unique, moderate])
    assert scores[unique] == 0.0
    assert scores[moderate] == pytest.approx(0.5)
    assert scores[promiscuous] == 1.0


def test_same_gene_isoforms_are_on_target(monkeypatch: pytest.MonkeyPatch) -> None:
    """A gene with many isoforms must not look promiscuous - every hit collapses to one gene."""
    gene_map = {f"BRCA1_iso{i}": "BRCA1" for i in range(40)}
    paf = "\n".join(_paf_line("s0", f"BRCA1_iso{i}") for i in range(40))
    _stub_paf(monkeypatch, paf, gene_map)

    assert score_offtarget(["ACGT"]) == {"ACGT": 0.0}


def test_short_seed_matches_are_ignored(monkeypatch: pytest.MonkeyPatch) -> None:
    """Loose seeding yields many 12-13bp hits that are too short to hybridize; they don't count."""
    gene_map = {f"T{i}": f"GENE{i}" for i in range(30)}
    paf = "\n".join(
        [_paf_line("s0", "T0", matching_bases=20)]
        + [_paf_line("s0", f"T{i}", matching_bases=12) for i in range(1, 30)]
    )
    _stub_paf(monkeypatch, paf, gene_map)

    assert score_offtarget(["ACGT"]) == {"ACGT": 0.0}


def test_unmapped_target_falls_back_to_transcript_id(monkeypatch: pytest.MonkeyPatch) -> None:
    """A missing gene map must not silently report zero risk."""
    paf = "\n".join(_paf_line("s0", f"T{i}") for i in range(6))
    _stub_paf(monkeypatch, paf, {})

    assert score_offtarget(["ACGT"]) == {"ACGT": pytest.approx(0.5)}


def test_score_offtarget_handles_no_sequences() -> None:
    assert score_offtarget([]) == {}
