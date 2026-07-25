"""Stage 6 tests. No network and no model weights: the LM loader and the per-sequence
log-likelihood are stubbed, while window extraction, normalisation, the AlphaMissense tabix
lookup and the evidence filter all run for real.
"""

from __future__ import annotations

import math
import random

import pysam
import pytest
from sqlalchemy import func, select
from sqlalchemy.orm import Session

from backend.database.models import AiPrediction, AnalysisRun, Annotation
from backend.settings import settings
from backend.tests.conftest import add_annotation, make_variant
from backend.variant_lm import predictor
from backend.variant_lm.predictor import DELTA_SCALE, predict_variants

# One AlphaMissense row per line, tab-separated, as in AlphaMissense_hg38.tsv.gz.
ALPHAMISSENSE_ROWS = [
    ("chr1", 1000, "A", "G", "Q8NH21", "ENST00000335137.4", "V2L", "0.2937", "likely_benign"),
    ("chr1", 2000, "C", "T", "Q8NH21", "ENST00000335137.4", "P7S", "0.4100", "ambiguous"),
    ("chr1", 2000, "C", "T", "Q8NH21", "ENST00000999999.1", "P9S", "0.8123", "likely_pathogenic"),
]


def _write_reference(tmp_path) -> str:
    """A single 3 kb contig "1" with known bases at the positions the tests use."""
    rng = random.Random(0)
    bases = [rng.choice("ACGT") for _ in range(3000)]
    bases[999] = "A"  # 1-based pos 1000
    bases[1999] = "C"  # 1-based pos 2000
    sequence = "".join(bases)
    path = tmp_path / "ref.fa"
    path.write_text(">1\n" + "\n".join(sequence[i : i + 60] for i in range(0, 3000, 60)) + "\n")
    pysam.faidx(str(path))
    return str(path)


def _write_alphamissense(tmp_path) -> str:
    """A bgzip+tabix AlphaMissense fixture, comment header and chr-prefixed contigs included."""
    path = tmp_path / "am.tsv"
    header = "#copyright notice\n#CHROM\tPOS\tREF\tALT\tgenome\tuniprot_id\ttranscript_id\tprotein_variant\tam_pathogenicity\tam_class\n"
    body = "".join(
        "\t".join([chrom, str(pos), ref, alt, "hg38", uniprot, tx, pv, score, cls]) + "\n"
        for chrom, pos, ref, alt, uniprot, tx, pv, score, cls in ALPHAMISSENSE_ROWS
    )
    path.write_text(header + body)
    pysam.tabix_compress(str(path), str(path) + ".gz", force=True)
    pysam.tabix_index(str(path) + ".gz", seq_col=0, start_col=1, end_col=1, zerobased=False, force=True)
    return str(path) + ".gz"


@pytest.fixture
def refs(tmp_path, monkeypatch):
    """Point settings at the tiny fixtures and disable the LM unless a test asks for it."""
    monkeypatch.setattr(settings, "reference_fasta", _write_reference(tmp_path))
    monkeypatch.setattr(settings, "alphamissense_tsv", _write_alphamissense(tmp_path))
    monkeypatch.setattr(predictor, "_load_lm", lambda: None)


def _fake_log_likelihood(sequence: str, lm: tuple) -> float:
    """Deterministic stand-in: every G costs 0.01, so an A>G substitution makes the alt
    sequence exactly 0.01 less likely than the ref regardless of call order.
    """
    return -0.5 - 0.01 * sequence.count("G")


@pytest.fixture
def fake_lm(refs, monkeypatch):
    monkeypatch.setattr(predictor, "_load_lm", lambda: ("tokenizer", "model", "cpu"))
    monkeypatch.setattr(predictor, "_sequence_log_likelihood", _fake_log_likelihood)


def _predictions(session: Session, source: str | None = None) -> list[AiPrediction]:
    query = select(AiPrediction)
    if source is not None:
        query = query.where(AiPrediction.source == source)
    return list(session.scalars(query))


def test_strong_evidence_variant_is_skipped(session, analysis_run: AnalysisRun, fake_lm):
    variant = make_variant(session, analysis_run, chrom="1", pos=1000, ref="A", alt="G")
    add_annotation(session, variant, evidence_strength=1.0)

    assert predict_variants(session, analysis_run.id) == 0
    assert _predictions(session) == []


def test_weak_evidence_variant_is_predicted(session, analysis_run: AnalysisRun, fake_lm):
    variant = make_variant(session, analysis_run, chrom="1", pos=1000, ref="A", alt="G")
    add_annotation(session, variant, evidence_strength=0.2)

    created = predict_variants(session, analysis_run.id)

    assert created == 2  # one variant_lm row + one alphamissense row
    lm_rows = _predictions(session, "variant_lm")
    assert len(lm_rows) == 1
    expected = 1.0 - math.exp(-0.01 / DELTA_SCALE)
    assert lm_rows[0].functional_impact == pytest.approx(expected)
    assert lm_rows[0].regulatory_impact == pytest.approx(expected)
    assert lm_rows[0].confidence == 1.0  # fixture window is all ACGT


def test_unannotated_variant_is_predicted(session, analysis_run: AnalysisRun, fake_lm):
    make_variant(session, analysis_run, chrom="1", pos=1000, ref="A", alt="G")

    predict_variants(session, analysis_run.id)

    assert len(_predictions(session, "variant_lm")) == 1


def test_alphamissense_score_from_tabix(session, analysis_run: AnalysisRun, refs):
    variant = make_variant(session, analysis_run, chrom="1", pos=1000, ref="A", alt="G")

    assert predict_variants(session, analysis_run.id) == 1  # LM disabled, AlphaMissense only
    rows = _predictions(session, "alphamissense")
    assert len(rows) == 1
    assert rows[0].variant_id == variant.id
    assert rows[0].pathogenicity_probability == pytest.approx(0.2937)


def test_alphamissense_keeps_worst_transcript(session, analysis_run: AnalysisRun, refs):
    make_variant(session, analysis_run, chrom="1", pos=2000, ref="C", alt="T")

    predict_variants(session, analysis_run.id)

    rows = _predictions(session, "alphamissense")
    assert len(rows) == 1
    assert rows[0].pathogenicity_probability == pytest.approx(0.8123)


def test_alphamissense_absent_for_unmatched_allele(session, analysis_run: AnalysisRun, refs):
    make_variant(session, analysis_run, chrom="1", pos=1000, ref="A", alt="T")

    assert predict_variants(session, analysis_run.id) == 0
    assert _predictions(session) == []


def test_lm_inference_failure_degrades_gracefully(
    session, analysis_run: AnalysisRun, refs, monkeypatch
):
    monkeypatch.setattr(predictor, "_load_lm", lambda: ("tokenizer", "model", "cpu"))

    def boom(sequence: str, lm: tuple) -> float:
        raise RuntimeError("MPS backend out of memory")

    monkeypatch.setattr(predictor, "_sequence_log_likelihood", boom)
    make_variant(session, analysis_run, chrom="1", pos=1000, ref="A", alt="G")

    assert predict_variants(session, analysis_run.id) == 1  # AlphaMissense row survives
    assert _predictions(session, "variant_lm") == []


def test_missing_reference_files_degrade_gracefully(
    session, analysis_run: AnalysisRun, monkeypatch, tmp_path
):
    monkeypatch.setattr(settings, "reference_fasta", str(tmp_path / "nope.fa"))
    monkeypatch.setattr(settings, "alphamissense_tsv", str(tmp_path / "nope.tsv.gz"))
    make_variant(session, analysis_run, chrom="1", pos=1000, ref="A", alt="G")

    assert predict_variants(session, analysis_run.id) == 0


def test_reference_mismatch_skips_lm(session, analysis_run: AnalysisRun, fake_lm):
    # Reference base at 1:1000 is A, so a variant claiming ref=T is a coordinate error.
    make_variant(session, analysis_run, chrom="1", pos=1000, ref="T", alt="G")

    predict_variants(session, analysis_run.id)

    assert _predictions(session, "variant_lm") == []


def test_unknown_contig_skips_variant(session, analysis_run: AnalysisRun, fake_lm):
    make_variant(session, analysis_run, chrom="17", pos=1000, ref="A", alt="G")

    assert predict_variants(session, analysis_run.id) == 0


def test_predictions_never_touch_annotation_table(session, analysis_run: AnalysisRun, fake_lm):
    make_variant(session, analysis_run, chrom="1", pos=1000, ref="A", alt="G")
    make_variant(session, analysis_run, chrom="1", pos=2000, ref="C", alt="T")

    assert predict_variants(session, analysis_run.id) > 0
    assert session.scalar(select(func.count()).select_from(Annotation)) == 0
