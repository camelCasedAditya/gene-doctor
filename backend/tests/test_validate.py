import gzip

from backend.genomics.validate import EXPECTED_CHROMOSOMES, validate_fasta


def _write_fasta(tmp_path, records: dict[str, str], gz: bool = False, extra: str = ""):
    body = "".join(f">{name}\n{seq}\n" for name, seq in records.items()) + extra
    path = tmp_path / ("genome.fa.gz" if gz else "genome.fa")
    if gz:
        with gzip.open(path, "wt") as f:
            f.write(body)
    else:
        path.write_text(body)
    return str(path)


def _full_genome_records() -> dict[str, str]:
    return {name: "ACGT" * 10 for name in sorted(EXPECTED_CHROMOSOMES)}


def test_valid_full_genome(tmp_path):
    path = _write_fasta(tmp_path, _full_genome_records())
    result = validate_fasta(path)
    assert result.valid
    assert not result.errors
    assert set(result.chromosomes_found) == EXPECTED_CHROMOSOMES


def test_valid_gzipped(tmp_path):
    path = _write_fasta(tmp_path, _full_genome_records(), gz=True)
    result = validate_fasta(path)
    assert result.valid


def test_missing_y_is_warning_only(tmp_path):
    records = _full_genome_records()
    del records["Y"]
    path = _write_fasta(tmp_path, records)
    result = validate_fasta(path)
    assert result.valid
    assert result.missing == ["Y"]
    assert any("Y" in w for w in result.warnings)
    assert not result.errors


def test_missing_autosome_fails(tmp_path):
    records = _full_genome_records()
    del records["7"]
    path = _write_fasta(tmp_path, records)
    result = validate_fasta(path)
    assert not result.valid
    assert "7" in result.missing
    assert any("Missing chromosomes" in e for e in result.errors)


def test_duplicate_chromosome_fails(tmp_path):
    body = ">1\nACGT\n>1\nACGT\n"
    path = tmp_path / "genome.fa"
    path.write_text(body)
    result = validate_fasta(str(path))
    assert not result.valid
    assert "1" in result.duplicates


def test_invalid_characters_fail(tmp_path):
    records = _full_genome_records()
    records["1"] = "ACGTXYZ"
    path = _write_fasta(tmp_path, records)
    result = validate_fasta(path)
    assert not result.valid
    assert any("Invalid sequence characters" in e for e in result.errors)


def test_chr_prefix_normalized(tmp_path):
    body = "".join(f">chr{name}\nACGT\n" for name in sorted(EXPECTED_CHROMOSOMES))
    path = tmp_path / "genome.fa"
    path.write_text(body)
    result = validate_fasta(str(path))
    assert result.valid
    assert set(result.chromosomes_found) == EXPECTED_CHROMOSOMES


def test_extra_scaffold_allowed(tmp_path):
    records = _full_genome_records()
    records["KI270728.1"] = "ACGT"
    path = _write_fasta(tmp_path, records)
    result = validate_fasta(path)
    assert result.valid
    assert "KI270728.1" in result.extra_contigs


def test_progress_callback_monotonic_and_ends_at_one(tmp_path):
    path = _write_fasta(tmp_path, _full_genome_records())
    seen = []
    validate_fasta(path, on_progress=seen.append)
    assert seen == sorted(seen)
    assert seen[-1] == 1.0
