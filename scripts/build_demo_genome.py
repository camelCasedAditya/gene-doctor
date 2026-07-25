"""Build a small but structurally complete synthetic human genome for an end-to-end demo run.

Real GRCh38 sequence is sliced at loci of genes behind approved or trialled ASO therapies, so
annotation hits real ClinVar/GWAS records and transcript discovery finds real exon structure. Every
chromosome 1-22,X,Y,MT is present so upload validation passes, with filler slices for chromosomes
that carry no target gene.

Variants are planted inside real exons of the target genes.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pysam
from sqlalchemy import create_engine, select
from sqlalchemy.orm import Session

from backend.database.models import RefEnsemblTranscript
from backend.settings import REPO_ROOT, settings

DEFAULT_OUT = REPO_ROOT / "data" / "demo" / "demo_genome.fa"

# chromosome -> (gene, locus start, locus end). Chosen for real ASO relevance.
TARGETS = {
    "5": ("SMN2", 70049638, 70078522),      # nusinersen (spinal muscular atrophy)
    "21": ("SOD1", 31659666, 31668931),     # tofersen (ALS)
    "18": ("TTR", 31557009, 31598833),      # inotersen (hATTR amyloidosis)
    "2": ("APOB", 21001429, 21044073),      # mipomersen (homozygous FH)
    "11": ("HBB", 5225464, 5229395),        # beta-globin disorders
    "17": ("BRCA1", 43044295, 43125364),    # hereditary breast/ovarian cancer
    "4": ("HTT", 3041363, 3120000),         # tominersen (Huntington's)
    "X": ("DMD", 31097677, 31200000),       # eteplirsen (Duchenne)
}
FILLER_LEN = 40_000
SUB = {"A": "G", "G": "A", "C": "T", "T": "C"}
ALL_CHROMS = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]


def exon_positions(session: Session, gene: str, lo: int, hi: int, limit: int = 3) -> list[int]:
    """Positions inside real exons of `gene` that fall within the sliced window."""
    transcript = session.scalars(
        select(RefEnsemblTranscript).where(RefEnsemblTranscript.gene == gene)
    ).first()
    if transcript is None:
        return []
    picks = []
    for exon in transcript.exons:
        if lo + 200 < exon["start"] and exon["end"] < hi - 200:
            picks.append((exon["start"] + exon["end"]) // 2)  # exon body
            if len(picks) >= limit:
                break
    return picks


def write_record(handle, name: str, sequence: str) -> None:
    handle.write(f">{name}\n")
    for i in range(0, len(sequence), 60):
        handle.write(sequence[i : i + 60] + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT, help="output FASTA path")
    out_path = parser.parse_args().out
    out_path.parent.mkdir(parents=True, exist_ok=True)

    fasta = pysam.FastaFile(str(settings.reference_fasta))
    engine = create_engine(settings.database_url.replace("sqlite+aiosqlite", "sqlite"))
    planted: list[tuple[str, int, str, str, str]] = []

    with Session(engine) as session, open(out_path, "w") as out:
        for chrom in ALL_CHROMS:
            if chrom in TARGETS:
                gene, lo, hi = TARGETS[chrom]
                seq = list(fasta.fetch(chrom, lo - 1, hi).upper())
                for pos in exon_positions(session, gene, lo, hi):
                    idx = pos - lo
                    if 0 <= idx < len(seq) and seq[idx] in SUB:
                        ref = seq[idx]
                        seq[idx] = SUB[ref]
                        planted.append((chrom, pos, ref, seq[idx], gene))
                write_record(out, chrom, "".join(seq))
            elif chrom == "MT":
                write_record(out, chrom, fasta.fetch("MT").upper())
            else:
                # Filler: real sequence, mid-chromosome, no planted variants.
                length = fasta.get_reference_length(chrom)
                start = length // 2
                write_record(out, chrom, fasta.fetch(chrom, start, start + FILLER_LEN).upper())

    fasta.close()
    print(f"wrote {out_path}")
    print(f"planted {len(planted)} variants in real exons:")
    for chrom, pos, ref, alt, gene in planted:
        print(f"  {gene:6} chr{chrom}:{pos} {ref}>{alt}")


if __name__ == "__main__":
    main()
