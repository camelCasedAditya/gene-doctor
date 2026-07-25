"""Download the Ensembl cDNA transcriptome, build a small-k minimap2 index for it, and write a
transcript -> gene_symbol map.

ASOs hybridize to RNA, so the off-target screen belongs against the transcriptome rather than the
genome. The index uses k=12/w=1 because a 20nt ASO cannot be seeded by the default k=19 used for
the genome index - with k=19 a 20nt query returns no alignments at all.

Run via: python -m backend.annotation.ingest_transcriptome
"""

from __future__ import annotations

import gzip
import re
import subprocess
import urllib.request

from backend.settings import settings

URL = (
    "https://ftp.ensembl.org/pub/release-114/fasta/homo_sapiens/cdna/"
    "Homo_sapiens.GRCh38.cdna.all.fa.gz"
)
_GENE_SYMBOL = re.compile(r"gene_symbol:(\S+)")


def ingest_transcriptome() -> int:
    """Download + index the transcriptome and write the transcript->gene map. Returns map size."""
    settings.transcriptome_fasta.parent.mkdir(parents=True, exist_ok=True)
    if not settings.transcriptome_fasta.exists():
        urllib.request.urlretrieve(URL, settings.transcriptome_fasta)

    if not settings.transcriptome_mmi.exists():
        subprocess.run(
            [
                "minimap2",
                "-k", "12",
                "-w", "1",
                "-d", str(settings.transcriptome_mmi),
                str(settings.transcriptome_fasta),
            ],
            check=True,
            capture_output=True,
        )

    count = 0
    with (
        gzip.open(settings.transcriptome_fasta, "rt") as fasta,
        open(settings.transcript_gene_map, "w") as out,
    ):
        for line in fasta:
            if not line.startswith(">"):
                continue
            transcript_id = line[1:].split(None, 1)[0]
            match = _GENE_SYMBOL.search(line)
            gene = match.group(1) if match else ""
            out.write(f"{transcript_id}\t{gene}\n")
            count += 1

    return count


if __name__ == "__main__":
    n = ingest_transcriptome()
    print(f"Transcriptome indexed; wrote {n} transcript->gene entries")
