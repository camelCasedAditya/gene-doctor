"""Download AlphaMissense's precomputed hg38 scores and tabix-index them. The file is shipped
already BGZF-compressed, so it's directly tabix-indexable with no re-compression step.

Run via: python -m backend.annotation.ingest_alphamissense
"""

from __future__ import annotations

import subprocess
import urllib.request

from backend.settings import settings

URL = "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz"


def ingest_alphamissense() -> None:
    settings.alphamissense_tsv.parent.mkdir(parents=True, exist_ok=True)
    if not settings.alphamissense_tsv.exists():
        urllib.request.urlretrieve(URL, settings.alphamissense_tsv)

    tbi_path = settings.alphamissense_tsv.with_suffix(settings.alphamissense_tsv.suffix + ".tbi")
    if not tbi_path.exists():
        subprocess.run(
            ["tabix", "-s", "1", "-b", "2", "-e", "2", "-c", "#", str(settings.alphamissense_tsv)],
            check=True,
        )


if __name__ == "__main__":
    ingest_alphamissense()
    print(f"AlphaMissense scores ready at {settings.alphamissense_tsv}")
