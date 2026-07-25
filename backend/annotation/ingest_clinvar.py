"""Download the ClinVar GRCh38 VCF + tabix index. NCBI ships it pre-sorted and pre-indexed, so
no parsing/ingestion step is needed - annotate.py queries it directly via pysam.VariantFile.

Run via: python -m backend.annotation.ingest_clinvar
"""

from __future__ import annotations

import urllib.request

from backend.settings import settings

VCF_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
TBI_URL = VCF_URL + ".tbi"


def ingest_clinvar() -> None:
    settings.clinvar_vcf.parent.mkdir(parents=True, exist_ok=True)
    if not settings.clinvar_vcf.exists():
        urllib.request.urlretrieve(VCF_URL, settings.clinvar_vcf)
    tbi_path = settings.clinvar_vcf.with_suffix(settings.clinvar_vcf.suffix + ".tbi")
    if not tbi_path.exists():
        urllib.request.urlretrieve(TBI_URL, tbi_path)


if __name__ == "__main__":
    ingest_clinvar()
    print(f"ClinVar VCF ready at {settings.clinvar_vcf}")
