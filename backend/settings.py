import os
import sys
from pathlib import Path

from pydantic_settings import BaseSettings

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data"


def _ensure_tools_on_path() -> None:
    """minimap2/samtools/bcftools are installed alongside the interpreter in the conda env, but
    PATH won't include that directory when uvicorn is launched via the env's python without the
    env being activated. Every subprocess call site (align, variant_call, offtarget) depends on
    these resolving, so fix it once here at import.
    """
    bin_dir = str(Path(sys.executable).parent)
    if bin_dir not in os.environ.get("PATH", "").split(os.pathsep):
        os.environ["PATH"] = bin_dir + os.pathsep + os.environ.get("PATH", "")


_ensure_tools_on_path()


class Settings(BaseSettings):
    database_url: str = f"sqlite+aiosqlite:///{REPO_ROOT / 'gene_doctor.db'}"
    reference_fasta: Path = DATA_DIR / "reference" / "Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
    reference_mmi: Path = DATA_DIR / "reference" / "Homo_sapiens.GRCh38.dna_sm.primary_assembly.mmi"
    clinvar_vcf: Path = DATA_DIR / "clinvar" / "clinvar.vcf.gz"
    alphamissense_tsv: Path = DATA_DIR / "alphamissense" / "AlphaMissense_hg38.tsv.gz"
    gwas_associations_tsv: Path = DATA_DIR / "gwas" / "gwas-catalog-download-associations-alt-full.tsv"
    ensembl_gtf: Path = DATA_DIR / "ensembl_gtf" / "Homo_sapiens.GRCh38.114.gtf.gz"
    # ASOs hybridize to RNA, so off-target screening runs against the transcriptome. Its index uses
    # a small k because 20nt queries cannot be seeded by the genome index's k=19.
    transcriptome_fasta: Path = DATA_DIR / "transcriptome" / "Homo_sapiens.GRCh38.cdna.all.fa.gz"
    transcriptome_mmi: Path = DATA_DIR / "transcriptome" / "cdna_k12.mmi"
    transcript_gene_map: Path = DATA_DIR / "transcriptome" / "transcript_gene_map.tsv"
    data_dir: Path = DATA_DIR


settings = Settings()
