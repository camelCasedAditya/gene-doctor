# ASO Genome Explorer (AGE)

Local research platform that takes a human genome FASTA and proposes candidate antisense
oligonucleotide (ASO) therapies, combining validated clinical evidence with AI-based prediction.

**This is not a clinical diagnostic device.** Every AI-derived score is a research hypothesis, and
every designed ASO is a candidate requiring experimental validation. Validated evidence (ClinVar,
GWAS Catalog) and AI hypotheses are stored in separate tables and are visually distinguished
everywhere they appear in the UI.

## Pipeline

```
FASTA → validate → align (minimap2) → call variants (bcftools) → annotate (ClinVar + GWAS)
      → predict novel variants (Nucleotide Transformer + AlphaMissense)
      → weighted evidence aggregation → disease ranking
      → transcript discovery → ASO design → dashboard
```

## Deviation from the original spec: Evo 2 → Nucleotide Transformer

Evo 2 hard-requires a CUDA GPU (`ArcInstitute/evo2` issue #67: `Expected a cuda device, but got:
cpu`); its smallest checkpoint needs an NVIDIA GPU with ~48GB VRAM, and there is no CPU or Apple
Silicon path. This machine is an Apple M3 Max, so Evo 2 cannot run here in any form.

Its architectural role is filled by:

- **Nucleotide Transformer** (`InstaDeepAI/nucleotide-transformer-500m-human-ref`) for
  functional/regulatory impact. Runs on Apple Silicon via MPS. The v2 checkpoints are *not* usable:
  they ship their architecture as HF remote code that imports
  `find_pruneable_heads_and_indices`, removed in transformers 5.x.
- **AlphaMissense** precomputed scores for missense pathogenicity (tabix lookup, no inference).

GATK is skipped; the spec marks it optional and bcftools covers variant calling without a JVM.

## Setup

```bash
conda create -y -n gene-doctor -c bioconda -c conda-forge \
  python=3.12 minimap2 samtools bcftools pysam
conda activate gene-doctor
pip install -r backend/requirements.txt
cd frontend && npm install && cd ..
```

### Reference data (one-time, ~11GB total)

Place the GRCh38 primary assembly at `data/reference/`, then prepare indices and datasets:

```bash
# Reference: decompress (bcftools/pysam need a seekable FASTA) and index
cd data/reference
gunzip -k Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
samtools faidx Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
minimap2 -x asm5 -d Homo_sapiens.GRCh38.dna_sm.primary_assembly.mmi \
  Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
cd ../..

# Datasets (each is idempotent - skips what already exists)
python -m backend.annotation.ingest_clinvar         # ClinVar VCF + tabix index
python -m backend.annotation.ingest_alphamissense   # AlphaMissense scores + tabix index
python -m backend.annotation.ingest_gwas            # ~974k GWAS associations → SQLite
python -m backend.annotation.ingest_ensembl_gtf     # ~388k transcripts → SQLite
python -m backend.annotation.ingest_transcriptome   # cDNA + k=12 index for off-target screening
```

The transcriptome index uses `k=12/w=1` deliberately: a 20nt ASO query cannot be seeded by the
genome index's default `k=19`, which returns zero alignments for short queries.

## Run

```bash
python -m uvicorn backend.api.main:app --port 8000   # backend
cd frontend && npm run dev                            # frontend on :5173
```

Then: **Upload** (validate a FASTA by local path) → **Dashboard** (start the analysis, watch stage
progress) → Variant Explorer / Diseases / Transcripts / ASO Designer.

### Demo run

```bash
./run_demo.sh          # build a demo genome, start both servers, run the whole pipeline
./run_demo.sh --path   # just print the upload file path
```

The demo genome is synthetic but structurally complete: all 25 chromosomes, real GRCh38 sequence at
the loci of genes behind real ASO therapies (SMN2/nusinersen, SOD1/tofersen, TTR/inotersen,
APOB/mipomersen, DMD/eteplirsen, HTT/tominersen, HBB, BRCA1), with variants planted in real exons.
It passes validation and the full pipeline finishes in well under a minute, unlike a real
whole-genome run.

## Target provenance: validated vs AI-only

A transcript reaches ASO design one of two ways, recorded in `Transcript.evidence_basis` and shown
on every transcript and ASO in the UI and in every export:

- **`validated`** — reached through a disease backed by ClinVar/GWAS evidence.
- **`ai_hypothesis`** — reached only through AI-scored variants with no validated evidence anywhere.

AI predictions never create a disease, so without the second path a novel variant could never
produce a candidate at all and the AI arm of the pipeline would be unable to contribute a design.
Validated genes always claim the gene budget first; AI-only leads fill what remains.

Transcript discovery consumes diseases in rank order until it has `MAX_TARGET_GENES` (25) distinct
genes, rather than truncating to a fixed number of diseases. A gene in a GWAS-dense locus can carry
over a hundred separate trait associations — HBB really does — and a disease-count cap let one such
gene fill every slot and starve every other implicated gene of a design.

## Scoring

`overall_score` is a weighted mean over the sources that actually have data for a variant:

```
overall_score = Σ(wᵢ · sᵢ) / Σ(wᵢ)     over sources present
```

Renormalizing (rather than treating a missing source as zero) keeps scores comparable across
variants — a variant absent from ClinVar isn't penalized for it. Weights live in the
`evidence_weights` table, defaults 0.45 / 0.30 / 0.25, editable at runtime via
`GET`/`PUT /config/weights`.

## Tests

```bash
pytest backend/tests/ -q
```

Tests never touch the real reference genome, the real databases, or the network: they use synthetic
contigs, in-memory SQLite, and mocked model inference.

## Resetting

```bash
python -m backend.database.reset_results
```

Drops and recreates the per-analysis result tables while preserving the bulk-ingested reference
tables. Needed after changing a result model — SQLAlchemy's `create_all` never `ALTER`s an existing
table, so a stale column will otherwise surface as `no such column` at runtime.

## Known limitations

- **Variant-effect scoring uses a single unmasked forward pass.** The model can see the base it
  scores, which flattens ref-vs-alt deltas (observed range 4e-4 to 3e-2). Masked-marginal scoring of
  the variant token would be more discriminative and costs the same one pass.
- **`functional_impact` and `regulatory_impact` carry the same value.** A window-level LM delta
  cannot separate coding from regulatory disruption; they are not two independent signals.
- **Off-target screening covers near-full-length matches only** (≥18/20 bases), counted as distinct
  genes so isoforms of the intended gene don't register as off-target. It will not catch every
  1–2 mismatch near-match.
- **ASO efficacy is a heuristic** over GC%, Tm, and homopolymer/G-quadruplex penalties — not a
  trained predictor. Candidates passing all bands score identically (1.0), so efficacy does not
  finely rank them.
- **Disease-name matching is string normalization**, not EFO/MONDO ontology mapping, so the same
  disease under two vocabularies won't merge.
- **AI-only variants attach to diseases via shared gene**, since `AiPrediction` has no disease
  field. AI evidence never creates a disease on its own.
# gene-doctor
# gene-doctor
# gene-doctor
