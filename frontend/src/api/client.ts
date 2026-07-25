const API_BASE = import.meta.env.VITE_API_BASE ?? 'http://localhost:8000'

async function request<T>(path: string, init?: RequestInit): Promise<T> {
  const res = await fetch(`${API_BASE}${path}`, {
    headers: { 'Content-Type': 'application/json' },
    ...init,
  })
  if (!res.ok) {
    const body = await res.text()
    throw new Error(`${res.status} ${res.statusText}: ${body}`)
  }
  return res.json() as Promise<T>
}

export interface UploadStatus {
  id: string
  status: 'pending' | 'validating' | 'valid' | 'invalid' | 'error'
  progress: number
  result: {
    valid: boolean
    chromosomes_found: string[]
    duplicates: string[]
    missing: string[]
    extra_contigs: string[]
    errors: string[]
    warnings: string[]
  } | null
}

export function createUpload(filePath: string): Promise<UploadStatus> {
  return request<UploadStatus>('/upload', {
    method: 'POST',
    body: JSON.stringify({ file_path: filePath }),
  })
}

export function getUpload(id: string): Promise<UploadStatus> {
  return request<UploadStatus>(`/upload/${id}`)
}

function qs(params: Record<string, string | number | undefined>): string {
  const search = new URLSearchParams()
  for (const [key, value] of Object.entries(params)) {
    if (value !== undefined && value !== '') search.set(key, String(value))
  }
  const str = search.toString()
  return str ? `?${str}` : ''
}

export interface AnalysisCounts {
  variants: number
  known_variants: number
  unknown_variants: number
  diseases: number
  transcripts: number
  asos: number
}

export interface Analysis {
  id: string
  status: string
  current_stage: string | null
  error: string | null
  genome_upload_id: string
  counts: AnalysisCounts
}

/** Validated clinical evidence (ClinVar / GWAS). */
export interface VariantAnnotation {
  source: 'clinvar' | 'gwas'
  disease: string
  clinical_significance: string
  publications: string[]
  evidence_strength: number
}

/** Model-generated hypothesis — NOT clinical evidence. */
export interface AiPrediction {
  source: 'variant_lm' | 'alphamissense'
  functional_impact: string
  regulatory_impact: string
  pathogenicity_probability: number
  confidence: number
}

export interface Variant {
  id: string
  chrom: string
  pos: number
  ref: string
  alt: string
  gene: string
  transcript_id: string | null
  overall_score: number
  annotations: VariantAnnotation[]
  ai_predictions: AiPrediction[]
}

export interface Disease {
  id: string
  name: string
  description: string
  rank_score: number
  confidence: number
  known_variant_count: number
  predicted_variant_count: number
  genes: string[]
}

export interface Exon {
  start: number
  end: number
  exon_number: number
}

export interface Transcript {
  id: string
  gene: string
  ensembl_transcript_id: string
  chrom: string
  strand: string
  exon_structure: Exon[]
  candidate_exons: Exon[]
  opportunities: string[]
  rank_score: number
  evidence_basis: EvidenceBasis
}

/** Whether a target was reached via validated clinical evidence or via AI evidence alone. */
export type EvidenceBasis = 'validated' | 'ai_hypothesis'

export type AsoMechanism =
  | 'splice_switch'
  | 'exon_skip'
  | 'exon_include'
  | 'knockdown'
  | 'allele_specific'

export interface Aso {
  id: string
  transcript_id: string
  ensembl_transcript_id: string
  gene: string
  sequence: string
  genomic_position: number | string
  target_exon: number | null
  mechanism: AsoMechanism
  gc_pct: number
  tm: number
  predicted_efficacy: number
  predicted_specificity: number
  off_target_score: number
  confidence: number
  evidence_basis: EvidenceBasis
}

export interface ScoringWeights {
  clinvar_weight: number
  gwas_weight: number
  ai_weight: number
}

export function getAnalysis(id: string): Promise<Analysis> {
  return request<Analysis>(`/analysis/${id}`)
}

export function startAnalysis(genomeUploadId: string): Promise<{ id: string; status: string }> {
  return request<{ id: string; status: string }>('/analyze', {
    method: 'POST',
    body: JSON.stringify({ genome_upload_id: genomeUploadId }),
  })
}

export interface VariantFilters {
  analysis_id: string
  gene?: string
  disease?: string
  chrom?: string
  evidence?: 'known' | 'predicted'
  min_score?: number
  limit?: number
  offset?: number
}

export function listVariants(filters: VariantFilters): Promise<{ total: number; items: Variant[] }> {
  return request<{ total: number; items: Variant[] }>(`/variants${qs({ ...filters })}`)
}

export function listDiseases(analysisId: string): Promise<{ items: Disease[] }> {
  return request<{ items: Disease[] }>(`/diseases${qs({ analysis_id: analysisId })}`)
}

export function getDisease(id: string): Promise<Disease & { supporting_variants: Variant[] }> {
  return request<Disease & { supporting_variants: Variant[] }>(`/diseases/${id}`)
}

export function listTranscripts(analysisId: string): Promise<{ items: Transcript[] }> {
  return request<{ items: Transcript[] }>(`/transcripts${qs({ analysis_id: analysisId })}`)
}

export function listAsos(params: {
  analysis_id: string
  mechanism?: string
  transcript_id?: string
}): Promise<{ items: Aso[] }> {
  return request<{ items: Aso[] }>(`/asos${qs({ ...params })}`)
}

export function getWeights(): Promise<ScoringWeights> {
  return request<ScoringWeights>('/config/weights')
}

export function updateWeights(weights: ScoringWeights): Promise<ScoringWeights> {
  return request<ScoringWeights>('/config/weights', {
    method: 'PUT',
    body: JSON.stringify(weights),
  })
}
