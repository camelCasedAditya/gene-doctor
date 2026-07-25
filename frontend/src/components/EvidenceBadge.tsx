import { FlaskConical, Sparkles, TriangleAlert } from 'lucide-react'
import type { AiPrediction, VariantAnnotation } from '@/api/client'
import { Alert, AlertDescription, AlertTitle } from '@/components/ui/alert'
import { pct } from '@/lib/format'
import { cn } from '@/lib/utils'

const SOURCE_LABELS: Record<string, string> = {
  clinvar: 'ClinVar',
  gwas: 'GWAS',
  variant_lm: 'Variant-LM',
  alphamissense: 'AlphaMissense',
}

function sourceLabel(source: string): string {
  return SOURCE_LABELS[source] ?? source
}

/**
 * The single visual contract for evidence provenance.
 * `validated` = curated clinical evidence (solid, blue). `ai` = model hypothesis (dashed, amber).
 */
export function EvidenceBadge({
  kind,
  source,
  detail,
  className,
}: {
  kind: 'validated' | 'ai'
  source: string
  detail?: string
  className?: string
}) {
  const isAi = kind === 'ai'
  const Icon = isAi ? Sparkles : FlaskConical
  return (
    <span
      title={
        isAi
          ? 'AI prediction — research hypothesis, not a diagnosis or validated clinical finding'
          : 'Validated clinical evidence from a curated database'
      }
      className={cn(
        'inline-flex max-w-full items-center gap-1 rounded px-1.5 py-0.5 text-xs leading-tight',
        isAi
          ? 'border border-dashed border-amber-500/70 bg-amber-500/10 text-amber-700 dark:text-amber-400'
          : 'border border-sky-600/70 bg-sky-600/15 text-sky-800 dark:text-sky-300',
        className,
      )}
    >
      <Icon className="size-3 shrink-0" aria-hidden />
      <span className="font-medium">
        {isAi ? 'AI hypothesis' : 'Validated'} · {sourceLabel(source)}
      </span>
      {detail ? <span className="truncate opacity-80">{detail}</span> : null}
    </span>
  )
}

/** Both evidence classes for one variant, always rendered validated-first. */
export function EvidenceList({
  annotations,
  aiPredictions,
}: {
  annotations?: VariantAnnotation[]
  aiPredictions?: AiPrediction[]
}) {
  const validated = annotations ?? []
  const ai = aiPredictions ?? []
  if (validated.length === 0 && ai.length === 0) {
    return <span className="text-muted-foreground text-xs">No evidence</span>
  }
  return (
    <div className="flex flex-wrap gap-1">
      {validated.map((a, i) => (
        <EvidenceBadge
          key={`v-${i}`}
          kind="validated"
          source={a.source}
          detail={[a.clinical_significance, a.disease].filter(Boolean).join(' — ')}
        />
      ))}
      {ai.map((p, i) => (
        <EvidenceBadge
          key={`a-${i}`}
          kind="ai"
          source={p.source}
          detail={`p=${pct(p.pathogenicity_probability)}, conf ${pct(p.confidence)}`}
        />
      ))}
    </div>
  )
}

export function EvidenceLegend() {
  return (
    <div className="flex flex-wrap items-center gap-2 text-xs">
      <EvidenceBadge kind="validated" source="clinvar" />
      <span className="text-muted-foreground">curated clinical evidence</span>
      <EvidenceBadge kind="ai" source="variant_lm" />
      <span className="text-muted-foreground">model output, unvalidated</span>
    </div>
  )
}

export function ResearchDisclaimer({ scope }: { scope: 'ai' | 'aso' }) {
  return (
    <Alert>
      <TriangleAlert className="text-amber-600" />
      <AlertTitle>Research software — not a clinical diagnostic device</AlertTitle>
      <AlertDescription>
        {scope === 'aso'
          ? 'Every ASO below is a computational candidate hypothesis. Sequences, efficacy and specificity values are model predictions and require experimental validation before any use.'
          : 'Items marked as AI hypotheses are model predictions, not diagnoses or validated clinical findings. Treat them as leads for further investigation only.'}
      </AlertDescription>
    </Alert>
  )
}
