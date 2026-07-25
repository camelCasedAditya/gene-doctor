import { useQuery } from '@tanstack/react-query'
import { useState } from 'react'
import { listVariants } from '@/api/client'
import { EvidenceLegend, EvidenceList, ResearchDisclaimer } from '@/components/EvidenceBadge'
import { PageShell, QueryState } from '@/components/PageShell'
import { Button } from '@/components/ui/button'
import { Card, CardContent } from '@/components/ui/card'
import { bp, num } from '@/lib/format'

const LIMIT = 50

interface Filters {
  gene: string
  disease: string
  chrom: string
  evidence: '' | 'known' | 'predicted'
  minScore: string
}

const EMPTY: Filters = { gene: '', disease: '', chrom: '', evidence: '', minScore: '' }

export default function VariantExplorer() {
  return (
    <PageShell
      title="Variant Explorer"
      description="Search and filter variants by gene, disease, chromosome and score."
    >
      {(analysisId) => <VariantTable analysisId={analysisId} />}
    </PageShell>
  )
}

function VariantTable({ analysisId }: { analysisId: string }) {
  const [draft, setDraft] = useState<Filters>(EMPTY)
  const [filters, setFilters] = useState<Filters>(EMPTY)
  const [offset, setOffset] = useState(0)

  const query = useQuery({
    queryKey: ['variants', analysisId, filters, offset],
    queryFn: () =>
      listVariants({
        analysis_id: analysisId,
        gene: filters.gene || undefined,
        disease: filters.disease || undefined,
        chrom: filters.chrom || undefined,
        evidence: filters.evidence || undefined,
        min_score: filters.minScore ? Number(filters.minScore) : undefined,
        limit: LIMIT,
        offset,
      }),
  })

  const items = query.data?.items ?? []
  const total = query.data?.total ?? 0

  function apply(next: Filters) {
    setDraft(next)
    setFilters(next)
    setOffset(0)
  }

  return (
    <div className="flex flex-col gap-4">
      <ResearchDisclaimer scope="ai" />

      <Card>
        <CardContent className="flex flex-col gap-3">
          <form
            className="flex flex-wrap items-end gap-3"
            onSubmit={(e) => {
              e.preventDefault()
              apply(draft)
            }}
          >
            <TextFilter
              label="Gene"
              value={draft.gene}
              onChange={(gene) => setDraft({ ...draft, gene })}
            />
            <TextFilter
              label="Disease"
              value={draft.disease}
              onChange={(disease) => setDraft({ ...draft, disease })}
            />
            <TextFilter
              label="Chromosome"
              placeholder="chr7"
              value={draft.chrom}
              onChange={(chrom) => setDraft({ ...draft, chrom })}
            />
            <label className="text-muted-foreground flex flex-col gap-1 text-xs">
              Min score
              <input
                type="number"
                step="0.05"
                min="0"
                max="1"
                className="w-24 rounded-md border px-2 py-1 text-sm text-foreground"
                value={draft.minScore}
                onChange={(e) => setDraft({ ...draft, minScore: e.target.value })}
              />
            </label>
            <Button type="submit" size="sm">
              Apply
            </Button>
            <Button type="button" size="sm" variant="ghost" onClick={() => apply(EMPTY)}>
              Reset
            </Button>
          </form>

          <div className="flex flex-wrap items-center gap-2">
            <span className="text-muted-foreground text-xs">Evidence:</span>
            {(
              [
                ['', 'All'],
                ['known', 'Validated only'],
                ['predicted', 'AI predicted only'],
              ] as const
            ).map(([value, label]) => (
              <Button
                key={label}
                type="button"
                size="xs"
                variant={filters.evidence === value ? 'default' : 'outline'}
                onClick={() => apply({ ...draft, evidence: value })}
              >
                {label}
              </Button>
            ))}
            <Button
              type="button"
              size="xs"
              variant={filters.minScore === '0.7' ? 'default' : 'outline'}
              onClick={() =>
                apply({ ...draft, minScore: filters.minScore === '0.7' ? '' : '0.7' })
              }
            >
              High confidence (≥ 0.7)
            </Button>
          </div>

          <EvidenceLegend />
        </CardContent>
      </Card>

      <QueryState
        isLoading={query.isLoading}
        error={query.error}
        isEmpty={items.length === 0}
        emptyText="No variants match these filters."
      />

      {items.length > 0 ? (
        <>
          <div className="overflow-x-auto rounded-xl ring-1 ring-foreground/10">
            <table className="w-full text-left text-sm">
              <caption className="sr-only">Variants with clinical and AI evidence</caption>
              <thead className="bg-muted/50 text-muted-foreground text-xs">
                <tr>
                  <th scope="col" className="px-3 py-2 font-medium">
                    Position
                  </th>
                  <th scope="col" className="px-3 py-2 font-medium">
                    Change
                  </th>
                  <th scope="col" className="px-3 py-2 font-medium">
                    Gene
                  </th>
                  <th scope="col" className="px-3 py-2 font-medium">
                    Transcript
                  </th>
                  <th scope="col" className="px-3 py-2 font-medium">
                    Score
                  </th>
                  <th scope="col" className="px-3 py-2 font-medium">
                    Evidence
                  </th>
                </tr>
              </thead>
              <tbody>
                {items.map((v) => (
                  <tr key={v.id} className="border-t align-top">
                    <td className="px-3 py-2 font-mono text-xs whitespace-nowrap">
                      {v.chrom}:{bp(v.pos)}
                    </td>
                    <td className="px-3 py-2 font-mono text-xs">
                      {v.ref}&gt;{v.alt}
                    </td>
                    <td className="px-3 py-2 whitespace-nowrap">{v.gene}</td>
                    <td className="text-muted-foreground px-3 py-2 font-mono text-xs">
                      {v.transcript_id ?? '—'}
                    </td>
                    <td className="px-3 py-2">{num(v.overall_score, 2)}</td>
                    <td className="px-3 py-2">
                      <EvidenceList annotations={v.annotations} aiPredictions={v.ai_predictions} />
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>

          <div className="flex items-center gap-3 text-sm">
            <Button
              size="sm"
              variant="outline"
              disabled={offset === 0}
              onClick={() => setOffset(Math.max(0, offset - LIMIT))}
            >
              Previous
            </Button>
            <span className="text-muted-foreground">
              {offset + 1}–{offset + items.length} of {total.toLocaleString()}
            </span>
            <Button
              size="sm"
              variant="outline"
              disabled={offset + items.length >= total}
              onClick={() => setOffset(offset + LIMIT)}
            >
              Next
            </Button>
          </div>
        </>
      ) : null}
    </div>
  )
}

function TextFilter({
  label,
  value,
  placeholder,
  onChange,
}: {
  label: string
  value: string
  placeholder?: string
  onChange: (value: string) => void
}) {
  return (
    <label className="text-muted-foreground flex flex-col gap-1 text-xs">
      {label}
      <input
        className="w-36 rounded-md border px-2 py-1 text-sm text-foreground"
        placeholder={placeholder}
        value={value}
        onChange={(e) => onChange(e.target.value)}
      />
    </label>
  )
}
