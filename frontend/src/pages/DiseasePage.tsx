import { useQuery } from '@tanstack/react-query'
import { useState } from 'react'
import { getDisease, listDiseases } from '@/api/client'
import { EvidenceBadge, EvidenceList, ResearchDisclaimer } from '@/components/EvidenceBadge'
import { PageShell, QueryState } from '@/components/PageShell'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import { bp, num, pct } from '@/lib/format'
import { cn } from '@/lib/utils'

export default function DiseasePage() {
  return (
    <PageShell
      title="Diseases"
      description="Ranked candidate diseases, their supporting variants and the split between validated and AI evidence."
    >
      {(analysisId) => <Diseases analysisId={analysisId} />}
    </PageShell>
  )
}

function Diseases({ analysisId }: { analysisId: string }) {
  const [selected, setSelected] = useState<string | null>(null)
  const listQuery = useQuery({
    queryKey: ['diseases', analysisId],
    queryFn: () => listDiseases(analysisId),
  })
  const diseases = [...(listQuery.data?.items ?? [])].sort((a, b) => b.rank_score - a.rank_score)

  return (
    <div className="flex flex-col gap-4">
      <ResearchDisclaimer scope="ai" />
      <QueryState
        isLoading={listQuery.isLoading}
        error={listQuery.error}
        isEmpty={diseases.length === 0}
        emptyText="No diseases ranked for this analysis yet."
      />

      {diseases.length > 0 ? (
        <div className="grid gap-4 lg:grid-cols-[20rem_1fr]">
          <ul className="flex max-h-[70vh] flex-col gap-1 overflow-y-auto">
            {diseases.map((d) => (
              <li key={d.id}>
                <button
                  type="button"
                  onClick={() => setSelected(d.id)}
                  className={cn(
                    'w-full rounded-lg border px-3 py-2 text-left text-sm focus-visible:ring-3 focus-visible:ring-ring/50 focus-visible:outline-none',
                    selected === d.id ? 'bg-muted' : 'hover:bg-muted/50',
                  )}
                >
                  <span className="font-medium">{d.name}</span>
                  <span className="text-muted-foreground block text-xs">
                    rank {num(d.rank_score, 2)} · {d.known_variant_count} validated /{' '}
                    {d.predicted_variant_count} AI
                  </span>
                </button>
              </li>
            ))}
          </ul>

          {selected ? (
            <DiseaseDetail id={selected} />
          ) : (
            <p className="text-muted-foreground text-sm">Select a disease to see its evidence.</p>
          )}
        </div>
      ) : null}
    </div>
  )
}

function DiseaseDetail({ id }: { id: string }) {
  const query = useQuery({ queryKey: ['disease', id], queryFn: () => getDisease(id) })
  const disease = query.data

  return (
    <div className="flex flex-col gap-4">
      <QueryState isLoading={query.isLoading} error={query.error} />
      {disease ? (
        <>
          <Card>
            <CardHeader>
              <CardTitle>{disease.name}</CardTitle>
            </CardHeader>
            <CardContent className="flex flex-col gap-3 text-sm">
              <p className="text-muted-foreground">{disease.description}</p>
              <div className="flex flex-wrap gap-2 text-xs">
                <span className="rounded border px-1.5 py-0.5">
                  rank {num(disease.rank_score, 2)}
                </span>
                <span className="rounded border px-1.5 py-0.5">
                  confidence {pct(disease.confidence)}
                </span>
              </div>
              <div className="flex flex-wrap gap-2">
                <EvidenceBadge
                  kind="validated"
                  source="clinvar"
                  detail={`${disease.known_variant_count} curated variants`}
                />
                <EvidenceBadge
                  kind="ai"
                  source="variant_lm"
                  detail={`${disease.predicted_variant_count} predicted variants`}
                />
              </div>
              <div>
                <h3 className="text-muted-foreground text-xs">Affected genes</h3>
                <div className="mt-1 flex flex-wrap gap-1">
                  {(disease.genes ?? []).length === 0 ? (
                    <span className="text-muted-foreground text-xs">None listed</span>
                  ) : (
                    disease.genes.map((g) => (
                      <span key={g} className="rounded bg-muted px-1.5 py-0.5 font-mono text-xs">
                        {g}
                      </span>
                    ))
                  )}
                </div>
              </div>
            </CardContent>
          </Card>

          <Card>
            <CardHeader>
              <CardTitle>Supporting variants</CardTitle>
            </CardHeader>
            <CardContent>
              {(disease.supporting_variants ?? []).length === 0 ? (
                <p className="text-muted-foreground text-sm">No supporting variants.</p>
              ) : (
                <div className="overflow-x-auto">
                  <table className="w-full text-left text-sm">
                    <caption className="sr-only">
                      Variants supporting {disease.name}, by evidence class
                    </caption>
                    <thead className="text-muted-foreground text-xs">
                      <tr>
                        <th scope="col" className="py-2 pr-3 font-medium">
                          Position
                        </th>
                        <th scope="col" className="py-2 pr-3 font-medium">
                          Gene
                        </th>
                        <th scope="col" className="py-2 pr-3 font-medium">
                          Score
                        </th>
                        <th scope="col" className="py-2 font-medium">
                          Evidence
                        </th>
                      </tr>
                    </thead>
                    <tbody>
                      {disease.supporting_variants.map((v) => (
                        <tr key={v.id} className="border-t align-top">
                          <td className="py-2 pr-3 font-mono text-xs whitespace-nowrap">
                            {v.chrom}:{bp(v.pos)} {v.ref}&gt;{v.alt}
                          </td>
                          <td className="py-2 pr-3 whitespace-nowrap">{v.gene}</td>
                          <td className="py-2 pr-3">{num(v.overall_score, 2)}</td>
                          <td className="py-2">
                            <EvidenceList
                              annotations={v.annotations}
                              aiPredictions={v.ai_predictions}
                            />
                          </td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
              )}
            </CardContent>
          </Card>
        </>
      ) : null}
    </div>
  )
}
