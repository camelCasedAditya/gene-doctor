import { useMutation, useQuery, useQueryClient } from '@tanstack/react-query'
import { useState } from 'react'
import { Link } from 'react-router-dom'
import {
  getAnalysis,
  getWeights,
  listAsos,
  listDiseases,
  startAnalysis,
  updateWeights,
  type Analysis,
} from '@/api/client'
import { EvidenceBadge, ResearchDisclaimer } from '@/components/EvidenceBadge'
import { PageShell, QueryState } from '@/components/PageShell'
import { Alert, AlertDescription, AlertTitle } from '@/components/ui/alert'
import { Button } from '@/components/ui/button'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import { num, pct } from '@/lib/format'
import { useAnalysisStore } from '@/store/analysis'

export default function Dashboard() {
  return (
    <PageShell
      title="ASO Genome Explorer"
      description="Genome analysis status, variant counts, ranked diseases and candidate ASOs."
      beforeGate={<StartAnalysis />}
    >
      {(analysisId) => <DashboardBody analysisId={analysisId} />}
    </PageShell>
  )
}

function StartAnalysis() {
  const setAnalysisId = useAnalysisStore((s) => s.setAnalysisId)
  const [uploadId, setUploadId] = useState('')
  const mutation = useMutation({
    mutationFn: startAnalysis,
    onSuccess: (data) => setAnalysisId(data.id),
  })

  return (
    <Card className="max-w-xl">
      <CardHeader>
        <CardTitle>Start an analysis</CardTitle>
      </CardHeader>
      <CardContent className="flex flex-col gap-3">
        <form
          className="flex items-end gap-2"
          onSubmit={(e) => {
            e.preventDefault()
            mutation.mutate(uploadId.trim())
          }}
        >
          <label className="text-muted-foreground flex flex-1 flex-col gap-1 text-xs">
            Genome upload ID
            <input
              className="rounded-md border px-2 py-1 font-mono text-sm text-foreground"
              placeholder="from the Upload page"
              value={uploadId}
              onChange={(e) => setUploadId(e.target.value)}
            />
          </label>
          <Button type="submit" size="sm" disabled={!uploadId.trim() || mutation.isPending}>
            Analyze
          </Button>
        </form>
        {mutation.isError ? (
          <Alert variant="destructive">
            <AlertTitle>Could not start analysis</AlertTitle>
            <AlertDescription>{(mutation.error as Error).message}</AlertDescription>
          </Alert>
        ) : null}
      </CardContent>
    </Card>
  )
}

function DashboardBody({ analysisId }: { analysisId: string }) {
  const analysisQuery = useQuery({
    queryKey: ['analysis', analysisId],
    queryFn: () => getAnalysis(analysisId),
    // ponytail: fixed 5s poll — status vocabulary isn't pinned, so we don't guess terminal states.
    refetchInterval: 5000,
  })
  const diseasesQuery = useQuery({
    queryKey: ['diseases', analysisId],
    queryFn: () => listDiseases(analysisId),
  })
  const asosQuery = useQuery({
    queryKey: ['asos', analysisId, '', ''],
    queryFn: () => listAsos({ analysis_id: analysisId }),
  })

  const analysis = analysisQuery.data
  const topDiseases = [...(diseasesQuery.data?.items ?? [])]
    .sort((a, b) => b.rank_score - a.rank_score)
    .slice(0, 5)
  const topAsos = [...(asosQuery.data?.items ?? [])]
    .sort((a, b) => b.confidence - a.confidence)
    .slice(0, 5)

  return (
    <div className="flex flex-col gap-6">
      <ResearchDisclaimer scope="ai" />

      <Card>
        <CardHeader>
          <CardTitle>Analysis {analysisId}</CardTitle>
        </CardHeader>
        <CardContent className="flex flex-col gap-4">
          <QueryState isLoading={analysisQuery.isLoading} error={analysisQuery.error} />
          {analysis ? (
            <>
              <dl className="flex flex-wrap gap-x-8 gap-y-2 text-sm">
                <Field label="Status" value={analysis.status} />
                <Field label="Stage" value={analysis.current_stage ?? '—'} />
                <Field label="Genome upload" value={analysis.genome_upload_id} mono />
              </dl>
              {analysis.error ? (
                <Alert variant="destructive">
                  <AlertTitle>Pipeline error</AlertTitle>
                  <AlertDescription>{analysis.error}</AlertDescription>
                </Alert>
              ) : null}
              <Counts analysis={analysis} />
            </>
          ) : null}
        </CardContent>
      </Card>

      <div className="grid gap-6 lg:grid-cols-2">
        <Card>
          <CardHeader>
            <CardTitle>
              <Link className="underline" to="/diseases">
                Top ranked diseases
              </Link>
            </CardTitle>
          </CardHeader>
          <CardContent className="flex flex-col gap-3">
            <QueryState
              isLoading={diseasesQuery.isLoading}
              error={diseasesQuery.error}
              isEmpty={topDiseases.length === 0}
              emptyText="No diseases ranked yet."
            />
            {topDiseases.map((d) => (
              <div key={d.id} className="flex flex-col gap-1 border-b pb-2 last:border-0">
                <div className="flex items-baseline justify-between gap-2">
                  <span className="font-medium">{d.name}</span>
                  <span className="text-muted-foreground text-xs">
                    rank {num(d.rank_score, 2)} · confidence {pct(d.confidence)}
                  </span>
                </div>
                <div className="flex flex-wrap gap-1">
                  <EvidenceBadge
                    kind="validated"
                    source="clinvar"
                    detail={`${d.known_variant_count} known variants`}
                  />
                  <EvidenceBadge
                    kind="ai"
                    source="variant_lm"
                    detail={`${d.predicted_variant_count} predicted variants`}
                  />
                </div>
              </div>
            ))}
          </CardContent>
        </Card>

        <Card>
          <CardHeader>
            <CardTitle>
              <Link className="underline" to="/asos">
                Top candidate ASOs
              </Link>
            </CardTitle>
          </CardHeader>
          <CardContent className="flex flex-col gap-3">
            <QueryState
              isLoading={asosQuery.isLoading}
              error={asosQuery.error}
              isEmpty={topAsos.length === 0}
              emptyText="No ASO candidates yet."
            />
            {topAsos.map((a) => (
              <div key={a.id} className="flex flex-col gap-1 border-b pb-2 last:border-0">
                <div className="flex items-baseline justify-between gap-2">
                  <span className="font-medium">
                    {a.gene} · {a.mechanism}
                  </span>
                  <span className="text-muted-foreground text-xs">
                    efficacy {pct(a.predicted_efficacy)} · confidence {pct(a.confidence)}
                  </span>
                </div>
                <code className="text-muted-foreground truncate font-mono text-xs">
                  {a.sequence}
                </code>
              </div>
            ))}
          </CardContent>
        </Card>
      </div>

      <Weights />
    </div>
  )
}

function Field({ label, value, mono }: { label: string; value: string; mono?: boolean }) {
  return (
    <div>
      <dt className="text-muted-foreground text-xs">{label}</dt>
      <dd className={mono ? 'font-mono text-sm' : 'text-sm'}>{value}</dd>
    </div>
  )
}

function Counts({ analysis }: { analysis: Analysis }) {
  const c = analysis.counts
  const tiles: { label: string; value: number; tone?: 'validated' | 'ai' }[] = [
    { label: 'Variants', value: c?.variants ?? 0 },
    { label: 'Known (curated)', value: c?.known_variants ?? 0, tone: 'validated' },
    { label: 'Unknown (AI only)', value: c?.unknown_variants ?? 0, tone: 'ai' },
    { label: 'Diseases', value: c?.diseases ?? 0 },
    { label: 'Transcripts', value: c?.transcripts ?? 0 },
    { label: 'ASO candidates', value: c?.asos ?? 0 },
  ]
  return (
    <div className="grid grid-cols-2 gap-3 sm:grid-cols-3 lg:grid-cols-6">
      {tiles.map((t) => (
        <div
          key={t.label}
          className={
            t.tone === 'validated'
              ? 'rounded-lg border border-sky-600/60 bg-sky-600/10 p-3'
              : t.tone === 'ai'
                ? 'rounded-lg border border-dashed border-amber-500/70 bg-amber-500/10 p-3'
                : 'rounded-lg border p-3'
          }
        >
          <div className="text-xl font-medium">{t.value.toLocaleString()}</div>
          <div className="text-muted-foreground text-xs">{t.label}</div>
        </div>
      ))}
    </div>
  )
}

function Weights() {
  const queryClient = useQueryClient()
  const query = useQuery({ queryKey: ['weights'], queryFn: getWeights })
  const mutation = useMutation({
    mutationFn: updateWeights,
    onSuccess: (data) => queryClient.setQueryData(['weights'], data),
  })
  const weights = query.data

  return (
    <Card>
      <CardHeader>
        <CardTitle>Scoring weights</CardTitle>
      </CardHeader>
      <CardContent className="flex flex-col gap-3">
        <QueryState isLoading={query.isLoading} error={query.error} />
        {weights ? (
          <form
            className="flex flex-wrap items-end gap-3"
            onSubmit={(e) => {
              e.preventDefault()
              const data = new FormData(e.currentTarget)
              mutation.mutate({
                clinvar_weight: Number(data.get('clinvar_weight')),
                gwas_weight: Number(data.get('gwas_weight')),
                ai_weight: Number(data.get('ai_weight')),
              })
            }}
          >
            {(
              [
                ['clinvar_weight', 'ClinVar'],
                ['gwas_weight', 'GWAS'],
                ['ai_weight', 'AI prediction'],
              ] as const
            ).map(([name, label]) => (
              <label key={name} className="text-muted-foreground flex flex-col gap-1 text-xs">
                {label}
                <input
                  name={name}
                  type="number"
                  step="0.05"
                  min="0"
                  defaultValue={weights[name]}
                  className="w-24 rounded-md border px-2 py-1 text-sm text-foreground"
                />
              </label>
            ))}
            <Button type="submit" size="sm" disabled={mutation.isPending}>
              Save weights
            </Button>
            {mutation.isError ? (
              <span className="text-destructive text-xs">{(mutation.error as Error).message}</span>
            ) : null}
          </form>
        ) : null}
      </CardContent>
    </Card>
  )
}
