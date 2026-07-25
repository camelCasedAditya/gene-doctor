import { useQuery } from '@tanstack/react-query'
import { useState } from 'react'
import {
  listAsos,
  listTranscripts,
  listVariants,
  type Aso,
  type Transcript,
  type Variant,
} from '@/api/client'
import { EvidenceBadge, ResearchDisclaimer } from '@/components/EvidenceBadge'
import { PageShell, QueryState } from '@/components/PageShell'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import { bp, num, pct } from '@/lib/format'
import { cn } from '@/lib/utils'

const W = 900
const PAD = 40
const EXON_Y = 78
const EXON_H = 26
const ASO_Y = 122

export default function TranscriptViewer() {
  return (
    <PageShell
      title="Transcript Viewer"
      description="Exon structure, splice sites, variant positions and ASO target regions. Click an exon for details."
    >
      {(analysisId) => <Viewer analysisId={analysisId} />}
    </PageShell>
  )
}

function Viewer({ analysisId }: { analysisId: string }) {
  const [selectedId, setSelectedId] = useState<string | null>(null)
  const query = useQuery({
    queryKey: ['transcripts', analysisId],
    queryFn: () => listTranscripts(analysisId),
  })
  const transcripts = [...(query.data?.items ?? [])].sort((a, b) => b.rank_score - a.rank_score)
  const selected = transcripts.find((t) => t.id === selectedId) ?? transcripts[0] ?? null

  return (
    <div className="flex flex-col gap-4">
      <ResearchDisclaimer scope="aso" />
      <QueryState
        isLoading={query.isLoading}
        error={query.error}
        isEmpty={transcripts.length === 0}
        emptyText="No transcripts for this analysis yet."
      />

      {selected ? (
        <div className="grid gap-4 lg:grid-cols-[16rem_1fr]">
          <ul className="flex max-h-[70vh] flex-col gap-1 overflow-y-auto">
            {transcripts.map((t) => (
              <li key={t.id}>
                <button
                  type="button"
                  onClick={() => setSelectedId(t.id)}
                  className={cn(
                    'w-full rounded-lg border px-3 py-2 text-left text-sm focus-visible:ring-3 focus-visible:ring-ring/50 focus-visible:outline-none',
                    selected.id === t.id ? 'bg-muted' : 'hover:bg-muted/50',
                  )}
                >
                  <span className="font-medium">{t.gene}</span>
                  <span className="text-muted-foreground block font-mono text-xs">
                    {t.ensembl_transcript_id}
                  </span>
                </button>
              </li>
            ))}
          </ul>
          <TranscriptDetail analysisId={analysisId} transcript={selected} />
        </div>
      ) : null}
    </div>
  )
}

/** Accepts a bare coordinate or a locus string such as "chr7:117559780-117559800". */
function asoStart(position: Aso['genomic_position']): number | null {
  if (typeof position === 'number') return Number.isFinite(position) ? position : null
  const coords = String(position ?? '').replace(/^[^:]*:/, '')
  const match = /\d[\d,]*/.exec(coords)
  return match ? Number(match[0].replace(/,/g, '')) : null
}

function normChrom(chrom: string): string {
  return chrom.replace(/^chr/i, '')
}

function TranscriptDetail({
  analysisId,
  transcript,
}: {
  analysisId: string
  transcript: Transcript
}) {
  const [activeExon, setActiveExon] = useState<number | null>(null)

  const variantsQuery = useQuery({
    queryKey: ['variants', analysisId, transcript.gene, 'viewer'],
    queryFn: () =>
      listVariants({ analysis_id: analysisId, gene: transcript.gene, limit: 200, offset: 0 }),
  })
  const asosQuery = useQuery({
    queryKey: ['asos', analysisId, '', transcript.id],
    queryFn: () => listAsos({ analysis_id: analysisId, transcript_id: transcript.id }),
  })

  const exons = transcript.exon_structure ?? []
  const candidates = new Set((transcript.candidate_exons ?? []).map((e) => e.exon_number))
  const min = exons.length ? Math.min(...exons.map((e) => e.start)) : 0
  const max = exons.length ? Math.max(...exons.map((e) => e.end)) : 0
  const span = Math.max(1, max - min)
  const x = (pos: number) => PAD + ((pos - min) / span) * (W - 2 * PAD)

  const variants = (variantsQuery.data?.items ?? []).filter(
    (v) =>
      normChrom(v.chrom) === normChrom(transcript.chrom) && v.pos >= min - 100 && v.pos <= max + 100,
  )
  const asos = asosQuery.data?.items ?? []
  const exonOf = (pos: number) => exons.find((e) => pos >= e.start && pos <= e.end) ?? null
  const active = activeExon === null ? null : (exons.find((e) => e.exon_number === activeExon) ?? null)

  return (
    <div className="flex flex-col gap-4">
      <Card>
        <CardHeader>
          <CardTitle className="flex flex-wrap items-center gap-2">
            <span>
              {transcript.gene}{' '}
              <span className="text-muted-foreground font-mono text-xs">
                {transcript.ensembl_transcript_id} · {transcript.chrom} · strand{' '}
                {transcript.strand} · rank {num(transcript.rank_score, 2)}
              </span>
            </span>
            {transcript.evidence_basis === 'validated' ? (
              <EvidenceBadge kind="validated" source="disease-backed target" />
            ) : (
              <EvidenceBadge kind="ai" source="AI-only target" />
            )}
          </CardTitle>
        </CardHeader>
        <CardContent className="flex flex-col gap-3">
          {exons.length === 0 ? (
            <p className="text-muted-foreground text-sm">No exon structure available.</p>
          ) : (
            <div className="overflow-x-auto">
              <svg
                viewBox={`0 0 ${W} 160`}
                className="min-w-[720px] w-full"
                role="img"
                aria-label={`Exon diagram for ${transcript.gene} with variant and ASO target positions`}
              >
                {/* intron line */}
                <line
                  x1={x(min)}
                  x2={x(max)}
                  y1={EXON_Y + EXON_H / 2}
                  y2={EXON_Y + EXON_H / 2}
                  className="stroke-foreground/40"
                  strokeWidth={2}
                />
                {/* exons + splice-site ticks */}
                {exons.map((e) => {
                  const isCandidate = candidates.has(e.exon_number)
                  const left = x(e.start)
                  const width = Math.max(3, x(e.end) - left)
                  return (
                    <g key={e.exon_number}>
                      <rect
                        x={left}
                        y={EXON_Y}
                        width={width}
                        height={EXON_H}
                        rx={2}
                        tabIndex={0}
                        role="button"
                        aria-label={`Exon ${e.exon_number}${isCandidate ? ', ASO candidate exon' : ''}`}
                        onClick={() => setActiveExon(e.exon_number)}
                        onKeyDown={(ev) => {
                          if (ev.key === 'Enter' || ev.key === ' ') {
                            ev.preventDefault()
                            setActiveExon(e.exon_number)
                          }
                        }}
                        className={cn(
                          'cursor-pointer outline-none',
                          isCandidate
                            ? 'fill-sky-600/70 stroke-sky-800'
                            : 'fill-foreground/25 stroke-foreground/50',
                          activeExon === e.exon_number && 'stroke-[3px]',
                        )}
                      />
                      {[e.start, e.end].map((site) => (
                        <line
                          key={site}
                          x1={x(site)}
                          x2={x(site)}
                          y1={EXON_Y - 8}
                          y2={EXON_Y}
                          className="stroke-emerald-600"
                          strokeWidth={1.5}
                        />
                      ))}
                      <text
                        x={left + width / 2}
                        y={EXON_Y + EXON_H + 12}
                        textAnchor="middle"
                        className="fill-muted-foreground text-[9px]"
                      >
                        {e.exon_number}
                      </text>
                    </g>
                  )
                })}
                {/* variants, coloured by evidence class */}
                {variants.map((v) => {
                  const validated = (v.annotations ?? []).length > 0
                  return (
                    <g key={v.id}>
                      <line
                        x1={x(v.pos)}
                        x2={x(v.pos)}
                        y1={44}
                        y2={EXON_Y}
                        className={validated ? 'stroke-sky-600' : 'stroke-amber-500'}
                        strokeWidth={1}
                        strokeDasharray={validated ? undefined : '2 2'}
                      />
                      <circle
                        cx={x(v.pos)}
                        cy={40}
                        r={4}
                        className={
                          validated
                            ? 'fill-sky-600 stroke-sky-800'
                            : 'fill-amber-400 stroke-amber-600'
                        }
                      >
                        <title>
                          {`${v.chrom}:${v.pos} ${v.ref}>${v.alt} — ${
                            validated ? 'validated clinical evidence' : 'AI hypothesis only'
                          }`}
                        </title>
                      </circle>
                    </g>
                  )
                })}
                {/* ASO target regions */}
                {asos.map((a) => {
                  const start = asoStart(a.genomic_position)
                  if (start === null || start < min - span || start > max + span) return null
                  const left = x(start)
                  const width = Math.max(3, x(start + (a.sequence?.length ?? 20)) - left)
                  return (
                    <rect
                      key={a.id}
                      x={left}
                      y={ASO_Y}
                      width={width}
                      height={10}
                      rx={2}
                      className="fill-amber-400/70 stroke-amber-600"
                      strokeDasharray="3 2"
                    >
                      <title>{`ASO candidate ${a.mechanism} (predicted) — ${a.sequence}`}</title>
                    </rect>
                  )
                })}
                <text x={PAD} y={152} className="fill-muted-foreground text-[10px]">
                  {bp(min)}
                </text>
                <text
                  x={W - PAD}
                  y={152}
                  textAnchor="end"
                  className="fill-muted-foreground text-[10px]"
                >
                  {bp(max)}
                </text>
              </svg>
            </div>
          )}

          <div className="text-muted-foreground flex flex-wrap gap-4 text-xs">
            <LegendSwatch className="bg-foreground/25" label="Exon" />
            <LegendSwatch className="bg-sky-600/70" label="Candidate exon" />
            <LegendSwatch className="bg-emerald-600" label="Splice site" />
            <LegendSwatch className="bg-sky-600" label="Variant — validated evidence" round />
            <LegendSwatch className="bg-amber-400" label="Variant — AI hypothesis only" round />
            <LegendSwatch
              className="bg-amber-400/70 border border-dashed border-amber-600"
              label="ASO target (predicted)"
            />
          </div>

          <QueryState isLoading={variantsQuery.isLoading} error={variantsQuery.error} />
          <QueryState isLoading={asosQuery.isLoading} error={asosQuery.error} />

          {(transcript.opportunities ?? []).length > 0 ? (
            <p className="text-sm">
              <span className="text-muted-foreground">Opportunities: </span>
              {transcript.opportunities.join(', ')}
            </p>
          ) : null}
        </CardContent>
      </Card>

      {active ? (
        <Card>
          <CardHeader>
            <CardTitle>
              Exon {active.exon_number}
              {candidates.has(active.exon_number) ? ' — ASO candidate exon' : ''}
            </CardTitle>
          </CardHeader>
          <CardContent className="flex flex-col gap-3 text-sm">
            <p className="text-muted-foreground font-mono text-xs">
              {transcript.chrom}:{bp(active.start)}–{bp(active.end)} ({active.end - active.start + 1}{' '}
              bp)
            </p>
            <ExonVariants variants={variants.filter((v) => exonOf(v.pos)?.exon_number === active.exon_number)} />
            <ExonAsos asos={asos.filter((a) => a.target_exon === active.exon_number)} />
          </CardContent>
        </Card>
      ) : (
        <p className="text-muted-foreground text-sm">Click an exon in the diagram for details.</p>
      )}
    </div>
  )
}

function LegendSwatch({
  className,
  label,
  round,
}: {
  className: string
  label: string
  round?: boolean
}) {
  return (
    <span className="inline-flex items-center gap-1.5">
      <span className={cn('inline-block size-3', round ? 'rounded-full' : 'rounded', className)} />
      {label}
    </span>
  )
}

function ExonVariants({ variants }: { variants: Variant[] }) {
  return (
    <div>
      <h3 className="text-muted-foreground text-xs">Variants in this exon</h3>
      {variants.length === 0 ? (
        <p className="text-muted-foreground text-xs">None.</p>
      ) : (
        <ul className="mt-1 flex flex-col gap-1">
          {variants.map((v) => (
            <li key={v.id} className="font-mono text-xs">
              {v.chrom}:{bp(v.pos)} {v.ref}&gt;{v.alt}{' '}
              <span className={(v.annotations ?? []).length > 0 ? 'text-sky-700' : 'text-amber-700'}>
                {(v.annotations ?? []).length > 0 ? 'validated' : 'AI hypothesis'}
              </span>
            </li>
          ))}
        </ul>
      )}
    </div>
  )
}

function ExonAsos({ asos }: { asos: Aso[] }) {
  return (
    <div>
      <h3 className="text-muted-foreground text-xs">
        ASO candidates targeting this exon (predicted, unvalidated)
      </h3>
      {asos.length === 0 ? (
        <p className="text-muted-foreground text-xs">None.</p>
      ) : (
        <ul className="mt-1 flex flex-col gap-1">
          {asos.map((a) => (
            <li key={a.id} className="text-xs">
              <code className="font-mono">{a.sequence}</code>{' '}
              <span className="text-muted-foreground">
                {a.mechanism} · efficacy {pct(a.predicted_efficacy)} · specificity{' '}
                {pct(a.predicted_specificity)}
              </span>
            </li>
          ))}
        </ul>
      )}
    </div>
  )
}
