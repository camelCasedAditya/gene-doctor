import { useQuery } from '@tanstack/react-query'
import { useState } from 'react'
import { listAsos, type Aso, type AsoMechanism } from '@/api/client'
import { EvidenceBadge, ResearchDisclaimer } from '@/components/EvidenceBadge'
import { PageShell, QueryState } from '@/components/PageShell'
import { Button } from '@/components/ui/button'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import { num, pct } from '@/lib/format'

const MECHANISMS: AsoMechanism[] = [
  'splice_switch',
  'exon_skip',
  'exon_include',
  'knockdown',
  'allele_specific',
]

const CSV_COLUMNS = [
  'id',
  'gene',
  'evidence_basis',
  'ensembl_transcript_id',
  'sequence',
  'genomic_position',
  'target_exon',
  'mechanism',
  'gc_pct',
  'tm',
  'predicted_efficacy',
  'predicted_specificity',
  'off_target_score',
  'confidence',
] as const

function download(filename: string, mime: string, text: string) {
  const url = URL.createObjectURL(new Blob([text], { type: mime }))
  const link = document.createElement('a')
  link.href = url
  link.download = filename
  link.click()
  URL.revokeObjectURL(url)
}

function toFasta(asos: Aso[]): string {
  return asos
    .map(
      (a) =>
        `>${a.id} gene=${a.gene} basis=${a.evidence_basis} transcript=${a.ensembl_transcript_id} mechanism=${a.mechanism} pos=${a.genomic_position} exon=${a.target_exon ?? 'NA'} predicted_efficacy=${a.predicted_efficacy} CANDIDATE_HYPOTHESIS\n${a.sequence}`,
    )
    .join('\n')
    .concat('\n')
}

function toCsv(asos: Aso[]): string {
  const cell = (value: unknown) => {
    const str = value === null || value === undefined ? '' : String(value)
    return /[",\n]/.test(str) ? `"${str.replace(/"/g, '""')}"` : str
  }
  const rows = asos.map((a) => CSV_COLUMNS.map((c) => cell(a[c])).join(','))
  return [CSV_COLUMNS.join(','), ...rows].join('\n') + '\n'
}

export default function AsoDesigner() {
  return (
    <PageShell
      title="ASO Designer"
      description="Candidate antisense oligonucleotides generated from the ranked transcripts."
    >
      {(analysisId) => <AsoTable analysisId={analysisId} />}
    </PageShell>
  )
}

function AsoTable({ analysisId }: { analysisId: string }) {
  const [mechanism, setMechanism] = useState<AsoMechanism | ''>('')
  const query = useQuery({
    queryKey: ['asos', analysisId, mechanism, ''],
    queryFn: () => listAsos({ analysis_id: analysisId, mechanism: mechanism || undefined }),
  })
  const asos = query.data?.items ?? []
  const stamp = new Date().toISOString().slice(0, 10)

  return (
    <div className="flex flex-col gap-4">
      <ResearchDisclaimer scope="aso" />

      <Card>
        <CardHeader>
          <CardTitle className="flex flex-wrap items-center gap-2">
            Candidates
            <EvidenceBadge
              kind="ai"
              source="ASO design model"
              detail="every row is a computational prediction"
            />
          </CardTitle>
        </CardHeader>
        <CardContent className="flex flex-col gap-3">
          <div className="flex flex-wrap items-center gap-2">
            <span className="text-muted-foreground text-xs">Mechanism:</span>
            <Button
              size="xs"
              variant={mechanism === '' ? 'default' : 'outline'}
              onClick={() => setMechanism('')}
            >
              All
            </Button>
            {MECHANISMS.map((m) => (
              <Button
                key={m}
                size="xs"
                variant={mechanism === m ? 'default' : 'outline'}
                onClick={() => setMechanism(m)}
              >
                {m}
              </Button>
            ))}
          </div>

          <div className="flex flex-wrap items-center gap-2">
            <span className="text-muted-foreground text-xs">Export:</span>
            <Button
              size="xs"
              variant="outline"
              disabled={asos.length === 0}
              onClick={() => download(`aso-candidates-${stamp}.fasta`, 'text/plain', toFasta(asos))}
            >
              FASTA
            </Button>
            <Button
              size="xs"
              variant="outline"
              disabled={asos.length === 0}
              onClick={() => download(`aso-candidates-${stamp}.csv`, 'text/csv', toCsv(asos))}
            >
              CSV
            </Button>
            <Button
              size="xs"
              variant="outline"
              disabled={asos.length === 0}
              onClick={() =>
                download(
                  `aso-candidates-${stamp}.json`,
                  'application/json',
                  JSON.stringify(asos, null, 2),
                )
              }
            >
              JSON
            </Button>
          </div>
        </CardContent>
      </Card>

      <QueryState
        isLoading={query.isLoading}
        error={query.error}
        isEmpty={asos.length === 0}
        emptyText="No ASO candidates for this analysis yet."
      />

      {asos.length > 0 ? (
        <div className="overflow-x-auto rounded-xl ring-1 ring-foreground/10">
          <table className="w-full text-left text-sm">
            <caption className="sr-only">
              Candidate ASO sequences with predicted properties — research hypotheses, not validated
              therapeutics
            </caption>
            <thead className="bg-muted/50 text-muted-foreground text-xs">
              <tr>
                {[
                  'Sequence',
                  'Gene',
                  'Basis',
                  'Transcript',
                  'Exon',
                  'Mechanism',
                  'GC',
                  'Tm',
                  'Efficacy*',
                  'Specificity*',
                  'Off-target*',
                  'Confidence*',
                ].map((h) => (
                  <th key={h} scope="col" className="px-3 py-2 font-medium whitespace-nowrap">
                    {h}
                  </th>
                ))}
              </tr>
            </thead>
            <tbody>
              {asos.map((a) => (
                <tr key={a.id} className="border-t">
                  <td className="px-3 py-2 font-mono text-xs">{a.sequence}</td>
                  <td className="px-3 py-2 whitespace-nowrap">{a.gene}</td>
                  <td className="px-3 py-2">
                    {a.evidence_basis === 'validated' ? (
                      <EvidenceBadge kind="validated" source="disease-backed target" />
                    ) : (
                      <EvidenceBadge kind="ai" source="AI-only target" />
                    )}
                  </td>
                  <td className="text-muted-foreground px-3 py-2 font-mono text-xs whitespace-nowrap">
                    {a.ensembl_transcript_id}
                  </td>
                  <td className="px-3 py-2">{a.target_exon ?? '—'}</td>
                  <td className="px-3 py-2 whitespace-nowrap">{a.mechanism}</td>
                  <td className="px-3 py-2">
                    {a.gc_pct > 1 ? `${num(a.gc_pct, 1)}%` : pct(a.gc_pct)}
                  </td>
                  <td className="px-3 py-2">{num(a.tm, 1)}</td>
                  <td className="px-3 py-2">{pct(a.predicted_efficacy)}</td>
                  <td className="px-3 py-2">{pct(a.predicted_specificity)}</td>
                  <td className="px-3 py-2">{num(a.off_target_score, 2)}</td>
                  <td className="px-3 py-2">{pct(a.confidence)}</td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      ) : null}

      {asos.length > 0 ? (
        <p className="text-muted-foreground text-xs">
          * Model-predicted values — research hypotheses requiring experimental validation.
        </p>
      ) : null}
    </div>
  )
}
