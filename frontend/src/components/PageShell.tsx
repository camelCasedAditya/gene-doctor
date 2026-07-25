import { useState, type ReactNode } from 'react'
import { Link } from 'react-router-dom'
import { Alert, AlertDescription, AlertTitle } from '@/components/ui/alert'
import { Button } from '@/components/ui/button'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import { useAnalysisStore } from '@/store/analysis'

/** Page frame + the shared analysis selector. Children only render once an analysis is chosen. */
export function PageShell({
  title,
  description,
  beforeGate,
  children,
}: {
  title: string
  description?: string
  /** Rendered whether or not an analysis is selected (e.g. the "start analysis" form). */
  beforeGate?: ReactNode
  children: (analysisId: string) => ReactNode
}) {
  const analysisId = useAnalysisStore((s) => s.analysisId)
  const setAnalysisId = useAnalysisStore((s) => s.setAnalysisId)
  const [draft, setDraft] = useState('')

  return (
    <div className="mx-auto w-full max-w-6xl py-8 text-left">
      <header className="mb-6 flex flex-wrap items-end justify-between gap-4">
        <div>
          <h1 className="text-2xl font-medium">{title}</h1>
          {description ? (
            <p className="text-muted-foreground mt-1 text-sm">{description}</p>
          ) : null}
        </div>
        <form
          className="flex items-end gap-2"
          onSubmit={(e) => {
            e.preventDefault()
            setAnalysisId(draft.trim() || null)
          }}
        >
          <label className="text-muted-foreground flex flex-col gap-1 text-xs">
            Analysis ID
            <input
              className="w-56 rounded-md border px-2 py-1 font-mono text-sm text-foreground"
              placeholder={analysisId ?? 'paste analysis id'}
              value={draft}
              onChange={(e) => setDraft(e.target.value)}
            />
          </label>
          <Button type="submit" size="sm" disabled={!draft.trim()}>
            Load
          </Button>
          {analysisId ? (
            <Button
              type="button"
              size="sm"
              variant="ghost"
              onClick={() => {
                setDraft('')
                setAnalysisId(null)
              }}
            >
              Clear
            </Button>
          ) : null}
        </form>
      </header>

      {beforeGate ? <div className="mb-6">{beforeGate}</div> : null}
      {analysisId ? children(analysisId) : <NoAnalysis />}
    </div>
  )
}

function NoAnalysis() {
  return (
    <Card className="max-w-xl">
      <CardHeader>
        <CardTitle>No analysis selected</CardTitle>
      </CardHeader>
      <CardContent className="text-muted-foreground flex flex-col gap-3 text-sm">
        <p>
          Validate a genome FASTA on the <Link className="underline" to="/upload">Upload</Link> page
          and start an analysis from the <Link className="underline" to="/">Dashboard</Link>, or
          paste an existing analysis ID above.
        </p>
      </CardContent>
    </Card>
  )
}

/** Shared loading / error / empty rendering for a TanStack query. Renders nothing when data is ready. */
export function QueryState({
  isLoading,
  error,
  isEmpty,
  emptyText = 'No results.',
}: {
  isLoading: boolean
  error?: unknown
  isEmpty?: boolean
  emptyText?: string
}) {
  if (isLoading) return <p className="text-muted-foreground text-sm">Loading…</p>
  if (error) {
    return (
      <Alert variant="destructive">
        <AlertTitle>Request failed</AlertTitle>
        <AlertDescription>
          {error instanceof Error ? error.message : String(error)}
        </AlertDescription>
      </Alert>
    )
  }
  if (isEmpty) return <p className="text-muted-foreground text-sm">{emptyText}</p>
  return null
}
