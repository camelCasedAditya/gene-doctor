import { useMutation, useQuery } from '@tanstack/react-query'
import { useState } from 'react'
import { createUpload, getUpload } from '@/api/client'
import { Alert, AlertDescription, AlertTitle } from '@/components/ui/alert'
import { Button } from '@/components/ui/button'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import { Progress } from '@/components/ui/progress'

export default function UploadPage() {
  const [filePath, setFilePath] = useState('')
  const [uploadId, setUploadId] = useState<string | null>(null)

  const createMutation = useMutation({
    mutationFn: createUpload,
    onSuccess: (data) => setUploadId(data.id),
  })

  const statusQuery = useQuery({
    queryKey: ['upload', uploadId],
    queryFn: () => getUpload(uploadId!),
    enabled: uploadId !== null,
    refetchInterval: (query) => {
      const status = query.state.data?.status
      return status === 'valid' || status === 'invalid' || status === 'error' ? false : 1000
    },
  })

  const status = statusQuery.data

  return (
    <Card className="mx-auto mt-12 w-full max-w-xl text-left">
      <CardHeader>
        <CardTitle>Upload Genome</CardTitle>
      </CardHeader>
      <CardContent className="flex flex-col gap-4">
        <div className="flex gap-2">
          <input
            className="flex-1 rounded-md border px-3 py-2 text-sm"
            placeholder="/path/to/genome.fa.gz"
            value={filePath}
            onChange={(e) => setFilePath(e.target.value)}
          />
          <Button
            disabled={!filePath || createMutation.isPending}
            onClick={() => createMutation.mutate(filePath)}
          >
            Validate
          </Button>
        </div>

        {createMutation.isError && (
          <Alert variant="destructive">
            <AlertTitle>Request failed</AlertTitle>
            <AlertDescription>{(createMutation.error as Error).message}</AlertDescription>
          </Alert>
        )}

        {status && status.status !== 'valid' && status.status !== 'invalid' && (
          <div className="flex flex-col gap-2">
            <Progress value={status.progress * 100} />
            <p className="text-muted-foreground text-sm">Status: {status.status}</p>
          </div>
        )}

        {status?.status === 'valid' && (
          <Alert>
            <AlertTitle>Valid human genome FASTA</AlertTitle>
            <AlertDescription>
              Found chromosomes: {status.result?.chromosomes_found.join(', ')}
              {status.result?.warnings.length ? (
                <div className="mt-1">Warnings: {status.result.warnings.join('; ')}</div>
              ) : null}
            </AlertDescription>
          </Alert>
        )}

        {status?.status === 'invalid' && (
          <Alert variant="destructive">
            <AlertTitle>Validation failed</AlertTitle>
            <AlertDescription>{status.result?.errors.join('; ')}</AlertDescription>
          </Alert>
        )}
      </CardContent>
    </Card>
  )
}
