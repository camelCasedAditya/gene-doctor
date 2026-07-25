/** Format a 0-1 score as a percentage; tolerant of missing or already-scaled values. */
export function pct(value: number | null | undefined): string {
  if (typeof value !== 'number' || Number.isNaN(value)) return '—'
  if (value > 1) return value.toFixed(2)
  return `${Math.round(value * 100)}%`
}

export function num(value: number | null | undefined, digits = 1): string {
  if (typeof value !== 'number' || Number.isNaN(value)) return '—'
  return value.toFixed(digits)
}

export function bp(value: number | null | undefined): string {
  return typeof value === 'number' ? value.toLocaleString() : '—'
}
