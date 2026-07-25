import { create } from 'zustand'
import { persist } from 'zustand/middleware'

interface AnalysisState {
  analysisId: string | null
  setAnalysisId: (id: string | null) => void
}

/** The one genuinely cross-page piece of client state. Persisted so a reload keeps context. */
export const useAnalysisStore = create<AnalysisState>()(
  persist(
    (set) => ({
      analysisId: null,
      setAnalysisId: (id) => set({ analysisId: id }),
    }),
    { name: 'aso-explorer-analysis' },
  ),
)
