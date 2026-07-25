import { NavLink, Route, Routes } from 'react-router-dom'
import AsoDesigner from '@/pages/AsoDesigner'
import Dashboard from '@/pages/Dashboard'
import DiseasePage from '@/pages/DiseasePage'
import TranscriptViewer from '@/pages/TranscriptViewer'
import UploadPage from '@/pages/Upload'
import VariantExplorer from '@/pages/VariantExplorer'

const NAV_ITEMS = [
  { to: '/', label: 'Dashboard', end: true },
  { to: '/upload', label: 'Upload' },
  { to: '/variants', label: 'Variant Explorer' },
  { to: '/diseases', label: 'Diseases' },
  { to: '/transcripts', label: 'Transcripts' },
  { to: '/asos', label: 'ASO Designer' },
]

function App() {
  return (
    <div className="min-h-svh">
      <nav className="flex gap-4 border-b px-6 py-3 text-sm">
        {NAV_ITEMS.map((item) => (
          <NavLink
            key={item.to}
            to={item.to}
            end={item.end}
            className={({ isActive }) =>
              isActive ? 'font-medium text-foreground' : 'text-muted-foreground'
            }
          >
            {item.label}
          </NavLink>
        ))}
      </nav>
      <main className="px-6">
        <Routes>
          <Route path="/" element={<Dashboard />} />
          <Route path="/upload" element={<UploadPage />} />
          <Route path="/variants" element={<VariantExplorer />} />
          <Route path="/diseases" element={<DiseasePage />} />
          <Route path="/transcripts" element={<TranscriptViewer />} />
          <Route path="/asos" element={<AsoDesigner />} />
        </Routes>
      </main>
    </div>
  )
}

export default App
