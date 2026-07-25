#!/usr/bin/env bash
# End-to-end demo run for ASO Genome Explorer.
#
# Builds a small but structurally complete synthetic genome (all 25 chromosomes, real GRCh38
# sequence at the loci of genes behind real ASO therapies, variants planted in real exons), then
# drives the whole pipeline through the HTTP API and prints the results.
#
# Usage:  ./run_demo.sh            full run
#         ./run_demo.sh --path     print the upload file path and exit
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY=/opt/anaconda3/envs/gene-doctor/bin/python3
API=http://127.0.0.1:8000
UI=http://localhost:5173
GENOME="$REPO/data/demo/demo_genome.fa"

bold() { printf '\033[1m%s\033[0m\n' "$1"; }
dim()  { printf '\033[2m%s\033[0m\n' "$1"; }
ok()   { printf '\033[32m✓\033[0m %s\n' "$1"; }
warn() { printf '\033[33m!\033[0m %s\n' "$1"; }
die()  { printf '\033[31m✗\033[0m %s\n' "$1" >&2; exit 1; }

print_path() {
  echo
  bold "UPLOAD FILE PATH  (paste this into the Upload page)"
  printf '\033[36m%s\033[0m\n' "$GENOME"
  echo
}

[ -x "$PY" ] || die "python not found at $PY — activate the gene-doctor conda env or edit PY in this script"

# --path: just show the path (building it first if needed), then stop.
if [ "${1:-}" = "--path" ]; then
  if [ ! -f "$GENOME" ]; then
    warn "demo genome missing, building it first (~1 min)"
    mkdir -p "$(dirname "$GENOME")"
    PYTHONPATH="$REPO" "$PY" "$REPO/scripts/build_demo_genome.py" --out "$GENOME" >/dev/null
  fi
  print_path
  exit 0
fi

bold "═══ ASO Genome Explorer — end-to-end demo ═══"

# ---------------------------------------------------------------- prerequisites
REF="$REPO/data/reference/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
if [ ! -f "$REF" ]; then
  if [ -f "$REF.gz" ]; then
    warn "decompressed reference missing — restoring from .gz (needs ~3GB, takes a minute)"
    gunzip -k "$REF.gz"
  else
    die "reference genome not found at $REF(.gz)"
  fi
fi
ok "reference genome present"

# ---------------------------------------------------------------- demo genome
if [ ! -f "$GENOME" ]; then
  echo
  bold "① Building demo genome"
  mkdir -p "$(dirname "$GENOME")"
  PYTHONPATH="$REPO" "$PY" "$REPO/scripts/build_demo_genome.py" --out "$GENOME"
else
  ok "demo genome already built"
fi
print_path

# ---------------------------------------------------------------- backend
if ! curl -sf --max-time 3 "$API/health" >/dev/null 2>&1; then
  bold "② Starting backend"
  ( cd "$REPO" && nohup "$PY" -m uvicorn backend.api.main:app --port 8000 >/tmp/age-backend.log 2>&1 & )
  for _ in $(seq 1 30); do
    curl -sf --max-time 2 "$API/health" >/dev/null 2>&1 && break
    sleep 1
  done
  curl -sf --max-time 2 "$API/health" >/dev/null 2>&1 || die "backend failed to start — see /tmp/age-backend.log"
fi
ok "backend up on $API"

# ---------------------------------------------------------------- frontend
if ! curl -sf --max-time 3 "$UI" >/dev/null 2>&1; then
  bold "③ Starting frontend"
  ( cd "$REPO/frontend" && nohup npm run dev -- --port 5173 >/tmp/age-frontend.log 2>&1 & )
  for _ in $(seq 1 40); do
    curl -sf --max-time 2 "$UI" >/dev/null 2>&1 && break
    sleep 1
  done
fi
curl -sf --max-time 3 "$UI" >/dev/null 2>&1 && ok "frontend up on $UI" || warn "frontend not reachable (see /tmp/age-frontend.log)"

# ---------------------------------------------------------------- upload
echo
bold "④ Upload + validation"
UPLOAD=$(curl -sf -X POST "$API/upload" -H 'Content-Type: application/json' \
  -d "{\"file_path\": \"$GENOME\"}") || die "upload request failed"
UPLOAD_ID=$("$PY" -c 'import json,sys; print(json.load(sys.stdin)["id"])' <<<"$UPLOAD")
dim "   upload id: $UPLOAD_ID"

while :; do
  S=$(curl -sf "$API/upload/$UPLOAD_ID")
  read -r STATUS PROGRESS <<<"$("$PY" -c 'import json,sys; d=json.load(sys.stdin); print(d["status"], round(d["progress"]*100))' <<<"$S")"
  printf '\r   %s %s%%   ' "$STATUS" "$PROGRESS"
  case "$STATUS" in valid|invalid|error) break ;; esac
  sleep 1
done
echo
[ "$STATUS" = "valid" ] || { "$PY" -m json.tool <<<"$S"; die "validation failed"; }
"$PY" "$REPO/scripts/report_run.py" upload "$S"
ok "genome validated"

# ---------------------------------------------------------------- analyze
echo
bold "⑤ Running pipeline"
ANALYSIS=$(curl -sf -X POST "$API/analyze" -H 'Content-Type: application/json' \
  -d "{\"genome_upload_id\": \"$UPLOAD_ID\"}") || die "analyze request failed"
ANALYSIS_ID=$("$PY" -c 'import json,sys; print(json.load(sys.stdin)["id"])' <<<"$ANALYSIS")
dim "   analysis id: $ANALYSIS_ID"

LAST=""
while :; do
  STATUS=$(curl -sf "$API/analysis/$ANALYSIS_ID" | "$PY" -c 'import json,sys; print(json.load(sys.stdin)["status"])')
  if [ "$STATUS" != "$LAST" ]; then printf '   %s  %s\n' "$(date +%H:%M:%S)" "$STATUS"; LAST="$STATUS"; fi
  case "$STATUS" in done|failed) break ;; esac
  sleep 2
done
[ "$STATUS" = "done" ] || die "pipeline failed — $(curl -sf "$API/analysis/$ANALYSIS_ID" | "$PY" -c 'import json,sys; print(json.load(sys.stdin)["error"])')"

# ---------------------------------------------------------------- results
echo
bold "⑥ Results"
# Heredocs are used to keep the python readable; the JSON arrives as an argv payload rather than on
# stdin, because `python - <<EOF` would consume the heredoc as the program and drop piped input.
COUNTS_JSON=$(curl -sf "$API/analysis/$ANALYSIS_ID")
"$PY" "$REPO/scripts/report_run.py" counts "$COUNTS_JSON"

echo
ASOS_JSON=$(curl -sf "$API/asos?analysis_id=$ANALYSIS_ID")
"$PY" "$REPO/scripts/report_run.py" asos "$ASOS_JSON"

echo
bold "═══ done ═══"
echo "   Open:        $UI"
echo "   Analysis ID: $ANALYSIS_ID"
dim  "   Paste that id into the Analysis ID box to load these results."
print_path
