"""Pretty-print pipeline results for run_demo.sh.

Takes the JSON as an argv payload rather than on stdin, so the caller can keep using a plain
command substitution without fighting shell redirection.

Usage: report_run.py {upload|counts|asos} '<json>'
"""

from __future__ import annotations

import json
import sys
from collections import Counter

COUNT_KEYS = (
    "variants",
    "known_variants",
    "unknown_variants",
    "diseases",
    "transcripts",
    "asos",
)


def show_counts(payload: dict) -> None:
    counts = payload["counts"]
    for key in COUNT_KEYS:
        print(f"   {key.replace('_', ' '):18} {counts[key]}")


def show_upload(payload: dict) -> None:
    result = payload["result"]
    warnings = ", ".join(result["warnings"]) or "none"
    print(f"   chromosomes: {len(result['chromosomes_found'])}   warnings: {warnings}")


def show_asos(payload: dict) -> None:
    items = payload["items"]
    if not items:
        print("   no ASO candidates")
        return

    print(f"   ASO genes:  {dict(Counter(a['gene'] for a in items))}")
    print(f"   mechanisms: {dict(Counter(a['mechanism'] for a in items))}")
    print(f"   basis:      {dict(Counter(a['evidence_basis'] for a in items))}")
    print()
    print("   top candidate per gene:")

    seen: set[str] = set()
    for aso in items:
        if aso["gene"] in seen:
            continue
        seen.add(aso["gene"])
        flag = "" if aso["evidence_basis"] == "validated" else "   [AI-only target]"
        print(
            f"     {aso['gene']:7} {aso['sequence']}  {aso['mechanism']:15}"
            f" {aso['genomic_position']:24} GC={aso['gc_pct']:.0f}%"
            f" Tm={aso['tm']:.0f} spec={aso['predicted_specificity']:.2f}{flag}"
        )


def main() -> None:
    mode, raw = sys.argv[1], sys.argv[2]
    payload = json.loads(raw)
    modes = {"counts": show_counts, "asos": show_asos, "upload": show_upload}
    if mode not in modes:
        raise SystemExit(f"unknown mode: {mode}")
    modes[mode](payload)


if __name__ == "__main__":
    main()
