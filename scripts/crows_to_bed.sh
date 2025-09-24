#!/usr/bin/env bash
set -euo pipefail

# tsv_to_bed.sh — Convert TSVs with header (Chr Start End ...) to BED.
# Modes:
#   1) Per-file: tsv_to_bed.sh -T in.tsv[.gz] [-O out_prefix|out.bed] [-s Sample]
#   2) Batch   : tsv_to_bed.sh [-O out_dir/] [-s Sample]   # processes *.tsv / *.tsv.gz in CWD
#
# Behavior:
#   - BED columns: CHROM  (START-1)  END
#   - If -s is provided, add a 4th column (name) with the sample value.
#   - Skips header/blank lines.

SAMPLE=""
TSV_IN=""
OUT_OPT=""

print_help() {
  cat <<EOF
Usage:
  Per-file:
    $0 -T <in.tsv|in.tsv.gz> [-O <out_prefix|out.bed>] [-s <SampleName>]

  Batch (current directory):
    $0 [-O <out_dir/>] [-s <SampleName>]

Options:
  -T, --tsv-file   Input TSV (.tsv or .tsv.gz). Must contain: Chr, Start, End
  -O               Per-file: output prefix or exact .bed; Batch: output directory
  -s               Sample name (optional 4th BED column)
  -h, --help       Show help
EOF
}

# --- Parse args ---
if [[ $# -eq 0 ]]; then
  print_help; exit 1
fi

ARGS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    -T|--tsv-file) TSV_IN="$2"; shift 2 ;;
    -O)            OUT_OPT="${2:-}"; shift 2 ;;
    -s)            SAMPLE="${2:-}"; shift 2 ;;
    -h|--help)     print_help; exit 0 ;;
    --)            shift; break ;;
    -*)            echo "Unknown option: $1" >&2; print_help; exit 1 ;;
    *)             ARGS+=("$1"); shift ;;
  esac
done

# --- Helpers ---
tsv_stream_to_bed() {
  local sample="$1"  # may be empty
  awk -v SAMPLE="$sample" -F'\t' '
    BEGIN { OFS = "\t" }
    NR==1 {
      # Skip header if column 2 is not numeric or header-like
      if ($2 ~ /[^0-9]/ || tolower($1) ~ /chr|chrom/ || tolower($2) ~ /start/ || tolower($3) ~ /end/) { next }
    }
    NF>=3 && $1!="" && $2!="" && $3!="" {
      chr=$1
      start=$2+0
      end=$3+0
      bedStart=start-1
      if (bedStart<0) bedStart=0
      if (SAMPLE != "") {
        print chr, bedStart, end, SAMPLE
      } else {
        print chr, bedStart, end
      }
    }
  '
}

convert_one() {
  local in="$1" out="$2" sample="$3"
  if [[ "$in" == *.gz ]]; then
    if command -v zcat >/dev/null 2>&1; then
      zcat -- "$in" | tsv_stream_to_bed "$sample" > "$out"
    else
      gzip -cd -- "$in" | tsv_stream_to_bed "$sample" > "$out"
    fi
  else
    cat -- "$in" | tsv_stream_to_bed "$sample" > "$out"
  fi
  echo "✅ Wrote $out"
}

# --- Per-file mode ---
if [[ -n "$TSV_IN" ]]; then
  [[ -f "$TSV_IN" ]] || { echo "❌ Input not found: $TSV_IN"; exit 1; }

  if [[ -n "$OUT_OPT" ]]; then
    if [[ "$OUT_OPT" == *.bed ]]; then
      OUT_PATH="$OUT_OPT"
      mkdir -p "$(dirname "$OUT_PATH")"
    else
      OUT_PATH="${OUT_OPT}.bed"
      mkdir -p "$(dirname "$OUT_PATH")"
    fi
  else
    base="$(basename "$TSV_IN")"
    base="${base%.gz}"; base="${base%.tsv}"
    OUT_PATH="${base}.bed"
  fi

  echo "🔄 Converting $TSV_IN → $OUT_PATH"
  convert_one "$TSV_IN" "$OUT_PATH" "$SAMPLE"
  exit 0
fi

# --- Batch mode (current dir) ---
shopt -s nullglob
mapfile -t FILES < <(printf "%s\n" *.tsv *.tsv.gz 2>/dev/null || true)
if [[ ${#FILES[@]} -eq 0 ]]; then
  echo "❌ No .tsv or .tsv.gz files in $(pwd)"; exit 1
fi

OUT_DIR=""
if [[ -n "$OUT_OPT" ]]; then
  OUT_DIR="$OUT_OPT"
  mkdir -p "$OUT_DIR"
fi

for f in "${FILES[@]}"; do
  base="${f%.gz}"; base="${base%.tsv}"
  if [[ -n "$OUT_DIR" ]]; then
    out="${OUT_DIR%/}/${base}.bed"
  else
    out="${base}.bed"
  fi
  echo "🔄 Converting $f → $out"
  convert_one "$f" "$out" "$SAMPLE"
done

echo "✅ Done."
