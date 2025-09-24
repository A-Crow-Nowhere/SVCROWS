#!/usr/bin/env bash
set -euo pipefail

# vcf_to_bed.sh — Convert VCF (optionally .gz) to BED.
# Modes:
#   1) Per-file: vcf_to_bed.sh -V in.vcf[.gz] [-O out_prefix|out.bed] [-s Sample]
#   2) Batch   : vcf_to_bed.sh [-O out_dir/] [-s Sample]   # processes *.vcf / *.vcf.gz in CWD
#
# Output:
#   - 3-col BED: CHROM  (POS-1)  END
#   - If -s is provided, add 4th column (name) with the Sample.
#
# END determination:
#   - If INFO/END exists, use it.
#   - Else if INS (SVTYPE=INS or symbolic <INS> or REF len=1 & ALT longer): END=POS
#   - Else END = POS + len(REF) - 1   (SNP ⇒ END=POS; MNV/DEL spans REF length)

SAMPLE=""
VCF_IN=""
OUT_OPT=""

print_help() {
  cat <<EOF
Usage:
  Per-file:
    $0 -V <in.vcf|in.vcf.gz> [-O <out_prefix|out.bed>] [-s <SampleName>]

  Batch (current directory):
    $0 [-O <out_dir/>] [-s <SampleName>]

Options:
  -V, --vcf-file  Input VCF (.vcf or .vcf.gz)
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
    -V|--vcf-file) VCF_IN="$2"; shift 2 ;;
    -O)            OUT_OPT="${2:-}"; shift 2 ;;
    -s)            SAMPLE="${2:-}"; shift 2 ;;
    -h|--help)     print_help; exit 0 ;;
    --)            shift; break ;;
    -*)            echo "Unknown option: $1" >&2; print_help; exit 1 ;;
    *)             ARGS+=("$1"); shift ;;
  esac
done

# --- Converter (stdin VCF -> stdout BED) ---
vcf_stream_to_bed() {
  local sample="$1"  # may be empty
  awk -v SAMPLE="$sample" -F'\t' '
    BEGIN {
      OFS = "\t"
    }
    # Skip headers
    /^#/ { next }

    # Parse VCF columns
    {
      chrom=$1; pos=$2+0; id=$3; ref=$4; alt=$5; info=$8
      if (chrom=="" || pos==0 || ref=="" || alt=="") next

      bedStart = pos - 1
      if (bedStart < 0) bedStart = 0

      # Extract END and SVTYPE from INFO (simple parser)
      endTag = ""; svtype = ""
      n = split(info, arr, ";")
      for (i=1; i<=n; i++) {
        split(arr[i], kv, "=")
        key = kv[1]
        val = kv[2]
        if (key=="END")    endTag = val
        if (key=="SVTYPE") svtype = val
      }

      # Use first ALT allele only for span logic
      split(alt, alts, ","); alt1 = alts[1]

      # Symbolic ALT?
      isSymbolic = (alt1 ~ /^<.*>$/) ? 1 : 0

      # Compute END
      if (endTag ~ /^[0-9]+$/) {
        end = endTag + 0
      } else {
        rlen = length(ref)
        alen = length(alt1)

        # Heuristics for insertion
        isINS = (svtype=="INS" || alt1 ~ /^<INS>$/ || (rlen==1 && alen>1 && alt1 !~ /^[ACGTNacgtn]$/ ? 1:0))
        if (isINS) {
          end = pos
        } else if (rlen==1 && alen==1) {
          # SNP
          end = pos
        } else {
          # DEL/MNV/other without END: span consumed bases
          end = pos + rlen - 1
          if (end < bedStart) end = bedStart  # safety
        }
      }

      if (SAMPLE != "") {
        print chrom, bedStart, end, SAMPLE
      } else {
        print chrom, bedStart, end
      }
    }
  '
}

convert_one() {
  local in="$1" out="$2" sample="$3"
  if [[ "$in" == *.gz ]]; then
    if command -v zcat >/dev/null 2>&1; then
      zcat -- "$in" | vcf_stream_to_bed "$sample" > "$out"
    else
      gzip -cd -- "$in" | vcf_stream_to_bed "$sample" > "$out"
    fi
  else
    cat -- "$in" | vcf_stream_to_bed "$sample" > "$out"
  fi
  echo "✅ Wrote $out"
}

# --- Per-file mode ---
if [[ -n "$VCF_IN" ]]; then
  [[ -f "$VCF_IN" ]] || { echo "❌ Input not found: $VCF_IN"; exit 1; }

  if [[ -n "$OUT_OPT" ]]; then
    if [[ "$OUT_OPT" == *.bed ]]; then
      OUT_PATH="$OUT_OPT"
      mkdir -p "$(dirname "$OUT_PATH")"
    else
      OUT_PATH="${OUT_OPT}.bed"
      mkdir -p "$(dirname "$OUT_PATH")"
    fi
  else
    base="$(basename "$VCF_IN")"
    base="${base%.gz}"; base="${base%.vcf}"
    OUT_PATH="${base}.bed"
  fi

  echo "🔄 Converting $VCF_IN → $OUT_PATH"
  convert_one "$VCF_IN" "$OUT_PATH" "$SAMPLE"
  exit 0
fi

# --- Batch mode (current dir) ---
shopt -s nullglob
mapfile -t FILES < <(printf "%s\n" *.vcf *.vcf.gz 2>/dev/null || true)
if [[ ${#FILES[@]} -eq 0 ]]; then
  echo "❌ No .vcf or .vcf.gz files in $(pwd)"; exit 1
fi

OUT_DIR=""
if [[ -n "$OUT_OPT" ]]; then
  OUT_DIR="$OUT_OPT"
  mkdir -p "$OUT_DIR"
fi

for f in "${FILES[@]}"; do
  base="${f%.gz}"; base="${base%.vcf}"
  if [[ -n "$OUT_DIR" ]]; then
    out="${OUT_DIR%/}/${base}.bed"
  else
    out="${base}.bed"
  fi
  echo "🔄 Converting $f → $out"
  convert_one "$f" "$out" "$SAMPLE"
done

echo "✅ Done."
