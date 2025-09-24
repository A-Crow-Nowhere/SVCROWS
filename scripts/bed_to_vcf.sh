#!/usr/bin/env bash
set -euo pipefail

# bedtovcf.sh — Convert 3-col BED/TSV to minimal SV-VCF.
# Modes:
#   1) Recursive:  bedtovcf.sh -r ref.fa -i in_root -o out_root
#   2) Per-file:   bedtovcf.sh -r ref.fa in.bed out.vcf.gz
#
# Input: CHROM  START  END  (BED 0-based, END right-open)
# Output: Minimal CNV-style VCF (BGZF-compressed + .tbi index)

REF="${HOME}/genomes/hg38.fa"
INROOT="."
OUTROOT="./vcf_files"

print_help() {
  cat <<EOF
Usage:
  Recursive mode:
    $0 -r ref.fa -i <input_root> -o <output_root>

  Per-file mode:
    $0 -r ref.fa <in.bed|.tsv[.gz]> <out.vcf.gz>

Options:
  -r, --ref   Reference FASTA (default: ${HOME}/genomes/hg38.fa)
  -i, --in    Input root directory for recursive scan (default: .)
  -o, --out   Output root directory (default: ./vcf_files)

Notes:
  - BED/TSV must have at least 3 columns: CHROM START END
  - START is 0-based; converted to 1-based POS in VCF
EOF
}

# ---------- parse args ----------
while [[ $# -gt 0 ]]; do
  case "$1" in
    -r|--ref) REF="$2"; shift 2;;
    -i|--in)  INROOT="$2"; shift 2;;
    -o|--out) OUTROOT="$2"; shift 2;;
    -h|--help) print_help; exit 0;;
    --) shift; break;;
    -*) echo "Unknown option: $1" >&2; exit 2;;
    *)  break;;
  esac
done

_expand() { printf '%s' "${1/#\~/$HOME}"; }
REF=$(_expand "$REF"); INROOT=$(_expand "$INROOT"); OUTROOT=$(_expand "$OUTROOT")
if command -v realpath >/dev/null 2>&1; then
  REF="$(realpath -m "$REF")"
  INROOT="$(realpath -m "$INROOT")"
  OUTROOT="$(realpath -m "$OUTROOT")"
fi

[[ -r "$REF" ]] || { echo "❌ Reference FASTA not found/readable: $REF" >&2; exit 1; }

# ---------- header from reference ----------
FAI="${REF}.fai"
if [[ ! -s "$FAI" ]]; then
  echo "Indexing reference FASTA: $REF"
  samtools faidx "$REF"
fi

CONTIG_LIST="$(mktemp)"
cut -f1 "$FAI" > "$CONTIG_LIST"

HEADER_TXT="$(mktemp)"
{
  echo "##fileformat=VCFv4.2"
  echo "##source=bed_to_vcf_minimal"
  awk 'BEGIN{FS="\t"} {printf "##contig=<ID=%s,length=%s>\n",$1,$2}' "$FAI"
  echo '##ALT=<ID=CNV,Description="Copy number variant">'
  echo '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant (1-based, inclusive)">'
  echo '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">'
  echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
} > "$HEADER_TXT"

convert_one() {
  local infile="$1" outpath="$2"
  [[ -r "$infile" ]] || { echo "❌ Cannot read: $infile" >&2; return 1; }

  local -a reader=(cat --)
  [[ "$infile" == *.gz ]] && reader=(zcat --)

  # 1) header  2) body from BED→VCF  3) sort body by CHROM,POS  4) compress+index
  {
    # header (unsorted; printed as-is)
    cat "$HEADER_TXT"

    # body (sorted)
    awk -F'\t' -v OFS='\t' '
      NR==FNR { ok[$1]=1; next }       # load contigs from CONTIG_LIST
      $0 ~ /^#/ { next }               # skip comments
      NF >= 3 && $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ {
        chrom=$1; start=$2; end=$3;
        if (!(chrom in ok)) next;      # only contigs in reference .fai
        pos = start + 1;               # 0-based -> 1-based
        if (pos < 1 || end < pos) next;
        info = "END=" end ";SVTYPE=CNV";
        print chrom, pos, ".", "N", "<CNV>", ".", "PASS", info;
      }
    ' "$CONTIG_LIST" <("${reader[@]}" "$infile") \
    | sort -k1,1 -k2,2n
  } | bgzip -c > "$outpath"

  tabix -f -p vcf "$outpath"
}

# ---------- decide mode ----------
if (( $# == 0 )); then
  # Recursive mode (restricted to */per_cell/*)
  mkdir -p "$OUTROOT"
  echo "🔎 Scanning: $INROOT  →  $OUTROOT"
  while IFS= read -r -d '' infile; do
    rel="${infile#$INROOT/}"
    outdir="${OUTROOT}/$(dirname "$rel")"
    mkdir -p "$outdir"
    base="$(basename "$infile")"
    base="${base%.gz}"; base="${base%.bed}"; base="${base%.tsv}"
    outpath="${outdir}/${base}.vcf.gz"
    echo "🧪 Converting: $infile"
    convert_one "$infile" "$outpath"
    echo "✅ Wrote: $outpath"
  done < <(
    find "$INROOT" -type f -path '*/per_cell/*' \
      \( -name '*.bed' -o -name '*.bed.gz' -o -name '*.tsv' -o -name '*.tsv.gz' \) \
      -print0
  )
elif (( $# == 2 )); then
  # Per-file mode
  in="$1"; out="$2"
  echo "🧪 Converting: $in"
  convert_one "$in" "$out"
  echo "✅ Wrote: $out"
else
  echo "Unknown arg(s): $*" >&2
  print_help >&2
  exit 2
fi

# ---------- cleanup ----------
rm -f "$CONTIG_LIST" "$HEADER_TXT"
