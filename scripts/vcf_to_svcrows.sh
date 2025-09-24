#!/bin/bash

# Exit on error and undefined variables
set -euo pipefail

# === Default Config ===
BREAK_CALL=""
OUT_PREFIX="./sample"
SAMPLE_NAME="Sample"

# === Parse Arguments ===
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        -V|--vcf-file) BREAK_CALL="$2"; shift ;;
        -s) SAMPLE_NAME="$2"; shift ;;
        -O) OUT_PREFIX="$2"; shift ;;
        -h|--help)
        
        
                       echo "Usage: $0 [options]

  By default will output a file named 'sample.crows.tsv'
  This will automatically select for features that have a definite
  Start and Stop site - and only those. 

             --Required--
  -V, --vcf-file             VCF formatted file

             --Optional--
  -O         Output prefix   by default is './sample'
  -s         sample name     add a sample name into the Var3 column (used to adjust the output)
  						     essentially, give each sample/run a unique string ID here
  						     eg:'MyCoolSample'
  						     Default: 'Sample'
              --Misc.--  
  -h, --help               	Show this help message"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
    shift
done


if [[ -n "$BREAK_CALL" ]]; then

bcftools view -i '(INFO/SVTYPE=="DEL" || INFO/SVTYPE=="DUP" || INFO/SVTYPE=="DUP:TANDEM" || INFO/SVTYPE=="INVDUP" || INFO/SVTYPE=="INV" || INFO/SVTYPE=="CNV")' "$BREAK_CALL" > filtered.tmp.vcf

bcftools view -H filtered.tmp.vcf | \
awk -v sample_name="$SAMPLE_NAME" -F'\t' 'BEGIN { OFS="\t"; id=1 }
{
    info = $8

    chr = $1
    start = $2
    end = "."; svlength = "."; svtype = "."; svmethod = "."; mapq = "0"; rereads = "0"

    # Parse INFO fields and build filtered Var1 (exclude RNAMES)
    n = split(info, info_fields, ";")
    var1 = ""
    for (i = 1; i <= n; i++) {
        split(info_fields[i], kv, "=")
        key = kv[1]; val = kv[2]
        if (key == "END") end = val
        if (key == "SVLEN") svlength = val
        if (key == "SVTYPE") svtype = val
        if (key == "SVMETHOD") svmethod = val
        if (key == "MAPQ") mapq = val
        if (key == "SU") rereads = val
        if (key != "RNAMES") var1 = var1 info_fields[i] ";"
    }
    sub(/;$/, "", var1)

    # ---- Build Var2 from FORMAT/SAMPLE if present; else set "." ----
    var2 = "."
    if (NF >= 9 && $9 != "" && $9 != ".") {
        format = $9
        sample = (NF >= 10 ? $10 : "")
        if (sample != "" && sample != ".") {
            var2 = ""
            nf = split(format, fmt_keys, ":")
            ns = split(sample, fmt_vals, ":")
            for (j = 1; j <= nf && j <= ns; j++) {
                var2 = var2 fmt_keys[j] "=" fmt_vals[j] ";"
            }
            sub(/;$/, "", var2)
            if (var2 == "") var2 = "."   # safety
        }
    }

    # ---- Fallback length calculation: END - START ----
    if ((svlength == "." || svlength == "") && start ~ /^[0-9]+$/ && end ~ /^[0-9]+$/) {
        svlength = (end + 0) - (start + 0)
    }

    print chr, start, end, svlength, svtype, id, var1, var2, sample_name, "FALSE", mapq, rereads
    id++
}' > "$OUT_PREFIX.crows.nohead.tsv"

echo -e "Chr\tStart\tEnd\tLength\tType\tID\tVar1\tVar2\tVar3\tIsKnown\tQScore\tNumReads" | \
cat - "$OUT_PREFIX.crows.nohead.tsv" > "$OUT_PREFIX.crows.tsv"

rm filtered.tmp.vcf "$OUT_PREFIX.crows.nohead.tsv"

fi
	
