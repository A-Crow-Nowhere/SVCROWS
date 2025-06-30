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
        -C|--vcf-file) BREAK_CALL="$2"; shift ;;
        -s) SAMPLE_NAME="$2"; shift ;;
        -O) OUT_PREFIX="$2"; shift ;;
        -h|--help)
        
        
       echo "Usage: $0 [options]

  By default will output a file named 'sample.crows.vcf'
  This will automatically select for features that have a definite
  Start and Stop site - and only those. 

             --Required--
  -C, --crow-file             crow.FCL (merged output) formatted file

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


	


if [[ -z "$BREAK_CALL" ]]; then
  echo "Usage: bash convert_structural_to_vcf.sh <input.tsv> <output.vcf>"
  exit 1
fi

# Write VCF headers
cat <<EOF > "$OUT_PREFIX.crows.vcf"
##fileformat=VCFv4.2
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=INS,Description="Insertion">
##ALT=<ID=CNV,Description="Copy number variation">
##ALT=<ID=BND,Description="Breakend">
##ALT=<ID=TRA,Description="Translocation">
##ALT=<ID=CTX,Description="Inter-chromosomal translocation">
##ALT=<ID=ITX,Description="Intra-chromosomal translocation">
##ALT=<ID=CPX,Description="Complex variant">
##ALT=<ID=MCNV,Description="Multiallelic CNV">
##ALT=<ID=MEI,Description="Mobile element insertion">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variant">
##INFO=<ID=END,Number=1,Type=Integer,Description="End coordinate of the variant">
##INFO=<ID=IsKnown,Number=1,Type=String,Description="Whether variant is known">
##INFO=<ID=NumReads,Number=1,Type=Integer,Description="Number of supporting reads">
##INFO=<ID=ROPercentPass,Number=1,Type=Integer,Description="Read overlap % passing filters">
##INFO=<ID=ROCount,Number=1,Type=Integer,Description="Read overlap count">
##INFO=<ID=BPStartLeft,Number=1,Type=Integer,Description="Breakpoint start (left)">
##INFO=<ID=BPStartRight,Number=1,Type=Integer,Description="Breakpoint start (right)">
##INFO=<ID=BPStartCount,Number=1,Type=Integer,Description="Breakpoint start count">
##INFO=<ID=BPEndLeft,Number=1,Type=Integer,Description="Breakpoint end (left)">
##INFO=<ID=BPEndRight,Number=1,Type=Integer,Description="Breakpoint end (right)">
##INFO=<ID=BPEndCount,Number=1,Type=Integer,Description="Breakpoint end count">
##INFO=<ID=Abundance,Number=1,Type=Integer,Description="Read abundance estimate">
##INFO=<ID=Rarity,Number=1,Type=String,Description="Rarity class">
##INFO=<ID=Matches,Number=1,Type=String,Description="Matching variant IDs">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
EOF

	# Append data
	tail -n +2 "$BREAK_CALL" | awk -F'\t' 'BEGIN {OFS="\t"} {
	    chrom = $1;
	    pos = $2;
	    end = $3;
	    len = $4;
	    type = $5;
	    id = $6;
	    info = $7;
	    format_field = "GT";
	    sample_value = "./.";

	    # Columns 10–25
	    isKnown = $10;
	    qscore = $11;
	    numReads = $12;
	    roPercent = $15;
	    roCount = $16;
	    bpStartLeft = $17;
	    bpStartRight = $18;
	    bpStartCount = $19;
	    bpEndLeft = $20;
	    bpEndRight = $21;
	    bpEndCount = $22;
	    abundance = $23;
	    rarity = $24;
	    matches = $25;

	    

	    # Add extended info fields
	    info = info ";IsKnown=" isKnown ";NumReads=" numReads;
	    info = info ";ROPercentPass=" roPercent ";ROCount=" roCount;
	    info = info ";BPStartLeft=" bpStartLeft ";BPStartRight=" bpStartRight ";BPStartCount=" bpStartCount;
	    info = info ";BPEndLeft=" bpEndLeft ";BPEndRight=" bpEndRight ";BPEndCount=" bpEndCount;
	    info = info ";Abundance=" abundance ";Rarity=" rarity ";Matches=" matches;

	    # Print final VCF line
	    print chrom, pos, id, "N", type,  qscore, $9, info, format_field, sample_value;
	}' >> "$OUT_PREFIX.crows.vcf"

	echo "✅ VCF written to $OUT_PREFIX.crows.vcf"
