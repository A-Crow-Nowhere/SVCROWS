#!/bin/bash

# Exit on error and undefined variables
set -euo pipefail

# === Default Config ===
LIST_FILE=""
OUT_PREFIX="./sample"




# === Parse Arguments ===
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        -L) LIST_FILE="$2"; shift ;;
        -O) OUT_PREFIX="$2"; shift ;;
        -h|--help)


                       echo "Usage: $0 [options

By default, will output a file called sample.merged.crows.tsv


             --Required--
  -L,                        list of files that are in the SVCROWS input format, e.g.:
  								
  							  # Files to merge that 
							  /path/to/sample1.tsv
							  /path/to/sample2.tsv
							  ./sample3.tsv
							  relative/path/to/sample4.tsv
  

             --Optional--
  -O, --output-prefix        by default is './sample'

              --Misc.--
  -h, --help                    Show this help message"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
    shift
done

	
> "$OUT_PREFIX.merged.crows.tsv" 

first=1
while IFS= read -r file || [ -n "$file" ]; do
    [[ -z "$file" || "$file" == \#* ]] && continue

    if [[ $first -eq 1 ]]; then
        # Include header from the first file
        cat "$file" >> "$OUT_PREFIX.merged.crows.tsv"
        first=0
    else
        # Skip header (assume first line is the header)
        tail -n +2 "$file" >> "$OUT_PREFIX.merged.crows.tsv"
    fi
done < "$LIST_FILE"
