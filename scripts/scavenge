#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="svcrows-env"

if [[ $# -lt 2 ]] || [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  cat <<'USAGE'
Usage:
  scavenge INPUT_QUERY_LIST OUTPUT_DIRECTORY [ExpandRORegion] [BPfactor] [DefaultSizes] [xs] [xl] [y1s] [y1l] [y2s] [y2l]

Notes:
  • Booleans: TRUE/FALSE (also yes/no, 1/0).
  • Omit trailing args to use SVCROWS defaults.
Examples:
  scavenge ./in ./out TRUE TRUE TRUE
  scavenge ./in ./out FALSE TRUE FALSE 5000 25000 500 2500 50 80
USAGE
  exit 0
fi

# Stream R output directly (no conda buffering)
conda run --no-capture-output -n "$ENV_NAME" Rscript - "$@" <<'RS'
options(
  warn = 1,
  error = function(e){ message("R ERROR: ", conditionMessage(e)); traceback(2); quit(status=1) }
)

args <- commandArgs(trailingOnly=TRUE)

# helpers for type conversion
b <- function(x){
  if (length(x)==0) return(NULL)
  s <- tolower(as.character(x))
  if (s %in% c("true","t","1","yes","y")) return(TRUE)
  if (s %in% c("false","f","0","no","n")) return(FALSE)
  if (s %in% c("na","")) return(NA)
  as.logical(x)
}
n <- function(x){
  if (length(x)==0) return(NULL)
  s <- tolower(as.character(x))
  if (s %in% c("na","")) return(NA_real_)
  as.numeric(x)
}

if (length(args) < 2) stop("Need at least INPUT_QUERY_LIST and OUTPUT_DIRECTORY")

params <- list(
  InputQueryList = args[[1]],
  OutputDirectory = args[[2]]
)
if (length(args) >= 3)  params$ExpandRORegion <- b(args[[3]])
if (length(args) >= 4)  params$BPfactor       <- b(args[[4]])
if (length(args) >= 5)  params$DefaultSizes   <- b(args[[5]])
if (length(args) >= 6)  params$xs             <- n(args[[6]])
if (length(args) >= 7)  params$xl             <- n(args[[7]])
if (length(args) >= 8)  params$y1s            <- n(args[[8]])
if (length(args) >= 9)  params$y1l            <- n(args[[9]])
if (length(args) >= 10) params$y2s            <- n(args[[10]])
if (length(args) >= 11) params$y2l            <- n(args[[11]])

# Ensure output dir exists before calling
odir <- params$OutputDirectory
dir.create(odir, recursive=TRUE, showWarnings=FALSE)
if (!dir.exists(odir)) stop("Cannot create OutputDirectory: ", odir)

# Resolve function directly from namespace (avoids :: resolution weirdness)
if (!requireNamespace("SVCROWS", quietly=TRUE))
  stop("Package 'SVCROWS' not installed in this env")

fn <- tryCatch(
  getExportedValue("SVCROWS","Scavenge"),
  error = function(e) get("Scavenge", envir=asNamespace("SVCROWS"))
)
if (!is.function(fn)) stop("'Scavenge' not found in SVCROWS")

# Log some context
cat("pkg path: ", system.file(package="SVCROWS"), "\n", file=stderr())
cat("🔧 Calling SVCROWS::Scavenge...\n", file=stderr()); flush.console()

# Actual call — your prints inside Scavenge should appear now
invisible(do.call(fn, params))

cat("✅ Done\n", file=stderr()); flush.console()
RS
