#!/usr/bin/env bash
set -euo pipefail
ENV_NAME="${ENV_NAME:-svcrows-env}"

if [[ $# -lt 3 ]] || [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  cat <<'USAGE'
Usage:
  hunt INPUT_QUERY_LIST FEATURE_LIST OUTPUT_DIRECTORY [BPfactor] [DefaultSizes] [xs] [xl] [y1s] [y1l] [y2s] [y2l]

Notes:
  • Booleans: TRUE/FALSE (also yes/no, 1/0). Use NA for missing.
  • 'DefaultSizes' is the correct spelling (wrapper also accepts legacy 'DefualtSizes').
Examples:
  hunt ./in ./Featurelist.tsv ./out TRUE FALSE 3000 6000 300 600 30 60
  hunt ./in ./Featurelist.tsv ./out TRUE TRUE
USAGE
  exit 0
fi

conda run --no-capture-output -n "$ENV_NAME" Rscript - "$@" <<'RS'
options(
  warn = 1,
  error = function(e){ message("R ERROR: ", conditionMessage(e)); traceback(2); quit(status=1) }
)

args <- commandArgs(trailingOnly=TRUE)

# helpers
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
  suppressWarnings(as.numeric(x))
}

if (length(args) < 3)
  stop("Need INPUT_QUERY_LIST, FEATURE_LIST, OUTPUT_DIRECTORY")

params <- list(
  InputQueryList  = args[[1]],
  FeatureList     = args[[2]],
  OutputDirectory = args[[3]]
)
if (length(args) >= 4)  params$BPfactor     <- b(args[[4]])
if (length(args) >= 5)  params$DefaultSizes <- b(args[[5]])   # <- correct name
if (length(args) >= 6)  params$xs           <- n(args[[6]])
if (length(args) >= 7)  params$xl           <- n(args[[7]])
if (length(args) >= 8)  params$y1s          <- n(args[[8]])
if (length(args) >= 9)  params$y1l          <- n(args[[9]])
if (length(args) >= 10) params$y2s          <- n(args[[10]])
if (length(args) >= 11) params$y2l          <- n(args[[11]])

# Back-compat: if someone passes the legacy mis-spelling as a named arg,
# normalize it to DefaultSizes so Hunt() accepts it.
if (!is.null(params$DefualtSizes) && is.null(params$DefaultSizes)) {
  params$DefaultSizes <- params$DefualtSizes
  params$DefualtSizes <- NULL
}

# Make sure output dir exists
odir <- params$OutputDirectory
dir.create(odir, recursive=TRUE, showWarnings=FALSE)
if (!dir.exists(odir)) stop("Cannot create OutputDirectory: ", odir)

if (!requireNamespace("SVCROWS", quietly=TRUE))
  stop("Package 'SVCROWS' not installed in this env")

fn <- tryCatch(
  getExportedValue("SVCROWS","Hunt"),
  error = function(e) get("Hunt", envir=asNamespace("SVCROWS"))
)
if (!is.function(fn)) stop("'Hunt' not found in SVCROWS")

cat("pkg path: ", system.file(package="SVCROWS"), "\n", file=stderr())
cat("🔎 Hunt params:\n", file=stderr()); str(params)
cat("🔎 Calling SVCROWS::Hunt...\n", file=stderr()); flush.console()

invisible(do.call(fn, params))
cat("✅ Done\n", file=stderr()); flush.console()
RS
