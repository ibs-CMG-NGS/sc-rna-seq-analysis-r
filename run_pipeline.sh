#!/usr/bin/env bash
# run_pipeline.sh
# Deprecated shim. Use Snakemake directly:
#   snakemake --profile profiles/default --configfile <config> --cores <N>
#
# Usage: bash run_pipeline.sh <config_file> [snakemake_extra_args]

set -euo pipefail

CONFIG="${1:?Usage: $0 <config_file> [snakemake_args]}"
shift

if [[ ! -f "$CONFIG" ]]; then
  echo "ERROR: Config file not found: $CONFIG" >&2
  exit 1
fi

echo "WARNING: run_pipeline.sh is deprecated; use snakemake --profile profiles/default ..." >&2

exec snakemake \
  --profile profiles/default \
  --configfile "$CONFIG" \
  "$@"
