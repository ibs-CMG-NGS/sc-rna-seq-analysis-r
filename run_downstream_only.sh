#!/usr/bin/env bash
# run_downstream_only.sh
# Deprecated shim. Use Snakemake directly:
#   snakemake --profile profiles/default --configfile <config> \
#     --config start_from=loaded --cores <N>
#
# Usage: bash run_downstream_only.sh <config_file> [snakemake_extra_args]

set -euo pipefail

CONFIG="${1:?Usage: $0 <config_file> [snakemake_args]}"
shift

if [[ ! -f "$CONFIG" ]]; then
  echo "ERROR: Config file not found: $CONFIG" >&2
  exit 1
fi

echo "WARNING: run_downstream_only.sh is deprecated; use --config start_from=loaded ..." >&2

exec snakemake \
  --profile profiles/default \
  --configfile "$CONFIG" \
  "$@" \
  --config start_from=loaded
