#!/usr/bin/env bash
# run_downstream_only.sh
# Skip input loading (assume output/loaded/seurat_validated.rds exists).
# Useful for re-running DE / enrichment with different cutoffs.
# Usage: bash run_downstream_only.sh <config_file> [snakemake_extra_args]

set -euo pipefail

CONFIG="${1:?Usage: $0 <config_file> [snakemake_args]}"
shift
EXTRA_ARGS="${*:-}"

if [[ ! -f "$CONFIG" ]]; then
  echo "ERROR: Config file not found: $CONFIG" >&2
  exit 1
fi

OUTPUT_DIR=$(python3 -c "import yaml,sys; c=yaml.safe_load(open('$CONFIG')); print(c.get('output_dir','output'))")

if [[ ! -f "$OUTPUT_DIR/loaded/seurat_validated.rds" ]]; then
  echo "ERROR: Validated RDS not found at $OUTPUT_DIR/loaded/seurat_validated.rds" >&2
  echo "Run run_pipeline.sh first to generate it." >&2
  exit 1
fi

# Update config snapshot
mkdir -p "$OUTPUT_DIR"
cp "$CONFIG" "$OUTPUT_DIR/config_used.yml"

echo "[run_downstream] Starting downstream-only run with config: $CONFIG"

snakemake \
  --configfile "$CONFIG" \
  --cores 4 \
  --rerun-incomplete \
  --printshellcmds \
  --ignore-incomplete \
  $EXTRA_ARGS

echo "[run_downstream] Done. Results in: $OUTPUT_DIR"
