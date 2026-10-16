#!/usr/bin/env bash
# run_pipeline.sh
# Full pipeline run (all phases).
# Usage: bash run_pipeline.sh <config_file> [snakemake_extra_args]
# Example: bash run_pipeline.sh configs/config_exp1.yml --cores 8

set -euo pipefail

CONFIG="${1:?Usage: $0 <config_file> [snakemake_args]}"
shift
EXTRA_ARGS="${*:-}"

if [[ ! -f "$CONFIG" ]]; then
  echo "ERROR: Config file not found: $CONFIG" >&2
  exit 1
fi

# Save a config snapshot before run
OUTPUT_DIR=$(python3 -c "import yaml,sys; c=yaml.safe_load(open('$CONFIG')); print(c.get('output_dir','output'))")
mkdir -p "$OUTPUT_DIR"
cp "$CONFIG" "$OUTPUT_DIR/config_used.yml"
echo "[run_pipeline] Config snapshot saved to $OUTPUT_DIR/config_used.yml"

echo "[run_pipeline] Starting pipeline with config: $CONFIG"
echo "[run_pipeline] Output dir: $OUTPUT_DIR"

snakemake \
  --configfile "$CONFIG" \
  --cores 4 \
  --rerun-incomplete \
  --printshellcmds \
  $EXTRA_ARGS

echo "[run_pipeline] Pipeline complete. Results in: $OUTPUT_DIR"
