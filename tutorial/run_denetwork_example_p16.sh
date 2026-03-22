#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
APP_DIR="${REPO_DIR}/DENetwork"
PYTHON_BIN="${REPO_DIR}/.conda-env/bin/python"

cd "${APP_DIR}"

exec "${PYTHON_BIN}" run_denetwork.py \
  -n influenza_p16 \
  -d example_deseq2_output/DE_pos_influenza.tsv \
  -g example_deseq2_output/gene_log2fc_influenza.tsv \
  -r example/receptors_influenza.txt \
  -t de \
  -s 100 \
  -P 16 \
  "$@"
