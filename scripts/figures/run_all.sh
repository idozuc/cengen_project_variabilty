#!/usr/bin/env bash
# Rebuild every manuscript figure panel. Run from the repository root:
#   scripts/figures/run_all.sh config/figures.yaml
# Writes OUTPUT_DIR/data (intermediate tables) and OUTPUT_DIR/panels.
# Optional: PYTHON (default python3).
set -euo pipefail
config="${1:?usage: scripts/figures/run_all.sh CONFIG.yaml}"
output_dir="$("${PYTHON:-python3}" -c 'import sys, yaml; print(yaml.safe_load(open(sys.argv[1]))["output_dir"])' "$config")"

for script in scripts/figures/0*_export_*.R; do
  echo "== $script"
  Rscript "$script" --config "$config"
done

for script in scripts/figures/fig*.py; do
  echo "== $script"
  "${PYTHON:-python3}" "$script" --config "$config"
done

echo "== workflow diagrams (Figures 1-2)"
panels="$(cd "$output_dir" && pwd)/panels"
for name in methods_hvg_workflow methods_clustering_workflow; do
  (cd scripts/figures/workflows && pdflatex -interaction=nonstopmode -output-directory "$panels" "$name.tex" >/dev/null)
done
rm -f "$panels"/*.aux "$panels"/*.log
mv "$panels/methods_hvg_workflow.pdf" "$panels/fig01_methods_hvg_workflow.pdf"
mv "$panels/methods_clustering_workflow.pdf" "$panels/fig02_methods_clustering_workflow.pdf"
