#!/usr/bin/env bash
# H3c: re-run the H3 primary/H3a (script 05) and H3b (script 06) analyses on Sanger Project Score (Chronos-processed)
# Outputs: results/07ps_*  Usage: bash scripts/07_project_score_h3.sh <project_root>
set -euo pipefail
ROOT="${1:-.}"
python3 "$ROOT/scripts/05_h3_primary_h3a.py" "$ROOT" projectscore_chronos_gene_effect.csv projectscore_chronos_common_essentials.csv 07ps > "$ROOT/results/07ps_h3_log.txt"
python3 "$ROOT/scripts/06_h3b_activity_vs_dependency.py" "$ROOT" projectscore_chronos_gene_effect.csv 07ps > "$ROOT/results/07ps_h3b_log.txt"
