#!/usr/bin/env bash
set -euo pipefail

# Run while your main environment (normally rnaseq) is active.
# RSeQC is checked separately with conda run; environments are not activated here.
# Usage: bash check_dependencies.sh [rnaseq|rmats]
# Optional: RSEQC_ENV=another_name bash check_dependencies.sh rnaseq
if (( $# > 1 )); then
  echo 'Usage: bash check_dependencies.sh [rnaseq|rmats]' >&2
  exit 2
fi
mode=${1:-rnaseq}
RSEQC_ENV=${RSEQC_ENV:-rseqc}
missing=0

# These commands are checked in the currently active environment.
case "$mode" in
  rnaseq) commands=(python3 Rscript STAR samtools fastqc multiqc kallisto) ;;
  rmats) commands=(python3) ;;
  *) echo 'Usage: bash check_dependencies.sh [rnaseq|rmats]' >&2; exit 2 ;;
esac

for cmd in "${commands[@]}"; do
  if command -v "$cmd" >/dev/null; then
    echo "OK: $cmd"
  else
    echo "MISSING: $cmd"
    missing=1
  fi
done

# Use the same Rscript as the main pipeline, preserving your system R setup.
if [[ "$mode" == rnaseq ]] && command -v Rscript >/dev/null; then
  Rscript -e 'p <- c("tximport","DESeq2","fgsea","ggplot2","ggrepel","SummarizedExperiment"); ok <- vapply(p, requireNamespace, logical(1), quietly=TRUE); print(ok); if (!all(ok)) quit(status=1)' || missing=1
fi

# RSeQC lives in its own environment, not in the active rnaseq environment.
if [[ "$mode" == rnaseq ]]; then
  if command -v conda >/dev/null; then
    echo "Checking RSeQC in Conda environment: $RSEQC_ENV"
    rseqc_commands=(infer_experiment.py read_distribution.py geneBody_coverage.py junction_annotation.py)
    for cmd in "${rseqc_commands[@]}"; do
      if conda run --no-capture-output -n "$RSEQC_ENV" "$cmd" --help >/dev/null; then
        echo "OK: $cmd (environment: $RSEQC_ENV)"
      else
        echo "FAILED: $cmd (environment: $RSEQC_ENV; check that the environment and tool are installed)" >&2
        missing=1
      fi
    done
    # RSeQC also invokes Rscript to produce some of its QC plots.
    if conda run --no-capture-output -n "$RSEQC_ENV" Rscript --help >/dev/null; then
      echo "OK: Rscript for RSeQC plots (environment: $RSEQC_ENV)"
    else
      echo "FAILED: Rscript for RSeQC plots (environment: $RSEQC_ENV)" >&2
      missing=1
    fi
  else
    echo 'MISSING: conda executable; cannot check the separate RSeQC environment.' >&2
    missing=1
  fi
fi

# Keep the original rMATS check: custom wrappers are configured in config.json.
if [[ "$mode" == rmats ]]; then
  echo 'Check the rmats_command configured in config.json in this environment.'
  if command -v rmats.py >/dev/null; then
    echo 'OK: rmats.py'
  else
    echo 'rmats.py is not on PATH; configure its absolute path or an installed wrapper.'
  fi
fi

exit "$missing"
