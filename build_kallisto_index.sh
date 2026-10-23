#!/usr/bin/env bash
set -euo pipefail
if (( $# != 2 )); then
  echo 'Usage: bash build_kallisto_index.sh transcripts.fa[.gz] new_index.idx' >&2
  exit 2
fi
command -v kallisto >/dev/null || { echo 'Missing kallisto. Nothing was installed.' >&2; exit 1; }
[[ -f "$1" ]] || { echo "Missing transcript FASTA: $1" >&2; exit 1; }
[[ ! -e "$2" ]] || { echo "Index already exists: $2; choose a new name." >&2; exit 1; }
kallisto index -i "$2" "$1"
