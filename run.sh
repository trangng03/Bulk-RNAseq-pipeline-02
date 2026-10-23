#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
command -v python3 >/dev/null || { echo "Missing python3; nothing was installed." >&2; exit 1; }
exec python3 "$SCRIPT_DIR/scripts/pipeline.py" "$@"
