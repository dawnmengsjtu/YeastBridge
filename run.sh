#!/usr/bin/env bash
# No arguments: re-export the submitted candidates. Other modes are explicit.
set -euo pipefail
cd "$(dirname "$0")"
exec "${PY:-python}" main.py "$@"
