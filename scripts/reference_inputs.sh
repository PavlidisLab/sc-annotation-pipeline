#!/usr/bin/env bash
# Write or verify SHA-256 hashes of the files that define the Census reference build.
# Usage: scripts/reference_inputs.sh write | check
# Some inputs (the large author obs.tsv files) are gitignored; the hashes in assets/REFERENCE_INPUTS.sha256 are tracked.
set -euo pipefail
cd "$(dirname "$0")/.."
MANIFEST=assets/REFERENCE_INPUTS.sha256
files() {
  ls assets/author_cell_annotations/*.obs.tsv assets/original_celltype_columns.tsv assets/rename_cells_*.tsv | LC_ALL=C sort
}
case "${1:-}" in
  write) files | xargs sha256sum > "$MANIFEST"; echo "wrote $PWD/$MANIFEST" ;;
  check) sha256sum --check --quiet "$MANIFEST" && echo "reference inputs match $PWD/$MANIFEST" ;;
  *) echo "usage: $0 write|check" >&2; exit 2 ;;
esac
