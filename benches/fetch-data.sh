#!/usr/bin/env bash
# Idempotent fetch of the FIRE test-data into benches/data/fire-test-data/.
# The bucket is public over anonymous HTTPS, so no AWS CLI is required.
# Source: https://github.com/fiberseq/FIRE pixi.toml `test-data` task.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEST="$SCRIPT_DIR/data/fire-test-data"
BASE="https://s3.kopah.orci.washington.edu/stergachis/public/FIRE/test-data"

mkdir -p "$DEST"

FILES=(
  test.cram
  test.cram.crai
  test.fa.gz
  test.fa.gz.fai
  test.fa.gz.gzi
  test.tbl
  test.yaml
)

for f in "${FILES[@]}"; do
  if [ -s "$DEST/$f" ]; then
    echo "skip  $f (already present)"
  else
    echo "fetch $f"
    curl --fail --location --silent --show-error --output "$DEST/$f" "$BASE/$f"
  fi
done

echo
echo "Dataset ready at: $DEST"
du -sh "$DEST"
