#!/usr/bin/env bash
# Times degibbs at several thread counts and writes a CSV: wall seconds and peak memory per run.
#
#   bench/run.sh [-b BINARY] [-r REFERENCE_BINARY] [-t "1 2 4 8 16"] [-o results.csv] IMAGE.nii[.gz] [IMAGE ...]
#
# IMAGE is a magnitude image with a sidecar next to it (or with a known factor: set PF, e.g. PF="--pf 6/8").
# For each image and thread count it runs `degibbs IMAGE OUT --threads N --overwrite` and records wall time and
# peak resident memory. With -r it also runs the reference build once per image (all threads) and reports the
# largest difference and the number of differing voxels (needs python3 with numpy and nibabel).
# Not part of CI; meant for a workstation. Prints a summary table as well.
set -euo pipefail
BIN=target/release/degibbs; REF=; THREADS="1 2 4 8 16"; OUT=bench-results.csv
while getopts "b:r:t:o:" o; do
  case $o in b) BIN=$OPTARG;; r) REF=$OPTARG;; t) THREADS=$OPTARG;; o) OUT=$OPTARG;; *) exit 2;; esac
done
shift $((OPTIND - 1))
[ $# -ge 1 ] || { sed -n '2,11p' "$0"; exit 2; }
PF=${PF:-}
case "$(uname)" in
  Darwin) TIME=(/usr/bin/time -l); peak() { awk '/maximum resident set size/ {printf "%.0f", $1/1048576}' "$1"; };;
  *)      TIME=(/usr/bin/time -v); peak() { awk '/Maximum resident set size/ {printf "%.0f", $6/1024}' "$1"; };;
esac
work=$(mktemp -d); trap 'rm -rf "$work"' EXIT
echo "image,threads,wall_seconds,peak_mb" > "$OUT"
for img in "$@"; do
  ext=.nii; case $img in *.nii.gz) ext=.nii.gz;; esac
  for n in $THREADS; do
    start=$(date +%s.%N 2>/dev/null || python3 -c 'import time;print(time.time())')
    # shellcheck disable=SC2086
    "${TIME[@]}" -o "$work/time.txt" "$BIN" "$img" "$work/out$ext" --threads "$n" --overwrite $PF 2> "$work/err.txt"
    end=$(date +%s.%N 2>/dev/null || python3 -c 'import time;print(time.time())')
    wall=$(awk -v a="$start" -v b="$end" 'BEGIN{printf "%.2f", b-a}')
    echo "$(basename "$img"),$n,$wall,$(peak "$work/time.txt")" | tee -a "$OUT"
  done
  if [ -n "$REF" ]; then
    # shellcheck disable=SC2086
    "$REF" "$img" "$work/ref$ext" $PF 2> /dev/null
    "$BIN" "$img" "$work/new$ext" --overwrite $PF 2> /dev/null
    python3 - "$work/ref$ext" "$work/new$ext" <<'PY'
import sys, numpy as np, nibabel as nib
a, b = (np.asanyarray(nib.load(p).dataobj, dtype=np.float32) for p in sys.argv[1:3])
d = np.abs(a - b)
print(f"vs reference: largest difference {d.max():.3e}, {int((d > 0).sum())} of {d.size} voxels differ")
PY
  fi
done
echo "wrote $OUT"
