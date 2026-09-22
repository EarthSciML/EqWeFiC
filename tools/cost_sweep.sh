#!/bin/bash
# Re-measure the per-directory test cost table in CLAUDE.md.
#
# Measures in a PINNED git worktree, not in the working checkout, because three
# agents commit to this repo every few minutes and a sweep that takes half an
# hour will otherwise measure four different repositories.  The 2026-09-21 sweep
# that this replaced reported 903 assertions for wrf_wesely at 18:26 and 2048
# for the same directory at 19:30 -- the difference was a commit at 18:28, not
# the runner.
#
# Reports WALL, CPU and PEAK RSS per directory.  Trust the CPU number: the host
# is a 20-core Slurm allocation shared with other agents, and wall time has been
# seen to vary 1.8x with load on the same document.
#
# usage:  tools/cost_sweep.sh [output-dir]
set -u
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
OUT=${1:-/tmp/esm-cost-$$}
SNAP=${REPO}-costsnap.$$   # unique per run: two concurrent sweeps must not share one
ESM=$REPO/esm
export ESD_ROOT=${ESD_ROOT:-$(dirname "$REPO")/EarthSciDiscretizations-integ}   # carries ESD PR 42, which wsm6/sedimentation.esm imports

SHA=$(git -C "$REPO" rev-parse --short HEAD)
git -C "$REPO" worktree prune
git -C "$REPO" worktree add --detach "$SNAP" "$SHA" >/dev/null 2>&1 || exit 1
mkdir -p "$OUT"
cd "$SNAP" || exit 1

DIRS="lib
components/land_surface
components/atmospheric_deposition/wrf_wesely
components/gaschem
components/wildland_fire
components/atmospheric_dynamics/ysu
components/atmospheric_dynamics/sfclayrev
components/atmospheric_dynamics/wsm6
components/atmospheric_dynamics/wrf_arw
components/atmospheric_radiation/rrtm_lw
components/atmospheric_radiation/dudhia_sw
components/atmospheric_radiation/wrf_solar
components/atmospheric_radiation/cloud_fraction
components
couplings"

printf '%-46s %6s %5s %5s %8s %8s %8s\n' "directory @ $SHA" pass fail err "wall s" "cpu s" "peak MB"
for d in $DIRS; do
  n=$(echo "$d" | tr '/' '_')
  /usr/bin/time -f "%e %U %S %M" "$ESM" test "$d" > "$OUT/$n.out" 2> "$OUT/$n.time"
  read -r w u s m <<<"$(tail -1 "$OUT/$n.time")"
  read -r _ p f e <<<"$(grep -E '^  TOTAL' "$OUT/$n.out" | tr -s ' ')"
  printf '%-46s %6s %5s %5s %8s %8s %8.0f\n' "$d" "${p:-?}" "${f:-?}" "${e:-?}" "$w" \
         "$(echo "$u $s" | awk '{printf "%.1f", $1+$2}')" "$((m / 1024))"
done
echo "(load average now $(cut -d' ' -f1-3 /proc/loadavg))"
git -C "$REPO" worktree remove --force "$SNAP" 2>/dev/null
