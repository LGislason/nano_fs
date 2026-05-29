#!/usr/bin/env bash

set -eu

NP="${NP:-4}"
CTL="${CTL:-nanorod_wu_field_snapshot.ctl}"

# Snapshot defaults. Override these on the command line for smoke tests or
# final high-resolution runs.
RES="${RES:-300}"
GAP="${GAP:-0.006}"
NBG="${NBG:-1.28}"
THETA="${THETA:-80}"
RUN_TIME="${RUN_TIME:-300}"
POL="${POL:-0}"
CASES="${CASES:-40:0.672 120:0.950}"
LOGFILE="${LOGFILE:-run_field_snapshots_$(date +%Y%m%d_%H%M%S).out}"

if [ ! -f "$CTL" ]; then
  echo "Missing control file: $CTL" >&2
  echo "Run this script from simulations/nanodimer/wu-faststart." >&2
  exit 1
fi

(
set -eu

echo "[$(date "+%F %T")] Field snapshot run started"
echo "Using NP=${NP}"
echo "Control file: ${CTL}"
echo "RES=${RES} GAP=${GAP} NBG=${NBG} THETA=${THETA} RUN_TIME=${RUN_TIME} POL=${POL}"
echo "CASES=${CASES}"

for item in ${CASES}; do
  phi="${item%%:*}"
  wvl="${item##*:}"
  phi_tag="$(printf "%03d" "$phi")"
  wvl_tag="${wvl/./}"
  outdir="field_phi${phi_tag}_wvl${wvl_tag}"

  echo "[$(date "+%F %T")] running ${outdir}"
  mkdir -p "$outdir"

  (
    cd "$outdir"
    mpirun -np "$NP" meep-openmpi \
      res="$RES" gap="$GAP" nbg="$NBG" theta="$THETA" phi="$phi" \
      wvl="$wvl" run_time="$RUN_TIME" pol="$POL" \
      "../$CTL" > run.log
  )
done

echo "[$(date "+%F %T")] Field snapshot run completed"
) > "$LOGFILE" 2>&1 &

echo "Started background field snapshot job."
echo "Log: $LOGFILE"
echo "Cases: $CASES"
