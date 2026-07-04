#!/usr/bin/env bash
#
# Asymmetric nanobar metasurface reflectance: single-bar controls + the
# asymmetric pair, both polarizations. Uses asym_reflectance.ctl (corrected
# z-stack + field-decay stopping).
#
#   shape 0 = horizontal bar (control)
#   shape 1 = vertical bar (control)
#   shape 2 = asymmetric H+V pair
#   pol   0 = Ex, 1 = Ey
#
# Reference (mode 0) is a vacuum run that measures the incident power; it is
# saved per polarization and subtracted in each structure run.

set -eu

NP="${NP:-32}"
CTL="${CTL:-asym_reflectance.ctl}"
RES="${RES:-80}"
GAP="${GAP:-0.10}"
NFREQ="${NFREQ:-150}"
WVL_MIN="${WVL_MIN:-0.50}"
WVL_MAX="${WVL_MAX:-1.20}"
DECAY_TIME="${DECAY_TIME:-30}"
DECAY_TOL="${DECAY_TOL:-1e-4}"
SHAPES="${SHAPES:-0 1 2}"
POLS="${POLS:-0 1}"
RESULTS_DIR="${RESULTS_DIR:-results_res${RES}_gap$(printf '%03d' "$(awk "BEGIN{print int(${GAP}*1000+0.5)}")")nm}"
LOGFILE="${LOGFILE:-run_asym_$(date +%Y%m%d_%H%M%S).out}"

mkdir -p "$RESULTS_DIR"

common="res=${RES} gap=${GAP} nfreq=${NFREQ} wvl_min=${WVL_MIN} wvl_max=${WVL_MAX} \
  decay_time=${DECAY_TIME} decay_tol=${DECAY_TOL}"

nohup bash -lc "
set -eu
echo \"[\$(date '+%F %T')] asym metasurface sweep started (NP=${NP})\"
echo \"Results: ${RESULTS_DIR}  RES=${RES} GAP=${GAP} WVL=[${WVL_MIN},${WVL_MAX}]\"

# vacuum incident references (per polarization)
for pol in ${POLS}; do
  mpirun -np ${NP} meep-openmpi mode=0 pol=\$pol ${common} \
    ${CTL} > ${RESULTS_DIR}/reference_incident_pol_\${pol}.txt
done

# structure runs: each shape x polarization
for shape in ${SHAPES}; do
  for pol in ${POLS}; do
    tag=\$(printf 'shape%d_pol%d' \"\$shape\" \"\$pol\")
    echo \"[\$(date '+%F %T')] running \$tag\"
    mpirun -np ${NP} meep-openmpi mode=1 shape=\$shape pol=\$pol ${common} \
      ${CTL} > ${RESULTS_DIR}/\${tag}.txt
  done
done
echo \"[\$(date '+%F %T')] asym metasurface sweep completed\"
" > "$LOGFILE" 2>&1 &

echo "Started background job."
echo "Log: $LOGFILE"
echo "Results: $RESULTS_DIR"
