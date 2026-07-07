#!/usr/bin/env bash
#
# Asymmetric nanobar metasurface sweep. Each invocation of asym_reflectance.ctl
# is self-contained: it runs the bare-substrate reference pass and the full
# metasurface pass internally and prints the differential reflected + transmitted
# flux. So we just loop over shape x polarization — no separate reference run.
#
#   shape 0 = horizontal bar (control), 1 = vertical bar (control), 2 = asym pair
#   pol   0 = Ex (via 90-deg geometry rotation), 1 = Ey
#
set -eu

NP="${NP:-32}"
CTL="${CTL:-asym_reflectance.ctl}"
RES="${RES:-80}"
GAP="${GAP:-0.10}"
NFREQ="${NFREQ:-150}"
WVL_MIN="${WVL_MIN:-0.50}"
WVL_MAX="${WVL_MAX:-1.20}"
SHAPES="${SHAPES:-0 1 2}"
POLS="${POLS:-0 1}"
RESULTS_DIR="${RESULTS_DIR:-results_res${RES}_gap$(printf '%03d' "$(awk "BEGIN{print int(${GAP}*1000+0.5)}")")nm}"
LOGFILE="${LOGFILE:-run_asym_$(date +%Y%m%d_%H%M%S).out}"

mkdir -p "$RESULTS_DIR"

nohup bash -lc "
set -eu
echo \"[\$(date '+%F %T')] asym metasurface sweep started (NP=${NP})\"
echo \"Results: ${RESULTS_DIR}  RES=${RES} GAP=${GAP} WVL=[${WVL_MIN},${WVL_MAX}]\"
for shape in ${SHAPES}; do
  for pol in ${POLS}; do
    tag=\$(printf 'shape%d_pol%d' \"\$shape\" \"\$pol\")
    echo \"[\$(date '+%F %T')] running \$tag\"
    mpirun -np ${NP} meep-openmpi shape=\$shape pol=\$pol \
      res=${RES} gap=${GAP} nfreq=${NFREQ} wvl_min=${WVL_MIN} wvl_max=${WVL_MAX} \
      ${CTL} > ${RESULTS_DIR}/\${tag}.txt
  done
done
echo \"[\$(date '+%F %T')] asym metasurface sweep completed\"
" > "$LOGFILE" 2>&1 &

echo "Started background job."
echo "Log: $LOGFILE"
echo "Results: $RESULTS_DIR"
