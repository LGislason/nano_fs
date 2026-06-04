#!/usr/bin/env bash

set -eu

NP="${NP:-4}"
CTL="${CTL:-nanorod_wu_fast.ctl}"

# Normal-incidence fixed-theta phi sweep defaults. Override these on the command line for
# cheaper test runs or final refinement. GAP is the true metal-to-metal surface
# gap; the control file derives the internal tip-reference spacing.
RES="${RES:-200}"
GAP="${GAP:-0.006}"
NBG="${NBG:-1.28}"
NFREQ="${NFREQ:-121}"
WVL_MIN="${WVL_MIN:-0.55}"
WVL_MAX="${WVL_MAX:-0.95}"
DECAY_TIME="${DECAY_TIME:-30}"
DECAY_TOL="${DECAY_TOL:-1e-4}"
THETA_LIST="${THETA_LIST:-80}"
PHI_LIST="${PHI_LIST:-0 20 40 60 80 90 100 120 140 160 180}"
GAP_TAG="${GAP_TAG:-$(printf "%04d" "$(awk "BEGIN {print int(${GAP} * 1000 + 0.5)}")")}"
RESULTS_DIR="${RESULTS_DIR:-results_sym_gap${GAP_TAG}nm_res${RES}_nfreq${NFREQ}_theta080_phi_sweep}"
LOGFILE="${LOGFILE:-run_wu_fast_$(date +%Y%m%d_%H%M%S).out}"

mkdir -p "$RESULTS_DIR"

nohup bash -lc "
set -eu

echo \"[\$(date '+%F %T')] Wu fixed-theta phi-sweep run started\"
echo \"Using NP=${NP}\"
echo \"Control file: ${CTL}\"
echo \"Results dir: ${RESULTS_DIR}\"
echo \"RES=${RES} GAP=${GAP} NBG=${NBG} NFREQ=${NFREQ} WVL_MIN=${WVL_MIN} WVL_MAX=${WVL_MAX}\"
echo \"DECAY_TIME=${DECAY_TIME} DECAY_TOL=${DECAY_TOL}\"
echo \"THETA_LIST=${THETA_LIST}\"
echo \"PHI_LIST=${PHI_LIST}\"

mpirun -np ${NP} meep-openmpi \
  mode=0 \
  res=${RES} gap=${GAP} nbg=${NBG} nfreq=${NFREQ} \
  wvl_min=${WVL_MIN} wvl_max=${WVL_MAX} \
  decay_time=${DECAY_TIME} decay_tol=${DECAY_TOL} \
  ${CTL} > ${RESULTS_DIR}/reference_incident.txt

for theta in ${THETA_LIST}; do
  for phi in ${PHI_LIST}; do
    tag=\$(printf 'theta_%03d_phi_%03d' \"\$theta\" \"\$phi\")
    echo \"[\$(date '+%F %T')] running \$tag\"
    mpirun -np ${NP} meep-openmpi \
      mode=1 theta=\$theta phi=\$phi \
      res=${RES} gap=${GAP} nbg=${NBG} nfreq=${NFREQ} \
      wvl_min=${WVL_MIN} wvl_max=${WVL_MAX} \
      decay_time=${DECAY_TIME} decay_tol=${DECAY_TOL} \
      ${CTL} > ${RESULTS_DIR}/\${tag}.txt
  done
done

status=\$?
if [ \"\$status\" -eq 0 ]; then
  echo \"[\$(date '+%F %T')] Wu fixed-theta phi-sweep run completed\"
else
  echo \"[\$(date '+%F %T')] Wu fixed-theta phi-sweep run failed with status \$status\"
fi
exit \"\$status\"
" > "$LOGFILE" 2>&1 &

echo "Started background job."
echo "Log: $LOGFILE"
echo "Results: $RESULTS_DIR"
