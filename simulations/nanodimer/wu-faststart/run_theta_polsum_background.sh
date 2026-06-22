#!/usr/bin/env bash
#
# Theta sweep with a two-polarization (Ex + Ey) sum at phi = 0.
#
# Summing the two orthogonal polarizations gives an orientation-independent
# (polarization-averaged) spectrum, so every structure angle shows BOTH coupled
# modes. That lets you track each band CENTER vs theta and check Wu's result that
# the resonance wavelengths are independent of the structure angle - instead of
# the phi=0 / Ex-only theta sweep, which mixes geometry change with reorientation
# and makes the dominant peak appear to slide.
#
# Cost: 2 runs per theta + 2 references (cheap vs a full phi sweep per theta).

set -eu

NP="${NP:-32}"
CTL="${CTL:-nanorod_wu_fast.ctl}"
RES="${RES:-200}"
GAP="${GAP:-0.006}"
NBG="${NBG:-1.28}"
NFREQ="${NFREQ:-121}"
WVL_MIN="${WVL_MIN:-0.55}"
WVL_MAX="${WVL_MAX:-0.95}"
DECAY_TIME="${DECAY_TIME:-30}"
DECAY_TOL="${DECAY_TOL:-1e-4}"
THETA_LIST="${THETA_LIST:-30 60 80 90 120 150}"
GAP_TAG="${GAP_TAG:-$(printf "%04d" "$(awk "BEGIN {print int(${GAP} * 1000 + 0.5)}")")}"
RESULTS_DIR="${RESULTS_DIR:-results_sym_gap${GAP_TAG}nm_res${RES}_nfreq${NFREQ}_theta_polsum}"
LOGFILE="${LOGFILE:-run_theta_polsum_$(date +%Y%m%d_%H%M%S).out}"

mkdir -p "$RESULTS_DIR"

nohup bash -lc "
set -eu
echo \"[\$(date '+%F %T')] theta two-polarization sweep started (NP=${NP})\"
echo \"Results dir: ${RESULTS_DIR}\"
echo \"RES=${RES} GAP=${GAP} NBG=${NBG} NFREQ=${NFREQ} WVL=[${WVL_MIN},${WVL_MAX}] THETA_LIST=${THETA_LIST}\"

# one incident reference per polarization (saved flux is pol-specific in the ctl)
for pol in 0 1; do
  mpirun -np ${NP} meep-openmpi \
    mode=0 pol=\$pol \
    res=${RES} gap=${GAP} nbg=${NBG} nfreq=${NFREQ} \
    wvl_min=${WVL_MIN} wvl_max=${WVL_MAX} \
    decay_time=${DECAY_TIME} decay_tol=${DECAY_TOL} \
    ${CTL} > ${RESULTS_DIR}/reference_incident_pol_\${pol}.txt
done

for theta in ${THETA_LIST}; do
  for pol in 0 1; do
    tag=\$(printf 'theta_%03d_pol_%d' \"\$theta\" \"\$pol\")
    echo \"[\$(date '+%F %T')] running \$tag\"
    mpirun -np ${NP} meep-openmpi \
      mode=1 theta=\$theta phi=0 pol=\$pol \
      res=${RES} gap=${GAP} nbg=${NBG} nfreq=${NFREQ} \
      wvl_min=${WVL_MIN} wvl_max=${WVL_MAX} \
      decay_time=${DECAY_TIME} decay_tol=${DECAY_TOL} \
      ${CTL} > ${RESULTS_DIR}/\${tag}.txt
  done
done

echo \"[\$(date '+%F %T')] theta two-polarization sweep completed\"
" > "$LOGFILE" 2>&1 &

echo "Started background job."
echo "Log: $LOGFILE"
echo "Results: $RESULTS_DIR"
