#!/usr/bin/env bash

set -eu

NP="${NP:-4}"
CTL="${CTL:-nanorod_validation.ctl}"
LOGFILE="${LOGFILE:-run_updated_$(date +%Y%m%d_%H%M%S).out}"

nohup bash -c "
echo \"[$(date '+%F %T')] background run started\"
echo \"Using NP=${NP}\"
echo \"Control file: ${CTL}\"

mpirun -np ${NP} meep-openmpi mode=0 pol=0 geom=0 ${CTL} > inc_ex_single.txt &&
mpirun -np ${NP} meep-openmpi mode=1 pol=0 geom=0 ${CTL} > rod_ex_single.txt &&
mpirun -np ${NP} meep-openmpi mode=0 pol=1 geom=0 ${CTL} > inc_ey_single.txt &&
mpirun -np ${NP} meep-openmpi mode=1 pol=1 geom=0 ${CTL} > rod_ey_single.txt &&
mpirun -np ${NP} meep-openmpi mode=0 pol=0 geom=1 ${CTL} > inc_ex_dimer.txt &&
mpirun -np ${NP} meep-openmpi mode=1 pol=0 geom=1 ${CTL} > rod_ex_dimer.txt &&
mpirun -np ${NP} meep-openmpi mode=0 pol=1 geom=1 ${CTL} > inc_ey_dimer.txt &&
mpirun -np ${NP} meep-openmpi mode=1 pol=1 geom=1 ${CTL} > rod_ey_dimer.txt

status=\$?
if [ \"\$status\" -eq 0 ]; then
  echo \"[$(date '+%F %T')] background run completed\"
else
  echo \"[$(date '+%F %T')] background run failed with status \$status\"
fi
exit \"\$status\"
" > "$LOGFILE" 2>&1 &

echo "Started background job."
echo "Log: $LOGFILE"
