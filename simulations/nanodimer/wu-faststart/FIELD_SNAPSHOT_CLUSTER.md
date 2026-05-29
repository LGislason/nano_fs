# Cluster Field Snapshots

Use this when the cluster has old Scheme Meep but no Python bindings.

## Run on the cluster

From `simulations/nanodimer/wu-faststart`, launch the snapshot wrapper:

```bash
chmod +x run_field_snapshots_background.sh
./run_field_snapshots_background.sh
```

The script uses the same cluster launch style as `run_wu_fast_background.sh`:
`mpirun -np ${NP} meep-openmpi ...`. The raw command format is also shown in
the header of `nanorod_wu_field_snapshot.ctl`.

Recommended first cases:

- `field_phi040_wvl0672`: main 0.672 um mode near its angular maximum.
- `field_phi120_wvl0950`: weak 0.950 um candidate near its angular maximum.

If runtime is acceptable, also run the crossed cases:

- `field_phi040_wvl0950`
- `field_phi120_wvl0672`

Use `res=300`, `gap=0.006`, `theta=80`, and `run_time=300` unless you are doing
a quick smoke test.

For a smoke test, run:

```bash
RES=30 RUN_TIME=5 CASES="40:0.672" ./run_field_snapshots_background.sh
```

Expected HDF5 files include:

```text
xy_rod_ex*.h5
xy_rod_ey*.h5
xy_rod_ez*.h5
xy_rod_eps*.h5
xy_above_ex*.h5
xy_above_ey*.h5
xy_above_ez*.h5
xy_above_eps*.h5
xz_ex*.h5
xz_ey*.h5
xz_ez*.h5
xz_eps*.h5
```

The outputs are 2D slices, not the full 3D cell.

## Plot locally

Copy each output directory back to this repository, then run:

```bash
conda run -n meep-env python plot_field_snapshot_h5.py \
  field_phi040_wvl0672 \
  -o field_phi040_wvl0672.png

conda run -n meep-env python plot_field_snapshot_h5.py \
  field_phi120_wvl0950 \
  -o field_phi120_wvl0950.png
```

Use `--contrast total` if you want raw `|E|^2` instead of edge-normalized
relative contrast.

## Notes

These snapshots use a continuous source and output time-domain field slices at
the end of the run. They are meant for high-resolution visual inspection. If
your cluster Meep supports DFT field output from Scheme, that would be better
for phase-independent steady-state intensity, but this file sticks to older
Meep output functions for compatibility.
