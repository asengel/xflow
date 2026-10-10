# XFlow validation: analytical model versus OPM Flow

This benchmark compares the **closed-form two-reservoir analytical solution** with a separate OPM Flow simulation for the default `Case()` parameters. The OPM Flow input deck is `simulation_deck/CASE.DATA` in the parent repository.

## Time-matched comparison (day 100)

| Quantity | Analytical | OPM Flow | Analytical minus OPM |
|---|---:|---:|---:|
| Crossflow (Sm³/day) | 0.081253 | 0.085080 | −0.003827 |
| Pressure R1 (bar) | 256.923389 | 256.921570 | +0.001820 |
| Pressure R2 (bar) | 236.922374 | 236.922974 | −0.000600 |
| Cumulative transfer (Sm³) | 8,698.911352 | 8,702.685547 | −3.774195 |

Cumulative transfer differs by **−0.04337%** relative to OPM at day 100. This compares a finite-time OPM output to the analytical result at the same time; it does **not** confuse that with the analytical asymptotic equilibrium volume (8,699.795041 Sm³).

The analytical day-zero crossflow rate is 799.920484 Sm³/day; at day 0.1 it is 792.599161 Sm³/day, versus the first OPM rate of 792.633850 Sm³/day. The numerical time series in `opm_flow_timeseries.csv` is extracted simulator output.

To reproduce the time-matched comparison tables and plots from the bundled reference data, run `scripts/run_python_case.py --no-show`. The OPM reference is fixed to the bundled validation case; users do not need to process simulator summary files.

