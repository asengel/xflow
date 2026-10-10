from __future__ import annotations

from dataclasses import dataclass
from collections.abc import Iterable
import math

import pandas as pd


MD_TO_M2 = 9.869233e-16
CP_TO_PA_S = 1.0e-3
BAR_TO_PA = 1.0e5
DAY_TO_S = 86400.0


@dataclass(frozen=True)
class Case:
    """Input data for the synthetic two-reservoir crossflow problem."""

    # Reservoir R1 (arbitrary depth, pressure and flow direction)
    phi1: float = 0.15
    skin1: float = 0.0
    cw1_bar: float = 4.0e-5
    cf1_bar: float = 6.0e-5
    mu1_cp: float = 0.50
    p1_start_bar: float = 260.0
    depth1_m: float = 2500.0
    h1_m: float = 60.0
    k1_md: float = 500.0
    bw1: float = 1.0

    # Reservoir R2 (arbitrary depth, pressure and flow direction)
    phi2: float = 0.10
    skin2: float = 0.0
    cw2_bar: float = 4.0e-5
    cf2_bar: float = 6.0e-5
    mu2_cp: float = 0.50
    p2_start_bar: float = 230.0
    depth2_m: float = 2300.0
    h2_m: float = 40.0
    k2_md: float = 200.0
    bw2: float = 1.0

    # Shared geometry / hydraulics
    pressure_gradient_bar_m: float = 0.10
    re_m: float = 1000.0
    rw_m: float = 0.10

    # Analytical stopping threshold (not a numerical integration timestep)
    pressure_tolerance_bar: float = 0.01


def _productivity_index_surface(
    *,
    k_md: float,
    h_m: float,
    mu_cp: float,
    bw: float,
    skin: float,
    re_m: float,
    rw_m: float,
) -> float:
    """Pseudo-steady-state PI in standard m3/day/bar."""
    k_m2 = k_md * MD_TO_M2
    mu_pa_s = mu_cp * CP_TO_PA_S
    denominator = mu_pa_s * bw * (math.log(re_m / rw_m) - 0.75 + skin)
    j_m3s_pa = 2.0 * math.pi * k_m2 * h_m / denominator
    return j_m3s_pa * BAR_TO_PA * DAY_TO_S


def analytical_benchmark(case: Case = Case()) -> dict[str, float]:
    """Return exact benchmark quantities for the linear lumped model."""
    pv1 = math.pi * case.re_m**2 * case.h1_m * case.phi1
    pv2 = math.pi * case.re_m**2 * case.h2_m * case.phi2

    ct1 = case.cw1_bar + case.cf1_bar
    ct2 = case.cw2_bar + case.cf2_bar

    j1 = _productivity_index_surface(
        k_md=case.k1_md,
        h_m=case.h1_m,
        mu_cp=case.mu1_cp,
        bw=case.bw1,
        skin=case.skin1,
        re_m=case.re_m,
        rw_m=case.rw_m,
    )
    j2 = _productivity_index_surface(
        k_md=case.k2_md,
        h_m=case.h2_m,
        mu_cp=case.mu2_cp,
        bw=case.bw2,
        skin=case.skin2,
        re_m=case.re_m,
        rw_m=case.rw_m,
    )
    if min(pv1, pv2, ct1, ct2, j1, j2, case.bw1, case.bw2) <= 0:
        raise ValueError("Pore volumes, compressibilities, volume factors and PIs must be positive")
    if case.pressure_tolerance_bar < 0:
        raise ValueError("pressure_tolerance_bar cannot be negative")

    j_eff = 1.0 / (1.0 / j1 + 1.0 / j2)

    delta_p_hydro = (
        case.depth1_m - case.depth2_m
    ) * case.pressure_gradient_bar_m

    delta_p0 = case.p1_start_bar - case.p2_start_bar - delta_p_hydro
    initial_q = j_eff * delta_p0

    storage_inv = (
        case.bw1 / (pv1 * ct1)
        + case.bw2 / (pv2 * ct2)
    )
    lambda_per_day = j_eff * storage_inv
    tau_days = 1.0 / lambda_per_day

    transferred_std_m3 = delta_p0 / storage_inv

    dv1_res = transferred_std_m3 * case.bw1
    dv2_res = transferred_std_m3 * case.bw2

    p1_eq = case.p1_start_bar - dv1_res / (pv1 * ct1)
    p2_eq = case.p2_start_bar + dv2_res / (pv2 * ct2)

    if abs(delta_p0) <= case.pressure_tolerance_bar:
        time_to_tolerance = 0.0
    elif case.pressure_tolerance_bar == 0:
        time_to_tolerance = math.inf
    else:
        time_to_tolerance = tau_days * math.log(
            abs(delta_p0) / case.pressure_tolerance_bar
        )

    return {
        "pv1_m3": pv1,
        "pv2_m3": pv2,
        "j1_m3d_bar": j1,
        "j2_m3d_bar": j2,
        "j_eff_m3d_bar": j_eff,
        "delta_p_hydro_bar": delta_p_hydro,
        "delta_p0_bar": delta_p0,
        "initial_q_std_m3d": initial_q,
        "flow_direction": ("R1_to_R2" if delta_p0 > 0 else "R2_to_R1" if delta_p0 < 0 else "equilibrium"),
        "tau_days": tau_days,
        "time_to_tolerance_days": time_to_tolerance,
        "equilibrium_transfer_std_m3": transferred_std_m3,
        "p1_equilibrium_bar": p1_eq,
        "p2_equilibrium_bar": p2_eq,
    }


def analytical_timeseries(
    case: Case = Case(), *, times_days: Iterable[float]
) -> pd.DataFrame:
    """Return exact crossflow states at requested nonnegative times (days).

    Preserves the input time order, including duplicates; it does not integrate
    over timesteps. Signed rates and transfers are positive from R1 to R2.
    """
    b = analytical_benchmark(case)
    a1 = case.bw1 / (b["pv1_m3"] * (case.cw1_bar + case.cf1_bar))
    a2 = case.bw2 / (b["pv2_m3"] * (case.cw2_bar + case.cf2_bar))
    rows = []

    for t in times_days:
        t = float(t)
        if not math.isfinite(t) or t < 0:
            raise ValueError("All times_days must be finite and nonnegative")

        exponent = -t / b["tau_days"]
        decay = math.exp(exponent)
        fraction_transferred = -math.expm1(exponent)
        transferred = b["equilibrium_transfer_std_m3"] * fraction_transferred
        rate = b["initial_q_std_m3d"] * decay
        p1 = case.p1_start_bar - a1 * transferred
        p2 = case.p2_start_bar + a2 * transferred

        rows.append({
            "time_days": t,
            "pressure_r1_bar": p1,
            "pressure_r2_bar": p2,
            "excess_potential_bar": b["delta_p0_bar"] * decay,
            "wellbore_pressure_r1_depth_bar": p1 - rate / b["j1_m3d_bar"],
            "wellbore_pressure_r2_depth_bar": p2 + rate / b["j2_m3d_bar"],
            "crossflow_std_m3d": rate,
            "cumulative_crossflow_std_m3": transferred,
        })

    if not rows:
        raise ValueError("times_days must contain at least one time")
    return pd.DataFrame.from_records(rows)
