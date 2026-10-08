from __future__ import annotations

from dataclasses import dataclass
import math

import pandas as pd


MD_TO_M2 = 9.869233e-16
CP_TO_PA_S = 1.0e-3
BAR_TO_PA = 1.0e5
DAY_TO_S = 86400.0


@dataclass(frozen=True)
class Case:
    """Input data for the synthetic two-reservoir crossflow problem."""

    # Lower / overpressured reservoir R1
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

    # Upper / receiving reservoir R2
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

    # Numerics
    dt_s: float = 3600.0
    max_iterations: int = 10000
    pressure_tolerance_bar: float = 0.01


@dataclass(frozen=True)
class SimulationResult:
    data: pd.DataFrame
    j1_m3d_bar: float
    j2_m3d_bar: float
    j_eff_m3d_bar: float
    pv1_m3: float
    pv2_m3: float
    delta_p_hydro_bar: float


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
    """Return closed-form benchmark values for the linear lumped model."""
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

    transferred_std_m3 = delta_p0 / (
        case.bw1 / (pv1 * ct1)
        + case.bw2 / (pv2 * ct2)
    )

    # Reservoir-volume changes in each tank
    dv1_res = transferred_std_m3 * case.bw1
    dv2_res = transferred_std_m3 * case.bw2

    p1_eq = case.p1_start_bar - dv1_res / (pv1 * ct1)
    p2_eq = case.p2_start_bar + dv2_res / (pv2 * ct2)

    if case.pressure_tolerance_bar > 0:
        time_to_tolerance = tau_days * math.log(
            delta_p0 / case.pressure_tolerance_bar
        )
    else:
        time_to_tolerance = math.inf

    return {
        "pv1_m3": pv1,
        "pv2_m3": pv2,
        "j1_m3d_bar": j1,
        "j2_m3d_bar": j2,
        "j_eff_m3d_bar": j_eff,
        "delta_p_hydro_bar": delta_p_hydro,
        "delta_p0_bar": delta_p0,
        "initial_q_std_m3d": initial_q,
        "tau_days": tau_days,
        "time_to_tolerance_days": time_to_tolerance,
        "equilibrium_transfer_std_m3": transferred_std_m3,
        "p1_equilibrium_bar": p1_eq,
        "p2_equilibrium_bar": p2_eq,
    }


def simulate(case: Case = Case()) -> SimulationResult:
    """Numerically integrate the simplified two-tank crossflow model."""
    bench = analytical_benchmark(case)

    pv1 = bench["pv1_m3"]
    pv2 = bench["pv2_m3"]
    j1 = bench["j1_m3d_bar"]
    j2 = bench["j2_m3d_bar"]
    j_eff = bench["j_eff_m3d_bar"]
    delta_p_hydro = bench["delta_p_hydro_bar"]

    ct1 = case.cw1_bar + case.cf1_bar
    ct2 = case.cw2_bar + case.cf2_bar

    p1 = case.p1_start_bar
    p2 = case.p2_start_bar
    time_s = 0.0
    cumulative_std_m3 = 0.0

    rows: list[dict[str, float]] = []

    for _ in range(case.max_iterations):
        delta_p_drive = p1 - p2 - delta_p_hydro
        if delta_p_drive <= case.pressure_tolerance_bar:
            break

        q_std_m3d = j_eff * delta_p_drive
        q_std_m3s = q_std_m3d / DAY_TO_S

        # Intermediate wellbore pressures for QC.
        pw1_bar = p1 - q_std_m3d / j1
        pw2_bar = p2 + q_std_m3d / j2

        dvol_std = q_std_m3s * case.dt_s
        dv1_res = -dvol_std * case.bw1
        dv2_res = +dvol_std * case.bw2

        dp1 = dv1_res / (pv1 * ct1)
        dp2 = dv2_res / (pv2 * ct2)

        p1 += dp1
        p2 += dp2
        time_s += case.dt_s
        cumulative_std_m3 += dvol_std

        rows.append(
            {
                "time_days": time_s / DAY_TO_S,
                "pressure_r1_bar": p1,
                "pressure_r2_bar": p2,
                "excess_potential_bar": p1 - p2 - delta_p_hydro,
                "wellbore_pressure_r1_depth_bar": pw1_bar,
                "wellbore_pressure_r2_depth_bar": pw2_bar,
                "crossflow_std_m3d": q_std_m3d,
                "cumulative_crossflow_std_m3": cumulative_std_m3,
            }
        )

    if not rows:
        raise RuntimeError("Simulation produced no timesteps.")

    return SimulationResult(
        data=pd.DataFrame(rows),
        j1_m3d_bar=j1,
        j2_m3d_bar=j2,
        j_eff_m3d_bar=j_eff,
        pv1_m3=pv1,
        pv2_m3=pv2,
        delta_p_hydro_bar=delta_p_hydro,
    )
