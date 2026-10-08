import math

from xflow import Case, analytical_benchmark, simulate


def test_synthetic_benchmark_values():
    b = analytical_benchmark(Case())

    assert math.isclose(b["delta_p_hydro_bar"], 20.0, rel_tol=0, abs_tol=1e-12)
    assert math.isclose(b["delta_p0_bar"], 10.0, rel_tol=0, abs_tol=1e-12)

    assert math.isclose(b["j1_m3d_bar"], 379.9622299, rel_tol=1e-8)
    assert math.isclose(b["j2_m3d_bar"], 101.3232613, rel_tol=1e-8)
    assert math.isclose(b["j_eff_m3d_bar"], 79.9920484, rel_tol=1e-8)

    assert math.isclose(b["initial_q_std_m3d"], 799.9204840, rel_tol=1e-8)
    assert math.isclose(b["equilibrium_transfer_std_m3"], 8699.795041, rel_tol=1e-8)
    assert math.isclose(b["p1_equilibrium_bar"], 256.9230769, rel_tol=1e-8)
    assert math.isclose(b["p2_equilibrium_bar"], 236.9230769, rel_tol=1e-8)


def test_wellbore_hydrostatic_consistency():
    result = simulate(Case())
    first = result.data.iloc[0]
    head = (
        first["wellbore_pressure_r1_depth_bar"]
        - first["wellbore_pressure_r2_depth_bar"]
    )
    assert math.isclose(head, 20.0, rel_tol=0, abs_tol=1e-10)


def test_numerical_solution_matches_analytical_tolerance_time():
    case = Case(dt_s=3600.0)
    b = analytical_benchmark(case)
    result = simulate(case)
    t_num = result.data.iloc[-1]["time_days"]

    # Forward-Euler 1 h timestep should be within 0.2 day.
    assert abs(t_num - b["time_to_tolerance_days"]) < 0.2


def test_pressure_difference_converges_to_hydrostatic_head():
    case = Case()
    result = simulate(case)
    last = result.data.iloc[-1]

    pressure_difference = (
        last["pressure_r1_bar"] - last["pressure_r2_bar"]
    )
    assert abs(pressure_difference - 20.0) < 0.02
