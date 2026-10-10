"""Analytical wellbore crossflow driven by a JSON case file.

Run from the repository root:
  python scripts/run_python_case.py
  python scripts/run_python_case.py cases/my_case.json
  python scripts/run_python_case.py cases/validation_case.json --no-show

Actual OPM Flow series are plotted only when running the bundled
cases/validation_case.json file explicitly. No interactive prompts or numerical integration.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, fields
import json
import math
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from xflow import Case, analytical_benchmark, analytical_timeseries

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OPM_CSV = ROOT / "validation" / "opm_flow_timeseries.csv"
DEFAULT_CASE_JSON = ROOT / "cases" / "validation_case.json"

COMPARISON_COLUMNS = (
    ("crossflow_std_m3d", "Crossflow rate", "Sm3/day"),
    ("pressure_r1_bar", "Reservoir 1 pressure", "bar"),
    ("pressure_r2_bar", "Reservoir 2 pressure", "bar"),
    ("cumulative_crossflow_std_m3", "Cumulative transfer", "Sm3"),
)


def case_from_json(path: Path) -> Case:
    """Load explicit reservoir inputs; require every parameter for reproducibility."""
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError("Input JSON must be an object of Case parameter names and numbers")
    expected = {f.name for f in fields(Case)}
    unexpected = sorted(set(data) - expected)
    missing = sorted(expected - set(data))
    if unexpected:
        raise ValueError(f"Unknown case inputs: {', '.join(unexpected)}")
    if missing:
        raise ValueError(f"Missing case inputs: {', '.join(missing)}. Copy cases/validation_case.json as a template")
    values = {}
    for key, value in data.items():
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
            raise ValueError(f"Case input '{key}' must be a finite number")
        values[key] = float(value)
    return Case(**values)


def time_samples(benchmark: dict[str, float], days: float | None,
                 n_points: int, opm: pd.DataFrame | None) -> np.ndarray:
    if n_points < 2:
        raise ValueError("--points must be at least 2")
    if days is not None:
        end_days = days
    elif opm is not None:
        # Default validation run MUST span the full supplied OPM benchmark.
        end_days = float(opm.time_days.iloc[-1])
    else:
        end_days = benchmark["time_to_tolerance_days"]
    if not math.isfinite(end_days):
        raise ValueError("Zero pressure tolerance has no finite cutoff; specify --days")
    if end_days < 0:
        raise ValueError("--days cannot be negative")
    return np.linspace(0, end_days, n_points) if end_days else np.array([0.0])


def load_opm_timeseries(path: Path) -> pd.DataFrame:
    """Read measured/extracted simulator values, never synthesize OPM curves."""
    df = pd.read_csv(path)
    ycols = [name for name, _title, _unit in COMPARISON_COLUMNS]
    present = [name for name in ycols if name in df.columns]
    if "time_days" not in df or not present:
        raise ValueError("OPM CSV needs time_days and at least one comparison column")
    cols = ["time_days", *present]
    for col in cols:
        df[col] = pd.to_numeric(df[col], errors="raise")
    if df[cols].isna().any().any() or not np.isfinite(df[cols].to_numpy(dtype=float)).all():
        raise ValueError("OPM CSV contains missing or nonfinite values")
    if (df.time_days < 0).any():
        raise ValueError("OPM CSV contains negative times")
    df = df[cols].sort_values("time_days").reset_index(drop=True)
    if df.time_days.duplicated().any():
        raise ValueError("OPM CSV contains duplicate timestamps")
    return df


def compare_at_opm_times(case: Case, opm: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return a time-aligned raw comparison and metric statistics."""
    analytical = analytical_timeseries(case, times_days=opm.time_days)
    aligned = pd.DataFrame({"time_days": opm.time_days.to_numpy()})
    rows = []
    for column, title, units in COMPARISON_COLUMNS:
        if column not in opm:
            continue
        actual = opm[column].to_numpy(dtype=float)
        predicted = analytical[column].to_numpy(dtype=float)
        error = predicted - actual
        aligned[f"analytical_{column}"] = predicted
        aligned[f"opm_{column}"] = actual
        aligned[f"error_{column}"] = error
        rows.append({
            "variable": title,
            "unit": units,
            "first_time_days": float(opm.time_days.iloc[0]),
            "analytical_first": float(predicted[0]),
            "opm_first": float(actual[0]),
            "last_time_days": float(opm.time_days.iloc[-1]),
            "analytical_last": float(predicted[-1]),
            "opm_last": float(actual[-1]),
            "last_difference_analytical_minus_opm": float(error[-1]),
            "mae": float(np.mean(np.abs(error))),
            "rmse": float(np.sqrt(np.mean(error**2))),
            "max_absolute_error": float(np.max(np.abs(error))),
        })
    return aligned, pd.DataFrame(rows)


def plot_results(df: pd.DataFrame, case_label: str, outdir: Path,
                 opm: pd.DataFrame | None,
                 aligned: pd.DataFrame | None = None) -> list[Path]:
    """Compare actual OPM observations (visible symbols) and analytical curves.

    Separate residual axes expose millibar-scale pressure differences that
    would be invisible on the absolute-pressure scale.
    """
    t = df.time_days

    def opm_markers(ax, column: str, label: str, marker: str = "o") -> None:
        # Display representative *actual* samples, never invented or fitted OPM data.
        assert opm is not None
        locs = np.unique(np.linspace(0, len(opm) - 1,
                                     min(len(opm), 22), dtype=int))
        subset = opm.iloc[locs]
        ax.scatter(subset.time_days, subset[column], s=44, marker=marker,
                   facecolors="white", edgecolors="black", linewidths=1.25,
                   label=label, zorder=6)

    fig_rate, ax_q = plt.subplots(figsize=(9.5, 5.7))
    ax_q.plot(t, df.crossflow_std_m3d, label="XFlow analytical", linewidth=2.1,
              color="tab:blue", zorder=2)
    if opm is not None and "crossflow_std_m3d" in opm:
        opm_markers(ax_q, "crossflow_std_m3d", "OPM Flow (recorded samples)")
    ax_q.set(xlabel="Time (days)", ylabel="Signed crossflow R1 → R2 (Sm³/day)",
             title=f"Wellbore crossflow rate | {case_label}")
    ax_q.grid(alpha=0.3)
    ax_q.legend()
    fig_rate.tight_layout()
    path_rate = outdir / "crossflow_rate_vs_time.png"
    fig_rate.savefig(path_rate, dpi=180, bbox_inches="tight")

    have_residuals = (aligned is not None and
                      {"error_pressure_r1_bar", "error_pressure_r2_bar"}
                      <= set(aligned.columns))
    if have_residuals:
        fig_p, (ax_p, ax_res) = plt.subplots(
            2, 1, figsize=(10, 7.6), sharex=True,
            gridspec_kw={"height_ratios": [2.4, 1.0]},
        )
    else:
        fig_p, ax_p = plt.subplots(figsize=(9.5, 5.7))
    ax_p.plot(t, df.pressure_r1_bar, label="XFlow R1",
              linewidth=2.1, color="tab:blue", zorder=2)
    ax_p.plot(t, df.pressure_r2_bar, label="XFlow R2",
              linewidth=2.1, color="tab:orange", zorder=2)
    if opm is not None:
        if "pressure_r1_bar" in opm:
            opm_markers(ax_p, "pressure_r1_bar", "OPM R1 (block 3)", "o")
        if "pressure_r2_bar" in opm:
            opm_markers(ax_p, "pressure_r2_bar", "OPM R2 (block 1)", "s")
    ax_p.set(ylabel="Reservoir pressure (bar)",
             title=f"Reservoir pressures | {case_label}")
    if have_residuals:
        ax_p.set_title(f"Reservoir pressures | {case_label}\n"
                       "Solid lines: XFlow · Hollow symbols: OPM Flow")
    ax_p.grid(alpha=0.3)
    ax_p.legend(ncol=2, fontsize=9, loc="center right")
    if have_residuals:
        ax_res.plot(aligned.time_days, aligned.error_pressure_r1_bar * 1000,
                    color="tab:blue", label="R1: XFlow − OPM", linewidth=1.9)
        ax_res.plot(aligned.time_days, aligned.error_pressure_r2_bar * 1000,
                    color="tab:orange", label="R2: XFlow − OPM", linewidth=1.9)
        ax_res.axhline(0, color="black", linewidth=0.8, linestyle=":")
        ax_res.set(xlabel="Time (days)", ylabel="Difference (mbar)")
        ax_res.grid(alpha=0.3)
        ax_res.legend(fontsize=8, loc="upper right")
        ax_res.set_xlim(float(t.iloc[0]), float(t.iloc[-1]))
    else:
        ax_p.set_xlabel("Time (days)")
    fig_p.tight_layout()
    path_p = outdir / "reservoir_pressures_vs_time.png"
    fig_p.savefig(path_p, dpi=180, bbox_inches="tight")

    fig_vol, ax_v = plt.subplots(figsize=(9.5, 5.7))
    ax_v.plot(t, df.cumulative_crossflow_std_m3, label="XFlow analytical",
              linewidth=2.1, color="tab:blue", zorder=2)
    if opm is not None and "cumulative_crossflow_std_m3" in opm:
        opm_markers(ax_v, "cumulative_crossflow_std_m3", "OPM Flow (recorded samples)")
    ax_v.set(xlabel="Time (days)", ylabel="Signed cumulative transfer R1 → R2 (Sm³)",
             title=f"Cumulative wellbore crossflow | {case_label}")
    ax_v.grid(alpha=0.3)
    ax_v.legend()
    fig_vol.tight_layout()
    path_vol = outdir / "cumulative_crossflow_vs_time.png"
    fig_vol.savefig(path_vol, dpi=180, bbox_inches="tight")
    return [path_rate, path_p, path_vol]


def run(case_file: Path, case: Case, output_dir: Path, days: float | None,
        n_points: int, show: bool) -> None:
    benchmark = analytical_benchmark(case)
    is_validation = (case_file.resolve() == DEFAULT_CASE_JSON.resolve() and case == Case())
    case_label = "Validation case" if is_validation else case_file.stem.replace("_", " ").title()

    opm = None
    if is_validation:
        source = DEFAULT_OPM_CSV
        if not source.is_file():
            raise FileNotFoundError(f"OPM benchmark CSV not found: {source}")
        opm = load_opm_timeseries(source)
        print(f"Using actual OPM Flow benchmark: {source}")
        print(f"OPM samples: {len(opm)}; {opm.time_days.iloc[0]:.2f} to {opm.time_days.iloc[-1]:.2f} days")
    else:
        print("Custom case: analytical curves only (OPM comparison is reserved for the bundled validation case)")

    times = time_samples(benchmark, days, n_points, opm)
    df = analytical_timeseries(case, times_days=times)
    output_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_dir / "analytical_timeseries.csv", index=False)
    (output_dir / "inputs.json").write_text(json.dumps(asdict(case), indent=2) + "\n", encoding="utf-8")
    (output_dir / "analytical_summary.json").write_text(json.dumps(benchmark, indent=2) + "\n", encoding="utf-8")

    print("\n=== XFlow analytical results ===")
    print(f"Inputs: {case_file.resolve()}")
    print(f"Flow direction: {benchmark['flow_direction'].replace('_', ' ')}")
    print(f"Initial crossflow (t=0; signed R1 to R2): {benchmark['initial_q_std_m3d']:,.4f} Sm3/day")
    print(f"Time constant: {benchmark['tau_days']:.4f} days")
    print(f"Equilibrium transfer (signed R1 to R2; t=infinity): {benchmark['equilibrium_transfer_std_m3']:,.4f} Sm3")
    print(f"Equilibrium R1 pressure: {benchmark['p1_equilibrium_bar']:.4f} bar")
    print(f"Equilibrium R2 pressure: {benchmark['p2_equilibrium_bar']:.4f} bar")
    print(f"Last plotted time: {df.time_days.iloc[-1]:.4f} days")

    aligned = None
    if opm is not None:
        opm = opm[opm.time_days <= df.time_days.iloc[-1] + 1e-9].copy()
        if opm.empty:
            print("No OPM timestamps within selected reporting horizon; skipping comparison")
            opm = None
        else:
            aligned, stats = compare_at_opm_times(case, opm)
            aligned.to_csv(output_dir / "opm_analytical_time_aligned.csv", index=False)
            stats.to_csv(output_dir / "opm_comparison_metrics.csv", index=False)
            print("\n=== Time-matched analytical vs actual OPM Flow ===")
            print(f"First OPM sample: day {opm.time_days.iloc[0]:.2f} (not t=0)")
            print(f"Last OPM sample:  day {opm.time_days.iloc[-1]:.2f}")
            display = stats[["variable", "analytical_last", "opm_last",
                             "last_difference_analytical_minus_opm", "rmse"]]
            print(display.to_string(index=False, float_format=lambda v: f"{v:,.5f}"))
            print("Difference = analytical - OPM; RMSE uses all matched samples")

    paths = plot_results(df, case_label, output_dir, opm, aligned)
    print(f"\nAnalytical CSV: {(output_dir / 'analytical_timeseries.csv').resolve()}")
    for path in paths:
        print(f"Figure: {path.resolve()}")
    if show:
        plt.show()
    plt.close("all")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Run analytical crossflow from a JSON case file (defaults to cases/validation_case.json)"
    )
    parser.add_argument("case_file", nargs="?", type=Path, default=DEFAULT_CASE_JSON,
                        help="Full case JSON; default is the bundled validation case")
    parser.add_argument("--days", type=float,
                        help="Reporting horizon; validation defaults to OPM last time, other cases to tolerance time")
    parser.add_argument("--points", type=int, default=500,
                        help="Analytical plotting samples (default 500)")
    parser.add_argument("--output-dir", type=Path,
                        help="Override the output directory")
    parser.add_argument("--no-show", action="store_true",
                        help="Save results and figures without opening plot windows")
    args = parser.parse_args(argv)
    try:
        case = case_from_json(args.case_file)
        name = "validation" if (args.case_file.resolve() == DEFAULT_CASE_JSON.resolve() and case == Case()) else args.case_file.stem
        output_dir = args.output_dir or ROOT / "outputs" / name
        run(args.case_file, case, output_dir, args.days, args.points,
            not args.no_show)
        return 0
    except (ValueError, OSError, json.JSONDecodeError) as exc:
        print(f"XFlow error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
