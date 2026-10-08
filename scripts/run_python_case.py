from pathlib import Path

import matplotlib.pyplot as plt

from xflow import Case, analytical_benchmark, simulate


def main() -> None:
    case = Case()
    benchmark = analytical_benchmark(case)
    result = simulate(case)
    df = result.data

    print("=== Analytical benchmark ===")
    for key, value in benchmark.items():
        print(f"{key}: {value:.6f}")

    print("\n=== Numerical result ===")
    last = df.iloc[-1]
    print(f"simulation_time_days: {last['time_days']:.6f}")
    print(f"final_r1_pressure_bar: {last['pressure_r1_bar']:.6f}")
    print(f"final_r2_pressure_bar: {last['pressure_r2_bar']:.6f}")
    print(f"final_excess_potential_bar: {last['excess_potential_bar']:.6f}")
    print(
        "cumulative_crossflow_std_m3: "
        f"{last['cumulative_crossflow_std_m3']:.6f}"
    )

    outdir = Path("outputs")
    outdir.mkdir(exist_ok=True)
    df.to_csv(outdir / "simulation_results.csv", index=False)

    fig, ax1 = plt.subplots(figsize=(10, 6))
    ax1.set_xlabel("Time (days)")
    ax1.set_ylabel("Crossflow rate (std m3/day)")
    ax1.plot(df["time_days"], df["crossflow_std_m3d"], label="Crossflow rate")
    ax1.grid(True)

    ax2 = ax1.twinx()
    ax2.set_ylabel("Pressure (bar)")
    ax2.plot(
        df["time_days"],
        df["pressure_r1_bar"],
        linestyle="--",
        label="Reservoir 1",
    )
    ax2.plot(
        df["time_days"],
        df["pressure_r2_bar"],
        linestyle=":",
        label="Reservoir 2",
    )

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="best")
    plt.title("Wellbore crossflow and reservoir pressures")
    plt.tight_layout()

    plt.figure(figsize=(10, 6))
    plt.plot(
        df["time_days"],
        df["cumulative_crossflow_std_m3"],
        label="Cumulative crossflow",
    )
    plt.xlabel("Time (days)")
    plt.ylabel("Cumulative crossflow (std m3)")
    plt.title("Cumulative wellbore crossflow")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    plt.show()


if __name__ == "__main__":
    main()
