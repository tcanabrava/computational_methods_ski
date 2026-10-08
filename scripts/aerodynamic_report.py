"""Check the Seo coefficient tables and plot Figures 4–6.

Run: python -m scripts.aerodynamic_report
Requires matplotlib.
"""

import argparse
import csv
from decimal import Decimal
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scripts.tomaz_aerodynamics import aerodynamic_properties, coefficient_rows

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results" / "aerodynamics"


def audit_coefficients():
    """Check the CSV against the printed values and precision in Tables 1–3."""
    # Seo et al. (2004), p. 36: i, j, k, a_ijk, b_ijk, c_ijk.
    # This reference copy stays separate from the CSV used by the model.
    table = [
        (0, 0, 0, "8.54e-2", "4.79e-2", "7.87e-3"),
        (0, 0, 1, "5.33e-4", "-5.78e-3", "-2.17e-3"),
        (0, 0, 2, "-1.24e-4", "2.63e-4", "1.64e-4"),
        (0, 1, 0, "7.90e-5", "5.91e-3", "2.52e-3"),
        (0, 1, 1, "7.05e-5", "2.57e-4", "4.41e-4"),
        (0, 1, 2, "9.07e-6", "-9.12e-6", "-1.58e-5"),
        (0, 2, 0, "1.50e-4", "2.79e-5", "1.36e-5"),
        (0, 2, 1, "-2.37e-6", "-3.78e-6", "-9.50e-6"),
        (0, 2, 2, "-1.67e-7", "-1.62e-8", "3.21e-7"),
        (1, 0, 0, "3.74e-4", "4.78e-3", "5.77e-5"),
        (1, 0, 1, "-3.82e-4", "-1.02e-3", "5.39e-4"),
        (1, 0, 2, "2.38e-5", "7.87e-5", "-2.63e-5"),
        (1, 1, 0, "2.71e-4", "6.50e-4", "6.90e-4"),
        (1, 1, 1, "-3.63e-5", "-1.41e-5", "-1.99e-4"),
        (1, 1, 2, "-6.02e-8", "-3.41e-6", "7.54e-6"),
        (1, 2, 0, "3.74e-6", "-1.51e-5", "-1.41e-5"),
        (1, 2, 1, "9.61e-7", "1.15e-6", "4.74e-6"),
        (1, 2, 2, "-9.33e-9", "5.33e-8", "-1.88e-7"),
        (2, 0, 0, "4.53e-5", "5.20e-4", "5.16e-5"),
        (2, 0, 1, "4.19e-5", "9.20e-5", "-3.86e-5"),
        (2, 0, 2, "-1.38e-6", "-5.02e-6", "1.41e-6"),
        (2, 1, 0, "4.25e-5", "4.39e-6", "-3.19e-5"),
        (2, 1, 1, "6.17e-7", "1.59e-6", "1.35e-5"),
        (2, 1, 2, "3.64e-8", "2.02e-7", "-5.61e-7"),
        (2, 2, 0, "-1.49e-6", "-6.47e-7", "3.82e-7"),
        (2, 2, 1, "-2.05e-8", "-9.23e-8", "-3.33e-7"),
        (2, 2, 2, "-1.54e-10", "-2.79e-9", "1.48e-8"),
        (3, 0, 0, "1.20e-5", "-8.85e-6", "-4.00e-6"),
        (3, 0, 1, "-1.56e-6", "-2.33e-6", "1.04e-6"),
        (3, 0, 2, "4.91e-8", "1.13e-7", "-3.20e-8"),
        (3, 1, 0, "-1.71e-6", "-1.12e-6", "4.26e-7"),
        (3, 1, 1, "6.70e-8", "-7.78e-9", "-3.22e-7"),
        (3, 1, 2, "-3.62e-9", "-5.30e-9", "1.39e-8"),
        (3, 2, 0, "4.73e-8", "3.64e-8", "8.00e-10"),
        (3, 2, 1, "-7.22e-10", "1.84e-9", "7.90e-9"),
        (3, 2, 2, "4.09e-11", "5.67e-11", "-3.73e-10"),
        (4, 0, 0, "-1.63e-7", "9.96e-9", "4.03e-8"),
        (4, 0, 1, "1.56e-8", "1.99e-8", "-7.89e-9"),
        (4, 0, 2, "-4.75e-10", "-9.03e-10", "2.08e-10"),
        (4, 1, 0, "1.67e-8", "1.36e-8", "-1.46e-9"),
        (4, 1, 1, "-9.27e-10", "-3.46e-10", "2.45e-9"),
        (4, 1, 2, "4.15e-11", "4.95e-11", "-1.07e-10"),
        (4, 2, 0, "-4.05e-10", "-3.58e-10", "-5.47e-11"),
        (4, 2, 1, "8.76e-12", "-8.53e-12", "-6.05e-11"),
        (4, 2, 2, "-3.92e-13", "-4.38e-13", "2.95e-12"),
    ]
    expected = {(i, j, k): (sd, sl, qm) for i, j, k, sd, sl, qm in table}
    rows = coefficient_rows()
    indices = [tuple(int(row[key]) for key in ("i", "j", "k")) for row in rows]
    if len(rows) != len(expected) or set(indices) != set(expected):
        raise ValueError("Coefficient CSV has missing, duplicate or unexpected indices")
    for index, row in zip(indices, rows):
        for column, value in zip(("drag_area", "lift_area", "moment_volume"), expected[index]):
            if Decimal(row[column]).as_tuple() != Decimal(value).as_tuple():
                raise ValueError(f"Coefficient or printed precision mismatch: {index}, {column}")
    return len(table) * 3


def evaluate(alpha_deg, theta_deg=10, opening_deg=25):
    """Evaluate the model with the degree units used in the plots."""
    return aerodynamic_properties(
        math.radians(alpha_deg), math.radians(theta_deg), math.radians(opening_deg),
    )


def save_figure(figure, output, stem):
    figure.savefig(output / f"{stem}.png", dpi=180, bbox_inches="tight")
    figure.savefig(output / f"{stem}.pdf", bbox_inches="tight")
    plt.close(figure)


def write_curve(output, alphas, properties):
    with (output / "figure4_curve.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("alpha_deg", "theta_deg", "ski_opening_deg",
                         "drag_area_m2", "lift_area_m2", "moment_volume_m3"))
        for alpha, values in zip(alphas, properties):
            writer.writerow((alpha, 10, 25, values.drag_area_m2,
                             values.lift_area_m2, values.moment_volume_m3))


def plot_fits(output, alphas, properties):
    figure, axis = plt.subplots(figsize=(9, 5))
    moment_axis = axis.twinx()
    lines = axis.plot(alphas, [p.drag_area_m2 for p in properties],
                      color="#2369a1", label="Drag area $S_D$")
    lines += axis.plot(alphas, [p.lift_area_m2 for p in properties],
                       color="#28845a", linestyle="--", label="Lift area $S_L$")
    lines += moment_axis.plot(alphas, [p.moment_volume_m3 for p in properties],
                              color="#a5444e", label="Moment volume $Q_M$")
    axis.set(xlabel="Angle of attack α (degrees)", ylabel="Area (m²)",
             xlim=(0, 50), ylim=(-0.5, 1),
             title="Published polynomial fits · θ = 10°, λ = 25°")
    moment_axis.set(ylabel="Moment volume (m³)", ylim=(-0.2, 0.4))
    moment_axis.axhline(0, color="gray", linewidth=0.7)
    axis.grid(alpha=0.2)
    axis.legend(lines, [line.get_label() for line in lines], loc="upper left")
    figure.text(0.1, -0.02, "Comparison with Seo et al. (2004), Fig. 4; curves use printed Tables 1–3.", fontsize=9)
    save_figure(figure, output, "figure4_polynomial_fits")


def plot_polar(output, alphas):
    figure, axis = plt.subplots(figsize=(8, 5))
    for opening in (0, 10, 25):
        curve = [evaluate(alpha, opening_deg=opening) for alpha in alphas]
        axis.plot([p.drag_area_m2 for p in curve], [p.lift_area_m2 for p in curve],
                  label=f"λ = {opening}°")
    axis.set(xlabel="Drag area $S_D$ (m²)", ylabel="Lift area $S_L$ (m²)",
             title="Drag–lift polar · θ = 10° · α = 0–50°")
    axis.grid(alpha=0.2)
    axis.legend()
    save_figure(figure, output, "figure5_drag_lift_polar")


def plot_forward_lean(output, alphas):
    figure, axes = plt.subplots(2, 1, figsize=(8, 8), sharex=True)
    for theta in (0, 10, 20, 30, 40):
        curve = [evaluate(alpha, theta_deg=theta, opening_deg=10) for alpha in alphas]
        axes[0].plot(alphas, [p.drag_area_m2 for p in curve], label=f"θ = {theta}°")
        axes[1].plot(alphas, [p.lift_area_m2 for p in curve], label=f"θ = {theta}°")
    axes[0].set(ylabel="Drag area $S_D$ (m²)", title="Forward-lean comparison · λ = 10°")
    axes[1].set(xlabel="Angle of attack α (degrees)", ylabel="Lift area $S_L$ (m²)")
    for axis in axes:
        axis.grid(alpha=0.2)
        axis.legend()
    save_figure(figure, output, "figure6_forward_lean")


def find_trim_angle():
    """Find the zero moment between 20° and 27° (theta=10°, lambda=25°)."""
    low, high = 20.0, 27.0
    for _ in range(60):
        middle = (low + high) / 2
        if evaluate(middle).moment_volume_m3 > 0:
            low = middle
        else:
            high = middle
    return (low + high) / 2


def write_report(output, checked, alphas, properties, trim):
    peak_alpha, peak = max(zip(alphas, properties), key=lambda row: row[1].lift_area_m2)
    endpoint = evaluate(50)
    report = f"""# Aerodynamic verification

All {checked} CSV coefficients match the reference copy of Tables 1–3,
including their printed precision (Decimal comparison).

At θ = 10° and λ = 25°:

- Lift peaks at α = {peak_alpha:.1f}°, with S_L = {peak.lift_area_m2:.6f} m²
  (the paper describes stall near 37°).
- Zero pitching moment occurs at α = {trim:.4f}° (the paper gives about 24°).
- At α = 50°, Q_M = {endpoint.moment_volume_m3:.6f} m³. Figure 4 appears to
  show about −0.12 m³; that value is a visual estimate.

These curves use the rounded polynomial coefficients, not the wind-tunnel
measurements. The [precision analysis](coefficient-investigation.md) shows
that coefficient rounding can explain the moment discrepancy. The actual
unrounded values are unknown, and regression error may also contribute.

The paper reports fit errors of ±4% for S_D and S_L and ±9% for Q_M.
We cannot check those bounds without the original measurements.

## Output

- [Figure 4 polynomial fits](figure4_polynomial_fits.png)
- [Figure 5 drag–lift polar](figure5_drag_lift_polar.png)
- [Figure 6 forward-lean comparison](figure6_forward_lean.png)
- [Figure 4 numerical samples](figure4_curve.csv)

Each plot is also saved as a PDF.
"""
    (output / "verification.md").write_text(report, encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT,
                        help="Directory for plots, CSV and report")
    args = parser.parse_args()

    checked = audit_coefficients()
    args.output.mkdir(parents=True, exist_ok=True)
    alphas = [index / 10 for index in range(501)]
    properties = [evaluate(alpha) for alpha in alphas]
    trim = find_trim_angle()

    plt.rcParams.update({"font.size": 11, "axes.spines.top": False})
    write_curve(args.output, alphas, properties)
    plot_fits(args.output, alphas, properties)
    plot_polar(args.output, alphas)
    plot_forward_lean(args.output, alphas)
    write_report(args.output, checked, alphas, properties, trim)

    print(f"Audited {checked} coefficients; trim {trim:.4f}°")
    print(f"Results: {args.output}")


if __name__ == "__main__":
    main()
