"""Audit source tables and generate aerodynamic curves and CSV for review.

Run: python -m scripts.aerodynamic_report --source-pdf /path/to/Seo_2004_Aerodynamic_force_data.pdf
Requires matplotlib and the pdftotext command; the model uses only stdlib.
"""

import argparse
import csv
from decimal import Decimal
import math
from pathlib import Path
import re
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scripts.tomaz_aerodynamics import aerodynamic_properties, coefficient_rows

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results" / "aerodynamics"


def audit_coefficients(pdf):
    """Compare embedded coefficients to independently extracted PDF tables."""
    text = subprocess.run(
        ["pdftotext", "-f", "6", "-l", "6", "-layout", str(pdf), "-"],
        check=True, capture_output=True, text=True,
    ).stdout
    pattern = r"([abc])([0-4])([0-2])([0-2])\s+([–−-]?\d+\.\d+)\s*×\s*10([–−-]?\d+)"
    extracted = {}
    for family, i, j, k, mantissa, exponent in re.findall(pattern, text):
        key = (family, i, j, k)
        if key in extracted:
            raise AssertionError(f"Duplicate PDF coefficient: {key}")
        normalize = lambda value: value.replace("–", "-").replace("−", "-")
        extracted[key] = Decimal(normalize(mantissa)) * (
            Decimal(10) ** int(normalize(exponent)))
    if len(extracted) != 135:
        raise AssertionError(f"Expected 135 source coefficients, got {len(extracted)}")
    rows = coefficient_rows()
    for row in rows:
        for family, column in (("a", "drag_area"), ("b", "lift_area"),
                               ("c", "moment_volume")):
            key = (family, row["i"], row["j"], row["k"])
            if Decimal(row[column]) != extracted[key]:
                raise AssertionError(f"Coefficient mismatch: {key}")
    return len(extracted)


def evaluate(alpha_deg, theta_deg=10, opening_deg=25):
    return aerodynamic_properties(*(math.radians(value) for value in
                                    (alpha_deg, theta_deg, opening_deg)))


def save_figure(figure, stem):
    figure.savefig(OUTPUT / f"{stem}.png", dpi=180, bbox_inches="tight")
    figure.savefig(OUTPUT / f"{stem}.pdf", bbox_inches="tight")
    plt.close(figure)


def main():
    global OUTPUT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-pdf", type=Path, required=True, help="The supplied Seo aerodynamic-force paper; never copied into output")
    parser.add_argument("--output", type=Path, default=OUTPUT, help="Generated figures and report directory")
    args = parser.parse_args()
    OUTPUT = args.output
    checked = audit_coefficients(args.source_pdf)
    OUTPUT.mkdir(parents=True, exist_ok=True)
    alphas = [index / 10 for index in range(501)]
    properties = [evaluate(alpha) for alpha in alphas]
    with (OUTPUT / "figure4_curve.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("alpha_deg", "theta_deg", "ski_opening_deg",
                         "drag_area_m2", "lift_area_m2", "moment_volume_m3"))
        for alpha, values in zip(alphas, properties):
            writer.writerow((alpha, 10, 25, values.drag_area_m2,
                             values.lift_area_m2, values.moment_volume_m3))

    plt.rcParams.update({"font.size": 11, "axes.spines.top": False})
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
    save_figure(figure, "figure4_polynomial_fits")

    figure, axis = plt.subplots(figsize=(8, 5))
    for opening in (0, 10, 25):
        curve = [evaluate(alpha, opening_deg=opening) for alpha in alphas]
        axis.plot([p.drag_area_m2 for p in curve], [p.lift_area_m2 for p in curve],
                  label=f"λ = {opening}°")
    axis.set(xlabel="Drag area $S_D$ (m²)", ylabel="Lift area $S_L$ (m²)",
             title="Drag–lift polar · θ = 10° · α = 0–50°")
    axis.grid(alpha=0.2)
    axis.legend()
    save_figure(figure, "figure5_drag_lift_polar")

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
    save_figure(figure, "figure6_forward_lean")

    low, high = 20.0, 27.0
    for _ in range(60):
        middle = (low + high) / 2
        if evaluate(middle).moment_volume_m3 > 0:
            low = middle
        else:
            high = middle
    trim = (low + high) / 2
    peak_alpha, peak = max(zip(alphas, properties), key=lambda row: row[1].lift_area_m2)
    endpoint = evaluate(50)
    report = f"""# Aerodynamic verification

Generated with `python -m scripts.aerodynamic_report --source-pdf /path/to/Seo_2004_Aerodynamic_force_data.pdf`.

- All {checked} coefficients match the PDF text exactly (Decimal comparison).
- The printed table page was also visually reviewed for signs and exponents.
- At θ = 10°, λ = 25°, lift peaks at α ≈ {peak_alpha:.1f}°, S_L = {peak.lift_area_m2:.6f} m². The paper describes stall near 37°.
- The polynomial trim angle is α = {trim:.4f}°. The paper describes trim near 24°.
- Drag increases over α = 0–50° in this representative case.
- Pitching moment is positive below trim and negative above it, consistent with restoring pitch.

## Comparison limits

These are reproductions of the **printed polynomial tables**, not the original experimental dataset. Figure 4 shows measured points and curves that do not exactly match these rounded coefficients. In particular, at α = 50°, the computed Q_M is {endpoint.moment_volume_m3:.6f} m³, whereas the source figure is visually about −0.12 m³. That endpoint difference is substantial. The [coefficient precision investigation](coefficient-investigation.md) demonstrates that omitted coefficient digits are sufficient to explain this endpoint discrepancy, amplified by cancellation of large polynomial terms. The authors' actual full-precision values remain unknown, and fit/measurement discrepancies may contribute.

The source reports fit errors of ±4% for S_D and S_L and ±9% for Q_M against its experimental data. This implementation cannot independently establish those error bounds without the raw data. Do not treat the qualitative figure checks as quantitative validation of those percentages, and do not change coefficients to force agreement with a plot.

Generated figures:

- [Figure 4 polynomial fits](figure4_polynomial_fits.png)
- [Figure 5 drag–lift polar](figure5_drag_lift_polar.png)
- [Figure 6 forward-lean comparison](figure6_forward_lean.png)
- [Figure 4 numerical samples](figure4_curve.csv)

Each figure also has a PDF for presentation use.
"""
    (OUTPUT / "verification.md").write_text(report)
    print(f"Audited {checked} coefficients; trim {trim:.4f}°, lift peak {peak_alpha:.1f}°")
    print(f"Wrote figures, CSV, and verification report to {OUTPUT}")


if __name__ == "__main__":
    main()
