"""Investigate printed-coefficient rounding separately from float arithmetic.

Run: python -m scripts.coefficient_precision_report --source-pdf /path/to/Seo_2004_Aerodynamic_force_data.pdf
No baseline coefficients are changed by this analysis.
"""

import argparse
import csv
from decimal import Decimal, localcontext
from math import radians
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scripts.aerodynamic_report import audit_coefficients
from scripts.tomaz_aerodynamics import aerodynamic_properties, coefficient_rows

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results" / "aerodynamics"
FIELDS = ("drag_area", "lift_area", "moment_volume")


def source_rows():
    return coefficient_rows()


def decimal_evaluation(rows, alpha, theta, opening, field):
    """60-digit direct sum and worst-case error from nearest-digit rounding.

    A printed 1.39e-8 represents an unknown value within 0.005e-8 of
    that value, if it was rounded to nearest. Decimal retains trailing zeros.
    This bound applies to loss of printed coefficient precision, not to the
    regression error, measurement uncertainty, or floating-point rounding.
    """
    with localcontext() as context:
        context.prec = 60
        value = Decimal(0)
        bound = Decimal(0)
        absolute_terms = Decimal(0)
        for row in rows:
            i, j, k = (int(row[key]) for key in ("i", "j", "k"))
            # Define x^0 = 1 explicitly; Decimal otherwise rejects 0^0.
            power = ((alpha ** i if i else Decimal(1)) *
                     (theta ** j if j else Decimal(1)) *
                     (opening ** k if k else Decimal(1)))
            coefficient = Decimal(row[field])
            term = coefficient * power
            last_digit = Decimal(10) ** coefficient.as_tuple().exponent
            value += term
            absolute_terms += abs(term)
            bound += last_digit / 2 * abs(power)
        return value, bound, absolute_terms


def main():
    global OUTPUT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-pdf", type=Path, required=True, help="The supplied Seo aerodynamic-force paper")
    parser.add_argument("--output", type=Path, default=OUTPUT, help="Generated figures and report directory")
    args = parser.parse_args()
    OUTPUT = args.output
    audit_coefficients(args.source_pdf)
    rows = source_rows()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    table = []
    max_error = 0.0
    # Cover the measured envelope, not only the representative Figure 4 slice.
    for alpha in range(51):
        for theta in (0, 10, 20, 30, 40):
            for opening in (0, 10, 25):
                model = aerodynamic_properties(*(radians(v) for v in
                                                  (alpha, theta, opening)))
                floats = (model.drag_area_m2, model.lift_area_m2, model.moment_volume_m3)
                for field, actual in zip(FIELDS, floats):
                    exact, _, _ = decimal_evaluation(
                        rows, Decimal(alpha), Decimal(theta), Decimal(opening), field)
                    max_error = max(max_error, abs(actual - float(exact)))
    if max_error > 1e-10:
        raise AssertionError(f"Unexpected arithmetic discrepancy: {max_error}")

    for index in range(501):
        alpha = Decimal(index) / 10
        for field in FIELDS:
            value, bound, terms = decimal_evaluation(
                rows, alpha, Decimal(10), Decimal(25), field)
            table.append((alpha, field, value, bound, terms))
    with (OUTPUT / "coefficient_precision.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("alpha_deg", "property", "decimal_value",
                         "worst_case_coefficient_rounding_bound", "sum_absolute_terms"))
        writer.writerows(table)

    # Illustrative endpoint witness, not a correction or fitted calibration.
    altered = [dict(row) for row in rows]
    for row in altered:
        if (row["i"], row["j"], row["k"]) == ("3", "1", "2"):
            row["moment_volume"] = "1.38501e-8"
    witness, _, _ = decimal_evaluation(
        altered, Decimal(50), Decimal(10), Decimal(25), "moment_volume")
    assert format(Decimal("1.38501e-8"), ".2E") == format(Decimal("1.39e-8"), ".2E")

    figure, axes = plt.subplots(3, 1, figsize=(9, 10), sharex=True)
    for axis, field, label in zip(axes, FIELDS,
                                 ("Drag area (m²)", "Lift area (m²)", "Moment volume (m³)")):
        values = [row for row in table if row[1] == field]
        x = [float(row[0]) for row in values]
        nominal = [float(row[2]) for row in values]
        lower = [float(row[2] - row[3]) for row in values]
        upper = [float(row[2] + row[3]) for row in values]
        axis.fill_between(x, lower, upper, color="#2369a1", alpha=0.18,
                          label="Worst-case printed-coefficient rounding envelope")
        axis.plot(x, nominal, color="#2369a1", label="Printed coefficients")
        axis.set_ylabel(label)
        axis.grid(alpha=0.2)
    axes[0].set_title("Coefficient precision sensitivity · θ = 10°, λ = 25°")
    axes[0].legend(fontsize=8)
    axes[2].axhline(0, color="gray", linewidth=0.7)
    axes[2].set_xlabel("Angle of attack α (degrees)")
    figure.text(0.1, 0.02, "Assumes rounding to nearest printed digit. This is not a confidence interval.", fontsize=9)
    figure.tight_layout(rect=(0, 0.04, 1, 1))
    figure.savefig(OUTPUT / "coefficient_precision.png", dpi=180)
    figure.savefig(OUTPUT / "coefficient_precision.pdf")
    plt.close(figure)

    selected = []
    for alpha in (0, 10, 24, 37, 50):
        value, bound, terms = decimal_evaluation(
            rows, Decimal(alpha), Decimal(10), Decimal(25), "moment_volume")
        selected.append(f"| {alpha} | {value:.8f} | {bound:.8f} | {terms:.8f} |")
    report = f"""# Why the printed coefficients differ from Figure 4

Reproduce with `python -m scripts.coefficient_precision_report --source-pdf /path/to/Seo_2004_Aerodynamic_force_data.pdf`.

## Finding

**Loss of coefficient precision during publication is quantitatively sufficient to explain the large high-angle moment discrepancy.** This is a supported explanation, not confirmation of the authors' actual full-precision coefficients. Ordinary floating-point error in our implementation is far too small to explain it. Differences between a global polynomial fit and the particular measured Figure 4 slice may contribute as well.

## Transcription, units, axes, and arithmetic

- All 135 embedded coefficients match the PDF's Tables 1–3, including signs, exponents, and exponent indices.
- Equations (4)–(6) explicitly use α^i θ^j λ^k with angles in degrees; the module implements this ordering and conversion.
- Figure 4 uses the left axis for S_D and S_L and the right axis for Q_M. The right-axis endpoint is visually about −0.12 m³, whereas our table evaluation is −0.08009125 m³. That graph value is approximate, not a digitized raw measurement.
- An independent 60-digit Decimal direct sum was compared with the module's float evaluation at 765 angle combinations covering α = 0–50°, θ = 0,10,20,30,40°, and λ = 0,10,25°. The maximum absolute difference across the three properties was {max_error:.3e} in the corresponding property units. This is many orders of magnitude below the observed discrepancy.
- Density and speed are irrelevant to this comparison: Q_M is already moment divided by dynamic pressure, and the reference height scale is one.

## Rounding amplified by cancellation

Tables 1–3 print three significant digits. Under rounding to nearest, a printed coefficient c has an unknown rounding error no greater than half its last printed place. Propagating that loss of precision gives the bound:

```text
|δF| ≤ sum(0.5 × last_printed_place(c_ijk) × |α^i θ^j λ^k|)
```

For Q_M at θ = 10°, λ = 25°:

| α (degrees) | Printed-table Q_M (m³) | Worst-case rounding bound (m³) | Sum of absolute term magnitudes (m³) |
| --- | --- | --- | --- |
{chr(10).join(selected)}

At α = 50°, terms of combined absolute magnitude **81.24153375 m³** cancel to produce **−0.08009125 m³**. Their ratio is about 1,014. Three significant digits are insufficient to guarantee an accurate small residual after this cancellation. Using a more accurate summation algorithm cannot recover digits that were omitted from the source coefficients.

The endpoint mismatch of about 0.04 m³ is inside the ±0.19051875 m³ worst-case coefficient-rounding bound. At α = 24°, the bound contains zero, so coefficient rounding can also move the trim away from the printed-table value of 22.5147°.

These are conservative bounds under independent rounding assumptions. They are **not probability distributions, confidence intervals, or an estimate that the actual error equals the bound**. Different points on the plotted envelope can correspond to different coefficient perturbations; the envelope is not one alternative trajectory or one calibrated fit.

## A concrete endpoint example

The paper prints c_312 = 1.39 × 10⁻⁸. A hypothetical full-precision value of **1.38501 × 10⁻⁸** rounds to the same printed value, but differs by only about 0.36%.

At α = 50°, θ = 10°, λ = 25°, its multiplier is 50³ × 10 × 25² = 781,250,000. Changing only that unprinted portion produces **Q_M = {witness:.9f} m³**, close to the approximately −0.12 m³ figure endpoint.

This demonstrates that publication rounding alone can account for that endpoint. It does **not** identify the actual original coefficient or reproduce the whole measured curve. The embedded baseline coefficients remain unchanged.

## Other contributions and remaining uncertainty

Figure 4 presents experimental data, whereas the polynomials are global fits over α, θ, and λ. The paper reports regression discrepancies, and says actual α and θ were measured from photographs to account for ski bending (p. 34). Measurement scatter and nominal versus actual angles can therefore contribute. No exact per-point dataset is supplied, so their contribution cannot be quantified here.

A targeted search for an erratum/correction and a check of the [publisher article page](https://link.springer.com/article/10.1007/BF02843971) did not identify a correction or a full-precision coefficient supplement. This limited search does not establish that none exists.

## Decision for the simulation

Retain the published coefficients as the reproducible baseline. Record their precision limitation, especially for pitching dynamics. Obtaining full-precision coefficients or raw wind-tunnel measurements is the way to establish the true fit. Do not silently adjust coefficients to match a visually estimated graph point. Later sensitivity analysis should use fixed perturbed coefficient sets throughout each run, rather than treating this pointwise envelope as independent time-varying noise.

- [Rounding sensitivity plot](coefficient_precision.png)
- [Numerical values and bounds](coefficient_precision.csv)
- [First-week presentation](../../tomaz-first-week.typ)
"""
    (OUTPUT / "coefficient-investigation.md").write_text(report)
    print(f"765 angle combinations checked; max float/Decimal difference {max_error:.3e}")
    print(f"Hypothetical same-rounded c312 gives endpoint Q_M = {witness}")
    print(f"Wrote precision analysis to {OUTPUT}")


if __name__ == "__main__":
    main()
