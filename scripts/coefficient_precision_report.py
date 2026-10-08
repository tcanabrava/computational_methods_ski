"""Compare float arithmetic with Decimal and estimate coefficient rounding.

Run: python -m scripts.coefficient_precision_report
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


def decimal_evaluation(rows, alpha, theta, opening, field):
    """Return the 60-digit sum, coefficient rounding bound and absolute sum.

    If rounded to nearest, 1.39e-8 is within 0.005e-8 of the original value.
    Decimal keeps the printed precision, including trailing zeros. The bound
    covers coefficient rounding, not measurement or regression error.
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


def check_arithmetic(rows):
    """Compare float and Decimal at 765 postures across the measured range."""
    max_error = 0.0
    for alpha in range(51):
        for theta in (0, 10, 20, 30, 40):
            for opening in (0, 10, 25):
                model = aerodynamic_properties(radians(alpha), radians(theta), radians(opening))
                floats = (model.drag_area_m2, model.lift_area_m2, model.moment_volume_m3)
                for field, actual in zip(FIELDS, floats):
                    exact, _, _ = decimal_evaluation(
                        rows, Decimal(alpha), Decimal(theta), Decimal(opening), field)
                    max_error = max(max_error, abs(actual - float(exact)))
    if max_error > 1e-10:
        raise AssertionError(f"Unexpected arithmetic discrepancy: {max_error}")
    return max_error


def precision_table(rows):
    table = []
    for index in range(501):
        alpha = Decimal(index) / 10
        for field in FIELDS:
            value, bound, terms = decimal_evaluation(
                rows, alpha, Decimal(10), Decimal(25), field)
            table.append((alpha, field, value, bound, terms))
    return table


def write_csv(output, table):
    with (output / "coefficient_precision.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("alpha_deg", "property", "decimal_value",
                         "worst_case_coefficient_rounding_bound", "sum_absolute_terms"))
        writer.writerows(table)


def endpoint_example(rows):
    """Try an unprinted value of c_312 that rounds to the published 1.39e-8."""
    altered = [dict(row) for row in rows]
    for row in altered:
        if (row["i"], row["j"], row["k"]) == ("3", "1", "2"):
            row["moment_volume"] = "1.38501e-8"
    witness, _, _ = decimal_evaluation(
        altered, Decimal(50), Decimal(10), Decimal(25), "moment_volume")
    assert format(Decimal("1.38501e-8"), ".2E") == format(Decimal("1.39e-8"), ".2E")
    return witness


def plot_precision(output, table):
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
    figure.savefig(output / "coefficient_precision.png", dpi=180)
    figure.savefig(output / "coefficient_precision.pdf")
    plt.close(figure)


def write_report(output, rows, max_error, witness):
    selected = []
    for alpha in (0, 10, 24, 37, 50):
        value, bound, terms = decimal_evaluation(
            rows, Decimal(alpha), Decimal(10), Decimal(25), "moment_volume")
        selected.append(f"| {alpha} | {value:.8f} | {bound:.8f} | {terms:.8f} |")
    selected = "\n".join(selected)
    report = f"""# Coefficient precision

Reproduce with `python -m scripts.coefficient_precision_report`.

The CSV matches the reference copy of Tables 1–3. The polynomial uses degrees
and the exponent order α^i θ^j λ^k from equations (4)–(6).

## Arithmetic check

Float evaluation was compared with a 60-digit Decimal sum at 765 postures:
α = 0–50° in 1° steps, θ = 0, 10, 20, 30, 40°, and λ = 0, 10, 25°.
The largest absolute difference across S_D, S_L and Q_M was {max_error:.3e}
in the respective property units. Float arithmetic cannot explain the
roughly 0.04 m³ moment discrepancy at the Figure 4 endpoint.

## Printed coefficient rounding

Assuming rounding to nearest, each coefficient can differ from its original
value by half its last printed place. The resulting worst-case bound is:

```text
|δF| ≤ sum(0.5 × last_printed_place(c_ijk) × |α^i θ^j λ^k|)
```

For Q_M at θ = 10°, λ = 25°:

| α (degrees) | Q_M (m³) | Rounding bound (m³) | Sum of absolute terms (m³) |
| --- | --- | --- | --- |
{selected}

At α = 50°, terms with a combined absolute magnitude of 81.24153375 m³
cancel to give −0.08009125 m³. The small result is sensitive to omitted
digits. The approximate Figure 4 value of −0.12 m³ lies within the rounding
bound of ±0.19051875 m³. At α = 24°, the bound includes zero, so rounding
could also explain a trim angle different from the computed 22.5147°.

These bounds assume independent coefficient rounding. They are not
confidence intervals or estimates of the actual error. Each point on the
envelope may require a different set of coefficient changes; the envelope
is not a single alternative curve.

## Example at α = 50°

The paper prints c_312 = 1.39 × 10⁻⁸. The hypothetical value
1.38501 × 10⁻⁸ rounds to the same three significant digits.
Its multiplier at θ = 10°, λ = 25° is 50³ × 10 × 25² = 781,250,000.
Changing only this coefficient gives Q_M = {witness:.9f} m³,
close to the visually estimated −0.12 m³ endpoint.

This shows that rounding is sufficient to explain the endpoint difference.
It does not recover the original coefficient or validate the whole measured
curve. The example works on a copy; the CSV coefficients remain unchanged.

## What remains unknown

Figure 4 includes measurements; the polynomial is a global fit. Regression
error, measurement scatter and actual versus nominal angles may contribute.
The paper also notes angles measured from photographs to account for ski
bending (p. 34). Without the original coefficients or measurements, these
contributions cannot be separated.

Keep the published values as the baseline and record the precision limit,
especially for pitching moment. Any later sensitivity runs should use a
fixed perturbed coefficient set throughout each run.

- [Rounding plot](coefficient_precision.png)
- [Values and bounds](coefficient_precision.csv)
- [First-week presentation](../../tomaz-first-week.typ)
"""
    (output / "coefficient-investigation.md").write_text(report, encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT,
                        help="Directory for plots, CSV and report")
    args = parser.parse_args()

    audit_coefficients()
    rows = coefficient_rows()
    max_error = check_arithmetic(rows)
    table = precision_table(rows)
    witness = endpoint_example(rows)

    args.output.mkdir(parents=True, exist_ok=True)
    write_csv(args.output, table)
    plot_precision(args.output, table)
    write_report(args.output, rows, max_error, witness)

    print(f"765 angle combinations checked; max float/Decimal difference {max_error:.3e}")
    print(f"Hypothetical same-rounded c312 gives endpoint Q_M = {witness}")
    print(f"Results: {args.output}")


if __name__ == "__main__":
    main()
