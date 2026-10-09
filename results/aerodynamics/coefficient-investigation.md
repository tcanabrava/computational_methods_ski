# Coefficient precision

Reproduce with `python -m scripts.coefficient_precision_report`.

The CSV matches the reference copy of Tables 1–3. The polynomial uses degrees
and the exponent order α^i θ^j λ^k from equations (4)–(6).

## Arithmetic check

Float evaluation was compared with a 60-digit Decimal sum at 765 postures:
α = 0–50° in 1° steps, θ = 0, 10, 20, 30, 40°, and λ = 0, 10, 25°.
The largest absolute difference across S_D, S_L and Q_M was 7.744e-15
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
| 0 | 0.09049250 | 0.00097875 | 0.44399250 |
| 10 | 0.06437555 | 0.00671535 | 3.40574695 |
| 24 | -0.00739017 | 0.03218394 | 14.85427359 |
| 37 | -0.05661884 | 0.08739661 | 38.40006816 |
| 50 | -0.08009125 | 0.19051875 | 81.24153375 |

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
Changing only this coefficient gives Q_M = -0.119075625 m³,
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
