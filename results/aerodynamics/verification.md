# Aerodynamic verification

All 135 CSV coefficients match the reference copy of Tables 1–3,
including their printed precision (Decimal comparison).

At θ = 10° and λ = 25°:

- Lift peaks at α = 36.7°, with S_L = 0.769157 m²
  (the paper describes stall near 37°).
- Zero pitching moment occurs at α = 22.5147° (the paper gives about 24°).
- At α = 50°, Q_M = -0.080091 m³. Figure 4 appears to
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
