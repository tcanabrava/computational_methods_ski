// First-week aerodynamic model report. All figures are static evaluations.
// Generate the figures from the repository root before compiling:
//   python -m scripts.tomaz_aerodynamics
//   python -m scripts.aerodynamic_report
//   python -m scripts.coefficient_precision_report
//   typst compile tomaz-first-week.typ /tmp/tomaz-first-week.pdf
// Python generators require matplotlib. Coefficients: scripts/seo2004_coefficients.csv.
// Figures and numerical reports are generated output, not versioned inputs.
#import "@preview/touying:0.8.0": *
#import themes.simple: *

#show: simple-theme.with(
  aspect-ratio: "16-9",
  footer: [Computational Methods · Week 1 · Tomaz Canabrava],
)
#set text(size: 18pt)
#set par(leading: 0.55em)
#set document(
  title: "First-week report: aerodynamic model and inputs",
  author: "Tomaz Canabrava",
)
#let source(body) = text(size: 10pt, fill: gray, body)
#let plot(path, height: 6cm) = align(center,
  image(path, width: 100%, height: height, fit: "contain"),
)
// Keep the discreet illustration inline so the report needs no extra asset.
#let pole-dance-svg = ```xml
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 100 190">
  <title>Pole dancer</title>
  <g fill="none" stroke="#ededed" stroke-linecap="round" stroke-linejoin="round">
    <path d="M72 10V181" stroke-width="2.5"/>
    <path d="M47 72Q43 88 48 96L43 112" stroke-width="6"/>
    <path d="M48 75L62 51L72 39M47 78L59 88L72 75" stroke-width="4.5"/>
    <path d="M43 111L24 128L10 163M44 112L62 133L80 115" stroke-width="5"/>
  </g>
  <circle cx="47" cy="58" r="8" fill="#ededed"/>
</svg>
```.text

#title-slide[
  #text(size: 30pt, weight: "bold")[Aerodynamic model for V-style ski jumping]
  #v(0.7em)
  #text(size: 22pt)[First-week report · Aerodynamics and model inputs]
  #v(1em)
  Tomaz Canabrava
  #v(1em)
  Formulation, input conventions, verification, and limitations.
]

== Objectives and scope

My first-week work concerns the aerodynamic component of the model:

- Define lift $L$, drag $D$, and pitching moment $M$ from published data.
- Establish their dependence on $theta$, $lambda$, and $alpha$.
- Specify how wind enters the air-relative velocity.

I implemented the published force fits and documented:

- The force-evaluation interface and coefficient source.
- Input/output units and sign conventions.
- Assumptions, supported ranges, and limitations.

#source[First-week scope: aerodynamic formulation, input conventions, verification, and limitations. All evaluations use prescribed aerodynamic inputs.]

== Posture and angle definitions

#table(columns: (0.5fr, 2.3fr), inset: 7pt,
  table.header([Angle], [Meaning]),
  [$theta$], [Body inclination relative to the ski plane: forward lean.],
  [$lambda$], [Full opening angle between the two skis: the V opening.],
  [$alpha$], [Ski-plane inclination relative to the air-relative motion.],
)

$ alpha = phi - beta_r - theta $

$phi$: body orientation; $beta_r$: air-relative path angle.\
Both are measured from the forward horizontal axis.

I use *radians at the interface* and *degrees in the published polynomials*.

#source[Seo, Aerodynamic force data, Fig. 2, printed p. 33; Seo, Optimal flight technique, printed p. 99. The model converts radians to degrees internally.]

== Representation of wind

I assume constant horizontal wind $U_w$ and calculate:

$ v_(r,x) = v_(g,x) - U_w, quad v_(r,y) = v_(g,y) $
$ U_r = sqrt(v_(r,x)^2 + v_(r,y)^2), quad beta_r = "atan2"(v_(r,y), v_(r,x)) $

- Positive $U_w$: tailwind; negative $U_w$: headwind.
- With ground velocity held fixed, headwind increases air-relative speed.
- Wind can also change $beta_r$, and therefore $alpha$.
- I account for wind through $U_r$ and $alpha$ in the force evaluation.

#source[Seo, Optimal flight technique, printed pp. 98–99. Constant horizontal wind is the selected model assumption.]

== Polynomial representation of aerodynamic properties

I evaluate the three fits published by Seo, Watanabe & Murakami:

$
  F(alpha, theta, lambda) =
  sum_(i=0)^4 sum_(j=0)^2 sum_(k=0)^2
  d_(i j k) alpha^i theta^j lambda^k
$

#table(columns: (0.65fr, 1.7fr, 0.7fr), inset: 7pt,
  table.header([Output $F$], [Meaning], [Coefficient family]),
  [$S_D$], [Drag area · m²], [$a_(i j k)$],
  [$S_L$], [Lift area · m²], [$b_(i j k)$],
  [$Q_M$], [Moment volume · m³], [$c_(i j k)$],
)

*45 terms per output; 135 published coefficients in total.*\
The arguments in this formula are numerical angles in degrees.

#source[Seo, Aerodynamic force data, Eqs. (4)–(6), printed p. 35; Tables 1–3, printed p. 36. d denotes the appropriate coefficient family.]

== Calculation of aerodynamic loads

$ q = frac(1, 2) rho U_r^2 $
$ D = q S_D, quad L = q S_L, quad M = q Q_M $

#table(columns: (1fr, 1fr, 1.2fr), inset: 7pt,
  table.header([Quantity], [Unit], [Interpretation]),
  [$rho$], [kg/m³], [Air-density input],
  [$q$], [Pa = N/m²], [Dynamic pressure],
  [$D, L$], [N], [Aerodynamic forces],
  [$M$], [N m], [Torque about the centre of gravity],
)

$S_D$, $S_L$, and $Q_M$ incorporate the reference geometry.\
I therefore apply dynamic pressure directly, without additional geometric factors.

#source[Seo, Aerodynamic force data, Eqs. (1)–(3), printed p. 34. Density is supplied explicitly; 1.18 kg/m³ on p. 38 belongs to the paper's simplified pitching example.]

== Force directions and moment convention

I adopt coordinates with $X$ forward and $Y$ upward.

$ bold(F)_D = -D (cos(beta_r), sin(beta_r)) $
$ bold(F)_L = L (-sin(beta_r), cos(beta_r)) $

- Drag opposes air-relative motion.
- Lift acts perpendicular to that motion, along the stated normal.
- Positive $M$: nose-up; negative $M$: nose-down.
- $M = 0$ defines a trim point for a prescribed posture.

#source[Seo, Aerodynamic force data, Fig. 1 and sign definition, printed p. 33; coordinate conventions from Optimal flight technique, printed pp. 98–99.]

== Experimental coverage and input policy

#table(columns: (0.6fr, 1.25fr, 1.5fr), inset: 7pt,
  table.header([Input], [Measured settings], [Model policy]),
  [$alpha$], [0–50°, nominal 5° spacing], [Accepted angle envelope],
  [$theta$], [0–40°, nominal 10° spacing], [Accepted angle envelope],
  [$lambda$], [0°, 10°, 25°], [Polynomial interpolation over 0–25°],
  [$U$], [20 and 25 m/s], [Separate speed-coverage flag],
)

I evaluate intermediate angles using the published fit.\
I reject angles outside the envelope by default and report speed coverage separately.

The authors' estimate at 35 m/s is based on an ellipsoid analogy.\
It is not a direct ski-jumper measurement at that speed.

#source[Seo, Aerodynamic force data, printed pp. 32–34. Actual α and θ were checked from photographs because the skis bend. Envelope membership alone does not establish accuracy.]

== Static evaluation at a fixed posture

#plot("results/aerodynamics/figure4_polynomial_fits.png", height: 5.8cm)

At $theta = 10°$, $lambda = 25°$: lift peaks near *36.7°*.\
Printed-fit trim: *22.5147°*; the paper describes approximately *24°*.

#source[Static evaluation of the printed polynomial fits. Seo, Aerodynamic force data, Fig. 4, printed p. 34; results/aerodynamics/verification.md. Original raw measurements are unavailable.]

== Influence of ski opening on lift and drag

#plot("results/aerodynamics/figure5_drag_lift_polar.png", height: 5.8cm)

I evaluate $alpha = 0$–$50°$ at fixed $theta = 10°$ for each opening.\
At larger $alpha$, V opening can give more lift for comparable drag.

#source[Static aerodynamic evaluations. Seo, Aerodynamic force data, Fig. 5 and wake explanation, printed p. 35. These results describe the lift–drag relation at prescribed postures.]

== Influence of forward lean

#plot("results/aerodynamics/figure6_forward_lean.png", height: 6.2cm)

I compare $theta$ at the same prescribed $alpha$, with $lambda = 10°$.

#source[Static evaluations of the printed fits, corresponding to Seo's Fig. 6, printed p. 35. Each curve holds the specified posture inputs fixed.]

== Verification of the aerodynamic implementation

- I checked *all 135 coefficients* against a reference transcription of the paper's tables.
- I verified signs, exponent indices, and degree conversion.
- *Nine aerodynamic reference tests pass*: scaling, directions, wind conversion,
  input validation, angle bounds, and representative curve behaviour.
- Independent high-precision arithmetic agrees with the reference evaluator
  within about *8 × 10⁻¹⁵* in the relevant property units on the checked grid.

These checks verify my evaluation of the published model.\
I cannot establish quantitative experimental agreement without the original data.

#source[Generated records: results/aerodynamics/verification.md and coefficient-investigation.md. Reference tests: scripts/tomaz_aerodynamics.py. Source-table audit: scripts/aerodynamic_report.py.]

== Sensitivity to published coefficient precision

At $alpha = 50°$, $theta = 10°$, $lambda = 25°$:

#table(columns: (1.8fr, 1fr), inset: 7pt,
  [Printed-table $Q_M$], [−0.080091 m³],
  [Approximate value read from Fig. 4], [−0.12 m³],
  [Sum of absolute moment-term magnitudes], [81.241534 m³],
)

Large positive and negative terms cancel to a small result.\
Missing coefficient digits can strongly affect that remainder.

My pointwise analysis shows that publication rounding can explain this difference.\
The authors' full-precision coefficients remain unavailable.\
*I retain the printed coefficients and record this uncertainty.*

#source[Pointwise coefficient audit, results/aerodynamics/coefficient-investigation.md. The graph value is approximate. High-precision arithmetic excludes ordinary floating-point error as the explanation at this scale.]

== Assumptions and limitations

// A discreet, non-textual easter egg in the lower-right background.
#place(top + right, dx: -1mm, dy: 76mm,
  image(bytes(pole-dance-svg), format: "svg", width: 12mm),
)

- Quasi-steady, free-flight aerodynamic fit.
- Reference athlete height *1.76 m*; ski length *2.52 m*.
- Straight torso/leg alignment; fixed arms at *170°*; ski tails touching.
- I supply density explicitly and retain the reference geometry of the fits.
- No angular-rate damping or general height-dependent ground correction.
- Printed coefficients and unavailable raw data limit quantitative agreement.

The paper reports fit errors of ±4% for areas and ±9% for moment volume.\
I have not independently verified those estimates against the original data.

#source[Seo, Aerodynamic force data, printed pp. 32–35 and ground-effect discussion pp. 37–38. The selected model reproduces its free-flight fits.]

== First-week deliverable and interface

*Inputs:* $U_r$, $alpha$, $theta$, $lambda$, and $rho$.\
When orientation and path angle are supplied, I calculate $alpha = phi - beta_r - theta$.

*Outputs:* $D$ in N, $L$ in N, $M$ in N m, plus validity information.

My deliverable consists of:

- Coefficient source and aerodynamic evaluator.
- Units, angle definitions, wind convention, and moment sign.
- Supported ranges, assumptions, and static verification figures.
- The recorded precision discrepancy and missing source data.

This week, I completed the aerodynamic formulation and its static verification.\
The model is prepared for connection to the flight solver in the following week.

#source[Reference evaluator: scripts/tomaz_aerodynamics.py. Coefficients: scripts/seo2004_coefficients.csv. Static report generators: scripts/aerodynamic_report.py and scripts/coefficient_precision_report.py. Sources: the two supplied Seo (2004) papers.]

== References

*Seo, Watanabe & Murakami (2004).*\
“Aerodynamic force data for a V-style ski jumping flight.”\
_Sports Engineering_ 7(1), 31–39.

Experimental inputs and angle definitions: pp. 32–34.\
Property definitions, fits, and coefficient tables: pp. 34–36.

#v(0.5em)
*Seo, Murakami & Yoshida (2004).*\
“Optimal flight technique for V-style ski jumping.”\
_Sports Engineering_ 7(2), 97–103, supplied scan.

I use pp. 98–99 for coordinates, air-relative wind convention, and angles.

#source[Source PDFs: Seo_2004_Aerodynamic_force_data.pdf and Seo_2004_Optimal_flight_technique.pdf. All page numbers above are printed journal pages. The scripts use the transcribed coefficient tables; the PDFs are not required to run them.]
