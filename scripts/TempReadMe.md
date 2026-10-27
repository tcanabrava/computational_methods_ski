# Ski-jump project: how the code is organised

## Files (in scripts/)

| File | What it is | Who edits it |
|---|---|---|
| `seo2004_coefficients.csv` | The paper's wind-tunnel coefficients | Tomaz |
| `tomaz_aerodynamics.py` | The air: angles in, lift/drag/moment out | Tomaz |
| `make_plots.py` | The solver: equations, RK4, hill, landing detection, figures | Shahmeer |
| `jump.ipynb` | Runs the solver: validation and experiments | Mukhamed (section 5) |
| `figures/` | Generated plots, not code |

Each file only uses the one above it: notebook → solver → aerodynamics.

## Running it

Open `jump.ipynb` and run all cells. Sections 1–4 check the solver (RK4, no-air
parabola, reference jump, time-step convergence); don't change them.

## Doing experiments (section 5)

One jump:

    mp.landing(theta=np.radians(7), lam=np.radians(25), W=0.0)

returns the horizontal landing distance in metres.

- Angles in **radians**, so wrap degrees in `np.radians(...)`
- Wind `W` in m/s, **positive = tail wind**, negative = head wind
- Measured range of the aero data: θ 0–40°, λ 0–25°. Outside it, results are extrapolated.

To change a global setting for your runs, set it from the notebook, e.g.
`mp.m = 70.0`. **Don't edit `make_plots.py`.** If you need a new function
(e.g. posture changing during flight), write it in the notebook first and we'll
move it into the solver together.

## Mukhamed's tasks

1. Effect of wind, λ and θ on landing distance: one plot each, with the physical reason
2. Best fixed posture (θ, λ): start from the grid in the notebook, then refine
3. Posture changing during flight (paper: 139 m by varying λ)
4. Compare our numbers with the paper (131 m; wind 135/131/127 m)
