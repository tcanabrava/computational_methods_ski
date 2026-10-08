"""Free-flight fits from Seo, Watanabe & Murakami (2004), pp. 34–36.

All public angles are radians. Only polynomial evaluation uses degrees.
Positive moment is nose-up. Run python -m scripts.tomaz_aerodynamics
to execute the nine reference checks. The coefficient data is embedded so
this first-week report does not depend on the flight application.
"""

import csv
from dataclasses import dataclass
from math import atan2, cos, degrees, fsum, hypot, isfinite, sin
from io import StringIO

REFERENCE_HEIGHT_M = 1.76
MEASURED_ANGLE_RANGES_DEG = ((0.0, 50.0), (0.0, 40.0), (0.0, 25.0))
MEASURED_SPEED_RANGE_M_S = (20.0, 25.0)


# Seo (2004), printed p. 36. Independently audited by aerodynamic_report.
COEFFICIENT_CSV = """i,j,k,drag_area,lift_area,moment_volume
0,0,0,8.54e-2,4.79e-2,7.87e-3
0,0,1,5.33e-4,-5.78e-3,-2.17e-3
0,0,2,-1.24e-4,2.63e-4,1.64e-4
0,1,0,7.90e-5,5.91e-3,2.52e-3
0,1,1,7.05e-5,2.57e-4,4.41e-4
0,1,2,9.07e-6,-9.12e-6,-1.58e-5
0,2,0,1.50e-4,2.79e-5,1.36e-5
0,2,1,-2.37e-6,-3.78e-6,-9.50e-6
0,2,2,-1.67e-7,-1.62e-8,3.21e-7
1,0,0,3.74e-4,4.78e-3,5.77e-5
1,0,1,-3.82e-4,-1.02e-3,5.39e-4
1,0,2,2.38e-5,7.87e-5,-2.63e-5
1,1,0,2.71e-4,6.50e-4,6.90e-4
1,1,1,-3.63e-5,-1.41e-5,-1.99e-4
1,1,2,-6.02e-8,-3.41e-6,7.54e-6
1,2,0,3.74e-6,-1.51e-5,-1.41e-5
1,2,1,9.61e-7,1.15e-6,4.74e-6
1,2,2,-9.33e-9,5.33e-8,-1.88e-7
2,0,0,4.53e-5,5.20e-4,5.16e-5
2,0,1,4.19e-5,9.20e-5,-3.86e-5
2,0,2,-1.38e-6,-5.02e-6,1.41e-6
2,1,0,4.25e-5,4.39e-6,-3.19e-5
2,1,1,6.17e-7,1.59e-6,1.35e-5
2,1,2,3.64e-8,2.02e-7,-5.61e-7
2,2,0,-1.49e-6,-6.47e-7,3.82e-7
2,2,1,-2.05e-8,-9.23e-8,-3.33e-7
2,2,2,-1.54e-10,-2.79e-9,1.48e-8
3,0,0,1.20e-5,-8.85e-6,-4.00e-6
3,0,1,-1.56e-6,-2.33e-6,1.04e-6
3,0,2,4.91e-8,1.13e-7,-3.20e-8
3,1,0,-1.71e-6,-1.12e-6,4.26e-7
3,1,1,6.70e-8,-7.78e-9,-3.22e-7
3,1,2,-3.62e-9,-5.30e-9,1.39e-8
3,2,0,4.73e-8,3.64e-8,8.00e-10
3,2,1,-7.22e-10,1.84e-9,7.90e-9
3,2,2,4.09e-11,5.67e-11,-3.73e-10
4,0,0,-1.63e-7,9.96e-9,4.03e-8
4,0,1,1.56e-8,1.99e-8,-7.89e-9
4,0,2,-4.75e-10,-9.03e-10,2.08e-10
4,1,0,1.67e-8,1.36e-8,-1.46e-9
4,1,1,-9.27e-10,-3.46e-10,2.45e-9
4,1,2,4.15e-11,4.95e-11,-1.07e-10
4,2,0,-4.05e-10,-3.58e-10,-5.47e-11
4,2,1,8.76e-12,-8.53e-12,-6.05e-11
4,2,2,-3.92e-13,-4.38e-13,2.95e-12
"""


def coefficient_rows():
    """Published Tables 1–3, retaining the printed significant digits."""
    return list(csv.DictReader(StringIO(COEFFICIENT_CSV)))


def _load_coefficients():
    rows = coefficient_rows()
    indices = [tuple(int(row[key]) for key in ("i", "j", "k")) for row in rows]
    expected = {(i, j, k) for i in range(5) for j in range(3) for k in range(3)}
    if len(indices) != 45 or set(indices) != expected:
        raise ValueError("Aerodynamic coefficient table is incomplete or duplicated")
    return tuple(
        (index, tuple(float(row[key]) for key in
                      ("drag_area", "lift_area", "moment_volume")))
        for index, row in zip(indices, rows)
    )


_COEFFICIENTS = _load_coefficients()


@dataclass(frozen=True)
class AerodynamicProperties:
    drag_area_m2: float
    lift_area_m2: float
    moment_volume_m3: float
    within_angle_domain: bool


@dataclass(frozen=True)
class AerodynamicForces:
    drag_n: float
    lift_n: float
    moment_nm: float
    properties: AerodynamicProperties
    within_measured_speed_range: bool


def _require_finite(**values):
    for name, value in values.items():
        if not isfinite(value):
            raise ValueError(f"{name} must be finite")


def aerodynamic_properties(alpha_rad, theta_rad, ski_opening_rad, *,
                           height_m=REFERENCE_HEIGHT_M,
                           allow_extrapolation=False):
    """Evaluate S_D [m²], S_L [m²], Q_M [m³] at a posture.

    Reject angles outside the experimental envelope by default. Explicit
    extrapolation returns within_angle_domain=False without clipping angles
    or outputs. Height scaling is the paper's approximate similarity rule.
    """
    _require_finite(alpha_rad=alpha_rad, theta_rad=theta_rad,
                    ski_opening_rad=ski_opening_rad, height_m=height_m)
    if height_m <= 0:
        raise ValueError("height_m must be positive")
    angles = tuple(degrees(value) for value in
                   (alpha_rad, theta_rad, ski_opening_rad))
    within_domain = all(lower <= angle <= upper for angle, (lower, upper)
                        in zip(angles, MEASURED_ANGLE_RANGES_DEG))
    if not within_domain and not allow_extrapolation:
        raise ValueError(
            f"Angles {angles} degrees outside measured envelope "
            "(alpha 0–50, theta 0–40, ski opening 0–25); "
            "set allow_extrapolation=True to evaluate and flag them"
        )
    alpha, theta, opening = angles
    terms = [(alpha ** i * theta ** j * opening ** k, coefficients)
             for (i, j, k), coefficients in _COEFFICIENTS]
    sd, sl, qm = (fsum(power * coefficients[column]
                      for power, coefficients in terms) for column in range(3))
    ratio = height_m / REFERENCE_HEIGHT_M
    return AerodynamicProperties(sd * ratio ** 2, sl * ratio ** 2,
                                 qm * ratio ** 3, within_domain)


def aerodynamic_forces(air_speed_m_s, alpha_rad, theta_rad, ski_opening_rad,
                       *, air_density_kg_m3, height_m=REFERENCE_HEIGHT_M,
                       allow_extrapolation=False):
    """Apply dynamic pressure to the polynomial properties.

    Density is explicit. Speeds outside 20–25 m/s are allowed and flagged:
    coefficients were measured at those speeds, but flight may exceed them.
    """
    _require_finite(air_speed_m_s=air_speed_m_s,
                    air_density_kg_m3=air_density_kg_m3)
    if air_speed_m_s < 0 or air_density_kg_m3 <= 0:
        raise ValueError("Air speed must be nonnegative and density positive")
    properties = aerodynamic_properties(
        alpha_rad, theta_rad, ski_opening_rad, height_m=height_m,
        allow_extrapolation=allow_extrapolation,
    )
    pressure = 0.5 * air_density_kg_m3 * air_speed_m_s ** 2
    lower, upper = MEASURED_SPEED_RANGE_M_S
    return AerodynamicForces(
        pressure * properties.drag_area_m2,
        pressure * properties.lift_area_m2,
        pressure * properties.moment_volume_m3,
        properties, lower <= air_speed_m_s <= upper,
    )


def angle_of_attack(body_angle_rad, air_path_angle_rad, theta_rad):
    """alpha = phi - beta_r - theta; do not wrap or clip physical angles."""
    _require_finite(body_angle_rad=body_angle_rad,
                    air_path_angle_rad=air_path_angle_rad, theta_rad=theta_rad)
    return body_angle_rad - air_path_angle_rad - theta_rad


def air_relative_velocity(ground_vx_m_s, ground_vy_m_s, wind_m_s):
    """Return (U_r, beta_r) for constant horizontal wind, positive tailwind."""
    _require_finite(ground_vx_m_s=ground_vx_m_s,
                    ground_vy_m_s=ground_vy_m_s, wind_m_s=wind_m_s)
    vx = ground_vx_m_s - wind_m_s
    speed = hypot(vx, ground_vy_m_s)
    if speed == 0:
        raise ValueError("Air-relative flight-path angle is undefined at zero speed")
    return speed, atan2(ground_vy_m_s, vx)


def ground_velocity(air_speed_m_s, air_path_angle_rad, wind_m_s):
    """Return (dX/dt, dY/dt), with X forward and Y upward."""
    _require_finite(air_speed_m_s=air_speed_m_s,
                    air_path_angle_rad=air_path_angle_rad, wind_m_s=wind_m_s)
    if air_speed_m_s < 0:
        raise ValueError("air_speed_m_s must be nonnegative")
    return (air_speed_m_s * cos(air_path_angle_rad) + wind_m_s,
            air_speed_m_s * sin(air_path_angle_rad))


def aerodynamic_force_components(forces, air_path_angle_rad):
    """Return aerodynamic (Fx, Fy) [N] for the solver or visual force arrows.

    Drag opposes air-relative motion; positive lift is the perpendicular
    direction (-sin(beta_r), cos(beta_r)). Gravity is separate.
    """
    _require_finite(air_path_angle_rad=air_path_angle_rad)
    c, s = cos(air_path_angle_rad), sin(air_path_angle_rad)
    return (-forces.drag_n * c - forces.lift_n * s,
            -forces.drag_n * s + forces.lift_n * c)


# The nine original static aerodynamic reference checks.
import math
import unittest

def posture(alpha=24, theta=10, opening=25):
    return tuple(math.radians(value) for value in (alpha, theta, opening))


class AerodynamicsTests(unittest.TestCase):
    def test_figure4_behaviour(self):
        # Independent physical observations: lift peaks near 37 degrees,
        # drag increases, and pitching moment restores towards trim.
        points = [aerodynamic_properties(*posture(alpha=a))
                  for a in range(0, 51)]
        self.assertTrue(all(right.drag_area_m2 > left.drag_area_m2
                            for left, right in zip(points, points[1:])))
        peak = max(range(51), key=lambda a: points[a].lift_area_m2)
        self.assertTrue(35 <= peak <= 39)
        self.assertGreater(points[20].moment_volume_m3, 0)
        self.assertLess(points[27].moment_volume_m3, 0)
        # Paper states a trim angle of about 24 degrees. Rounded-table
        # reproduction gives 22.51 degrees; a 3-degree band is qualitative.
        self.assertGreater(points[21].moment_volume_m3, 0)
        self.assertLess(points[25].moment_volume_m3, 0)

    def test_dynamic_pressure_scaling(self):
        first = aerodynamic_forces(20, *posture(), air_density_kg_m3=1.18)
        faster = aerodynamic_forces(40, *posture(), air_density_kg_m3=1.18)
        denser = aerodynamic_forces(20, *posture(), air_density_kg_m3=2.36)
        for name in ("drag_n", "lift_n", "moment_nm"):
            self.assertAlmostEqual(getattr(faster, name), 4 * getattr(first, name))
            self.assertAlmostEqual(getattr(denser, name), 2 * getattr(first, name))
        self.assertTrue(first.within_measured_speed_range)
        self.assertFalse(faster.within_measured_speed_range)

    def test_zero_speed_has_zero_aerodynamic_loads(self):
        loads = aerodynamic_forces(0, *posture(), air_density_kg_m3=1.18)
        self.assertEqual((loads.drag_n, loads.lift_n, loads.moment_nm), (0, 0, 0))

    def test_force_vectors_have_correct_directions(self):
        forces = aerodynamic_forces(25, *posture(), air_density_kg_m3=1.18)
        for beta in (0, math.radians(-16.2), math.radians(120)):
            fx, fy = aerodynamic_force_components(forces, beta)
            # Projections onto flight tangent and its positive normal.
            self.assertAlmostEqual(fx * math.cos(beta) + fy * math.sin(beta),
                                   -forces.drag_n)
            self.assertAlmostEqual(-fx * math.sin(beta) + fy * math.cos(beta),
                                   forces.lift_n)

    def test_paper_height_similarity_rule(self):
        first = aerodynamic_properties(*posture())
        scaled = aerodynamic_properties(*posture(), height_m=2 * 1.76)
        self.assertAlmostEqual(scaled.drag_area_m2, 4 * first.drag_area_m2)
        self.assertAlmostEqual(scaled.lift_area_m2, 4 * first.lift_area_m2)
        self.assertAlmostEqual(scaled.moment_volume_m3, 8 * first.moment_volume_m3)

    def test_range_policy_does_not_silently_clip(self):
        for angles in (posture(alpha=-1), posture(alpha=51),
                       posture(theta=41), posture(opening=26)):
            with self.assertRaises(ValueError):
                aerodynamic_properties(*angles)
            result = aerodynamic_properties(*angles, allow_extrapolation=True)
            self.assertFalse(result.within_angle_domain)
        for angles in (posture(0, 0, 0), posture(50, 40, 25)):
            self.assertTrue(aerodynamic_properties(*angles).within_angle_domain)

    def test_invalid_physical_inputs(self):
        for invalid in (math.nan, math.inf, -math.inf):
            with self.assertRaises(ValueError):
                aerodynamic_properties(invalid, *posture()[1:])
        for density in (0, -1, math.nan):
            with self.assertRaises(ValueError):
                aerodynamic_forces(25, *posture(), air_density_kg_m3=density)
        with self.assertRaises(ValueError):
            aerodynamic_forces(-1, *posture(), air_density_kg_m3=1.18)
        with self.assertRaises(ValueError):
            aerodynamic_properties(*posture(), height_m=0)

    def test_angle_relation_at_published_initial_state(self):
        alpha = angle_of_attack(math.radians(20), math.radians(-16.2),
                                math.radians(7))
        self.assertAlmostEqual(math.degrees(alpha), 29.2)

    def test_wind_sign_and_velocity_roundtrip(self):
        # For the same ground velocity, headwind increases air speed.
        headwind, _ = air_relative_velocity(25, -5, -1)
        calm, _ = air_relative_velocity(25, -5, 0)
        tailwind, _ = air_relative_velocity(25, -5, 1)
        self.assertGreater(headwind, calm)
        self.assertGreater(calm, tailwind)
        for vx, vy, wind in ((25, -5, -1), (25, -5, 1), (-4, 3, 2)):
            speed, angle = air_relative_velocity(vx, vy, wind)
            restored = ground_velocity(speed, angle, wind)
            self.assertAlmostEqual(restored[0], vx)
            self.assertAlmostEqual(restored[1], vy)
        with self.assertRaises(ValueError):
            air_relative_velocity(2, 0, 2)


if __name__ == "__main__":
    unittest.main()
