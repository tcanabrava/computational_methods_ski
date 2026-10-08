"""Seo, Watanabe & Murakami (2004), equations (4)–(6) and Tables 1–3.

Pass angles in radians; the paper's polynomials use degrees.
Run the reference checks with: python -m scripts.tomaz_aerodynamics
"""

import csv
import math
import unittest
from dataclasses import dataclass
from math import atan2, cos, degrees, fsum, hypot, isfinite, sin
from pathlib import Path

REFERENCE_HEIGHT_M = 1.76
MEASURED_ANGLE_RANGES_DEG = ((0.0, 50.0), (0.0, 40.0), (0.0, 25.0))
MEASURED_SPEED_RANGE_M_S = (20.0, 25.0)
COEFFICIENT_CSV = Path(__file__).with_name("seo2004_coefficients.csv")


def coefficient_rows():
    """Return the coefficients as strings, without losing printed precision."""
    with COEFFICIENT_CSV.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def _load_coefficients():
    rows = coefficient_rows()
    indices = [tuple(int(row[key]) for key in ("i", "j", "k")) for row in rows]
    expected = {(i, j, k) for i in range(5) for j in range(3) for k in range(3)}
    if len(indices) != 45 or set(indices) != expected:
        raise ValueError("Aerodynamic coefficient table is incomplete or duplicated")
    coefficients = []
    for index, row in zip(indices, rows):
        values = (float(row["drag_area"]), float(row["lift_area"]),
                  float(row["moment_volume"]))
        coefficients.append((index, values))
    return tuple(coefficients)


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
    """Return S_D (m²), S_L (m²) and Q_M (m³) for the given posture.

    Out-of-range angles raise ValueError unless allow_extrapolation is set.
    Extrapolated results are flagged; the angles are never clipped.
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
    # Equations (4)–(6): sum c_ijk * alpha^i * theta^j * lambda^k.
    terms = []
    for (i, j, k), coefficients in _COEFFICIENTS:
        power = alpha ** i * theta ** j * opening ** k
        terms.append((power, coefficients))
    sd = fsum(power * coefficients[0] for power, coefficients in terms)
    sl = fsum(power * coefficients[1] for power, coefficients in terms)
    qm = fsum(power * coefficients[2] for power, coefficients in terms)

    # The paper scales areas with height squared and moments with height cubed.
    ratio = height_m / REFERENCE_HEIGHT_M
    return AerodynamicProperties(
        sd * ratio ** 2, sl * ratio ** 2,
        qm * ratio ** 3, within_domain,
    )


def aerodynamic_forces(air_speed_m_s, alpha_rad, theta_rad, ski_opening_rad,
                       *, air_density_kg_m3, height_m=REFERENCE_HEIGHT_M,
                       allow_extrapolation=False):
    """Multiply S_D, S_L and Q_M by dynamic pressure (0.5 * rho * U²).

    Positive moment is nose-up. Speeds outside the measured 20–25 m/s
    range are allowed, but within_measured_speed_range will be False.
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
    """alpha = phi - beta_r - theta, in radians."""
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
    """Return aerodynamic (Fx, Fy) in newtons, with Y upward.

    Drag opposes air-relative motion; positive lift is the perpendicular
    direction (-sin(beta_r), cos(beta_r)). Gravity is separate.
    """
    _require_finite(air_path_angle_rad=air_path_angle_rad)
    c, s = cos(air_path_angle_rad), sin(air_path_angle_rad)
    return (-forces.drag_n * c - forces.lift_n * s,
            -forces.drag_n * s + forces.lift_n * c)


def posture(alpha=24, theta=10, opening=25):
    return tuple(math.radians(value) for value in (alpha, theta, opening))


class AerodynamicsTests(unittest.TestCase):
    def test_figure4_behaviour(self):
        # Figure 4: increasing drag, lift peak near 37°, restoring moment.
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
