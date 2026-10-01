= Project 2: Numerical Modelling of Ski-Jumping Flight

Simulate the flight of a V-style ski jumper and determine how posture and wind affect the landing distance.

The jumper is modelled as a body–ski system moving in a vertical plane. During flight, the system is subject to gravity, aerodynamic lift, drag, and pitching moment.

The main posture parameters are the forward-lean angle $theta$ and ski-opening angle $lambda$. The aerodynamic forces also depend on the angle of attack $alpha$ and on the velocity relative to the air.

The goal is to integrate the equations of motion from take-off until the jumper intersects the landing hill, at $X(t_f)$.

= State vector

The state vector contains six variables:

$bold(x) = (X, Y, U_r, beta_r, phi, omega)$

where $X$ and $Y$ are position coordinates, $U_r$ is velocity relative to the air, and $beta_r$ is the jumper's flight-path angle relative to the horizontal. The angle $phi$ describes body orientation, and $omega$ is angular velocity. The aerodynamic lift $L$, drag $D$, and pitching moment $M$ are obtained from wind-tunnel measurements.

= Equations of motion

$
  cases(
    frac(d X, d t) = U_j cos(beta_j),
    frac(d Y, d t) = U_j sin(beta_j),
    frac(d U_r, d t) = -g sin(beta_r) - frac(D, m),
    frac(d beta_r, d t) = (frac(L, m) - g cos(beta_r)) / U_r,
    frac(d phi, d t) = omega,
    frac(d omega, d t) = frac(M, I)
  )
$

The main numerical task is to solve this nonlinear ordinary differential equation system using a Runge–Kutta method.

= Suggested tasks

- Implement the aerodynamic model.
- Solve the six coupled equations with the fourth-order Runge–Kutta method (RK4).
- Perform a time-step convergence study.
- Determine the landing point numerically.
- Study the effect of $theta$, $lambda$, and wind speed.
- Optionally, optimize the posture to maximize flight distance.

= References

- Seo, K., Murakami, M., & Yoshida, K. (2004). “Optimal flight technique for V-style ski jumping.” *Sports Engineering, 7*(2), 97–103. https://doi.org/10.1007/BF02915921
- Seo, K., Watanabe, I., & Murakami, M. (2004). “Aerodynamic force data for a V-style ski jumping flight.” *Sports Engineering, 7*(1), 31–39. https://doi.org/10.1007/BF02843971
