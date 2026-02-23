=== Notes

Important equations:
$
  partial_t u_eta + u_eta dot nabla u_eta + nabla p_eta - v nabla^2 u_eta + (chi_Omega) / eta (u_eta - u_s) = f \
  chi_Omega = cases(
    1 #h(1em) "for" x in overline(Omega)_s,
    0 #h(1em) "for" x in Omega_f
  ) \
  partial_t omega_eta + u_eta dot nabla omega_eta - nu nabla^2 omega_eta + nabla times ((chi_Omega) / eta (u_eta - u_s)) = nabla times f #h(1em) \
  omega_eta = nabla times u_eta \
$

Discretization:
$
  omega_eta (x, t) = sum_(k in ZZ^2) hat(omega)_eta (k, t) e^(i k dot x) \
  hat(omega)_eta (k, t) = 1/(4 pi^2) integral_Omega omega_eta (x, t) e^(- i k dot x) dif Omega
$

Where $k = (k_x, k_y)$ and $x in [0, 2 pi]^2$. $k_x in [-N_x \/ 2, N_x \/ 2), k_y in [-N_y \/ 2, N_y \/ 2)$,
$
  u_eta (x, t) = - sum_(k in ZZ^2, k != 0) (i k^perp) / abs(k)^2 hat(omega)_eta (k, t) e^(i k dot x) + U_infinity #h(1em) (1)
$

where $k^perp = (- k_y, k_x)$; truncate with $2/3$ rule (zero highest frequency $1 \/ 3$ nodes)

Time discretization, omitting $eta$ subscript:
$
  partial_t omega - nu nabla^2 omega = N(omega), N(omega) = - u dot nabla omega - nabla times (chi_Omega / eta (u - u_s)) \
  partial_t hat(omega) + nu abs(k)^2 hat(omega) = hat(N) (omega) #h(1em) ("Fourier space") \
$

Simple multiplication by integrating factor to get:
$
  hat(omega) (k, t_(n + 1)) = e^(- nu Delta t_(n + 1) abs(k)^2) hat(omega) (k, t_n) + integral_(t_n)^(t_(n + 1)) e^(- nu (t_(n + 1) - s) abs(k)^2) hat(N) (hat(omega) (k, s)) dif s \
  hat(omega) (k, t_(n + 1)) = e^(- nu Delta t_(n + 1) abs(k)^2) (hat(omega) (k, t_n) + beta_10 hat(N)^n beta_11 e^(-nu Delta t_n abs(k)^2) hat(N)^(n - 1)) #h(1em) ("Discretized") \
  beta_10 = 1/2 (Delta t_(n + 1)) / (Delta t_n) (Delta t_(n + 1) + 2 Delta t_n), beta_11 = -1/2 (Delta t_(n + 1)^2) / (Delta t_n), Delta t_n = t_n - t_(n - 1)
$
For startup with no previous data we use a first-order scheme.. and for determining timestep:
$
  U_"max" = max_(x in"Grid") abs(u(x, t_n)) \
  Delta t_(n + 1) = C Delta x \/ U_"max", Delta x = min(L_x \/ N_x, L_y \/ N_y), Delta t_(n + 1) < eta
$

=== Implementation

Time discretization tells us how to propagate $omega$, for which we require $u$, which we get via (1)
- Start with initial $omega$
- Bring into Fourier space
- Propagate

Fluid-only simplification
$
  N(hat(omega)) = -u dot nabla hat(omega) \
  hat(omega) (k, t) = 1/(4 pi^2) integral_Omega omega (x, t) e^(- i k dot x) dif Omega \
  nabla hat(omega) = -1/(4 pi^2) integral_Omega omega (x, t) i k e^(- i k dot x) dif Omega \
$
