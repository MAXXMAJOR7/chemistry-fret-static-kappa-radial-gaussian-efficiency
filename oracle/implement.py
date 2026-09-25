"""Hidden oracle for Blackbox Task 50 (reference implementation, never shown to the solver).

Ensemble-averaged Foerster transfer efficiency (percent) for a donor-acceptor pair whose
orientations AND separations are frozen on the transfer time scale (static averaging regime).

Inputs (neutral names in the public schema):
    value_a  a      centre parameter of the radial-Gaussian D-A distance distribution, Angstrom
    value_b  sigma  width parameter of that distribution, Angstrom (b = 1 / (2 sigma^2)); 0 -> single distance
    value_c  Q_D    donor fluorescence quantum yield (no acceptor)
    value_d  J      spectral overlap integral, M^-1 cm^-1 nm^4

Pipeline:
    S1  Foerster prefactor from first principles:
            R0^6 = 9000 ln10 kappa^2 Q_D J / (128 pi^5 N_A n^4)          (Foerster 1948)
        with n = 1.4 (biomolecules in aqueous solution); written as R0^6 = kappa^2 X.
    S2  Static isotropic orientation average (TWIST 1). kappa^2 is frozen per molecule and
        distributed as p(k2) = ln(2+sqrt3) / (2 sqrt(3 k2))                        0 <= k2 <= 1
                              ln((2+sqrt3)/(sqrt k2 + sqrt(k2-1))) / (2 sqrt(3 k2))  1 <= k2 <= 4
        (Dale, Eisinger & Blumberg 1979; Loura 2012 eq 1). The efficiency, not kappa^2, is averaged:
            E_k(r) = int p(k2) k2 X / (k2 X + r^6) dk2.
        With u = |kappa| and c = r^6 / X, integration by parts over the log kink at u = 1 gives the
        smooth form  E_k(r) = (1/sqrt3) int_0^{arccosh 2} [cosh t - sqrt(c) atan(cosh t / sqrt(c))] dt.
    S3  Static distance average (TWIST 2) over the radial Gaussian
            P(r) = K r^2 exp(-b (r - a)^2),  int_0^inf P dr = 1        (Haas et al. 1975; Koren et al. eq 4)
    S4  output = 100 * <E>, rounded to 4 decimals.
"""
import math

# --- S1 constants -----------------------------------------------------------------------------
N_A = 6.02214076e23          # Avogadro constant, mol^-1 (exact, SI 2019)
FOERSTER_NUM = 9000.0        # 9000 ln10 / (128 pi^5 N_A n^4): Foerster (1948) numerator (1000 cm^3/L x 9)
FOERSTER_DEN = 128.0         # ... and denominator 128 pi^5
N_REFR = 1.4                 # refractive index assumed for biomolecules in aqueous solution
NM4_TO_CM4 = 1e-28           # J given in M^-1 cm^-1 nm^4; 1 nm^4 = 1e-28 cm^4 (exact)
CM6_TO_A6 = 1e48             # 1 cm^6 = 1e48 Angstrom^6 (exact)

# --- quadrature (fixed orders -> deterministic) ------------------------------------------------
N_GL_KAPPA = 48              # Gauss-Legendre nodes for the orientation integral (smooth integrand)
N_GL_R = 32                  # Gauss-Legendre nodes per panel for the distance integral
N_PANELS_R = 16              # panels spanning a +/- 10 sigma window (clipped at r = 0)
HALF_WIDTH_SIGMAS = 10.0     # radial Gaussian tail beyond 10 sigma < 1e-21 of the mass
N_NEWTON = 100               # Newton iterations for Legendre roots (converges in < 10)


def _gauss_legendre(n):
    """Nodes/weights on [-1, 1] by Newton iteration on P_n (Golub-Welsch equivalent, stdlib only)."""
    nodes, weights = [], []
    for i in range(1, n + 1):
        x = math.cos(math.pi * (i - 0.25) / (n + 0.5))
        for _ in range(N_NEWTON):
            p0, p1 = 1.0, x
            for k in range(2, n + 1):
                p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
            dp = n * (x * p1 - p0) / (x * x - 1)
            dx = p1 / dp
            x -= dx
            if abs(dx) < 1e-16:
                break
        p0, p1 = 1.0, x
        for k in range(2, n + 1):
            p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
        dp = n * (x * p1 - p0) / (x * x - 1)
        nodes.append(x)
        weights.append(2 / ((1 - x * x) * dp * dp))
    return nodes, weights


_GL_K = _gauss_legendre(N_GL_KAPPA)
_GL_R = _gauss_legendre(N_GL_R)
SQRT3 = math.sqrt(3)
T_MAX = math.acosh(2)        # = ln(2 + sqrt 3): upper limit of |kappa| = 2 (in-line collinear dipoles)


def foerster_x(q_d, j_nm):
    """S1: X = R0^6 / kappa^2 in Angstrom^6."""
    pref = FOERSTER_NUM * math.log(10) / (FOERSTER_DEN * math.pi ** 5 * N_A * N_REFR ** 4)
    return pref * q_d * j_nm * NM4_TO_CM4 * CM6_TO_A6


def kappa_averaged_efficiency(r, x):
    """S2: static isotropic kappa^2 average of k2 X / (k2 X + r^6)."""
    c = r ** 6 / x
    if c == 0:
        return 1.0
    sc = math.sqrt(c)
    nodes, weights = _GL_K
    half = 0.5 * T_MAX
    acc = 0.0
    for z, w in zip(nodes, weights):
        ch = math.cosh(half * (z + 1))
        acc += w * (ch - sc * math.atan(ch / sc))
    return half * acc / SQRT3


def distance_averaged_efficiency(a, sigma, x):
    """S3: average of S2 over P(r) = K r^2 exp(-(r - a)^2 / (2 sigma^2)) on r >= 0."""
    if sigma == 0:
        return kappa_averaged_efficiency(a, x)
    lo = max(0.0, a - HALF_WIDTH_SIGMAS * sigma)
    hi = a + HALF_WIDTH_SIGMAS * sigma
    step = (hi - lo) / N_PANELS_R
    nodes, weights = _GL_R
    num = den = 0.0
    for p in range(N_PANELS_R):
        mid = lo + (p + 0.5) * step
        for z, w in zip(nodes, weights):
            r = mid + 0.5 * step * z
            pr = w * r * r * math.exp(-((r - a) / sigma) ** 2 / 2)
            if pr > 0:
                num += pr * kappa_averaged_efficiency(r, x)
                den += pr
    return num / den


def compute(a, sigma, q_d, j_nm):
    if not (a > 0 and sigma >= 0 and 0 < q_d <= 1 and j_nm > 0):
        raise ValueError("input outside physical domain")
    x = foerster_x(q_d, j_nm)
    return round(100 * distance_averaged_efficiency(a, sigma, x), 4)


def oracle(inputs):
    return {"output": compute(float(inputs["value_a"]), float(inputs["value_b"]),
                              float(inputs["value_c"]), float(inputs["value_d"]))}


if __name__ == "__main__":
    import json
    import sys
    print(json.dumps(oracle(json.loads(sys.argv[1]))))
