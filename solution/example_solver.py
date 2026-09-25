"""Example solver: probe the black box, then fit a hypothesised model.

Strategy (what an expert would do after exploratory probing):
  1. Output lives in 0..100 and falls sigmoidally with value_a -> a percentage.
     value_c and value_d enter only through their product (checked by swapping factors
     at fixed product) and the curve shifts in value_a by (c*d)^(1/6): a Foerster
     efficiency with R0^6 proportional to Q_D * J, value_a a distance in Angstrom.
  2. With value_b = 0, the far tail goes as r^-6 (as expected) but 100 - output at short
     range goes as r^+3, not r^+6.  An r^3 tail is the signature of a kappa^2 density that
     diverges as 1/sqrt(kappa^2) at 0: static isotropic orientations, efficiency averaged
     over p(kappa^2) (Dale-Eisinger-Blumberg).  The absolute scale fixes the Foerster
     prefactor; it matches 9000 ln10 / (128 pi^5 N_A n^4) with n = 1.4.
  3. value_b broadens the response symmetrically in log-output only after accounting
     for an r^2 Jacobian: compare plain Gaussian vs radial (r^2-weighted) Gaussian in r
     centred at value_a with standard deviation value_b; the radial form wins.
  Model: fitted prefactor, static kappa^2 average, radial-Gaussian distance average.

The solver uses its own numerics (Simpson rules), independent of the oracle.
This solver is illustrative, not guaranteed optimal. Standard library only.
"""
import math

PROBE_BUDGET = 200

_params = None
_A = math.log(2 + math.sqrt(3))


def _e_kappa(r, x, static=True):
    """Efficiency at distance r for R0^6 = kappa^2 * x, averaged over static isotropic kappa^2."""
    if not static:
        return 1 / (1 + r ** 6 / (x * 2 / 3))
    c = r ** 6 / x
    if c == 0:
        return 1.0
    # integrate over |kappa| = u with density A/sqrt3 on [0,1] and (A - acosh u)/sqrt3 on [1,2];
    # the first piece is closed-form, the second uses u = cosh t (Simpson, 200 intervals)
    sc = math.sqrt(c)
    part1 = _A * (1 - sc * math.atan(1 / sc))
    n, tmax = 200, _A
    h = tmax / n
    s = 0.0
    for i in range(n + 1):
        t = i * h
        u = math.cosh(t)
        f = (_A - t) * u * u / (u * u + c) * math.sinh(t)
        s += f * (1 if i in (0, n) else (4 if i % 2 else 2))
    return (part1 + s * h / 3) / math.sqrt(3)


def _model(a, sigma, q, j, pref, static=True, radial=True):
    x = pref * q * j
    if sigma == 0:
        return 100 * _e_kappa(a, x, static)
    lo, hi, n = max(0.0, a - 9 * sigma), a + 9 * sigma, 400
    h = (hi - lo) / n
    num = den = 0.0
    for i in range(n + 1):
        r = lo + i * h
        w = (1 if i in (0, n) else (4 if i % 2 else 2)) * math.exp(-((r - a) / sigma) ** 2 / 2)
        if radial:
            w *= r * r
        if w > 0:
            num += w * _e_kappa(r, x, static)
            den += w
    return 100 * num / den


def fit(query):
    global _params

    def q(a, b, c, d):
        return query({"value_a": a, "value_b": b, "value_c": c, "value_d": d})["output"]

    # Step 1: scale. Fit the prefactor on a mid-efficiency probe for both kappa hypotheses.
    fits = {}
    target = q(40.0, 0.0, 0.5, 1e15)
    for static in (True, False):
        lo, hi = -8.0, -2.0            # log10 prefactor (Angstrom^6 per unit Q*J)
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if _model(40.0, 0.0, 0.5, 1e15, 10 ** mid, static) < target:
                lo = mid
            else:
                hi = mid
        fits[static] = 10 ** (0.5 * (lo + hi))
    # Step 2: short-range tail decides static vs dynamic averaging.
    y = q(20.0, 0.0, 1.0, 1e16)
    static = abs(_model(20.0, 0.0, 1.0, 1e16, fits[True], True) - y) < \
        abs(_model(20.0, 0.0, 1.0, 1e16, fits[False], False) - y)
    pref = fits[static]
    # Step 3: shape of the distance distribution (strongest contrast at small a, wide sigma).
    y = q(15.0, 12.0, 0.3, 1e14)
    radial = abs(_model(15.0, 12.0, 0.3, 1e14, pref, static, True) - y) < \
        abs(_model(15.0, 12.0, 0.3, 1e14, pref, static, False) - y)
    _params = (pref, static, radial)


def solve(inputs):
    if _params is None:
        raise RuntimeError("call fit(query) first")
    pref, static, radial = _params
    a, b, c, d = (float(inputs[k]) for k in ("value_a", "value_b", "value_c", "value_d"))
    return {"output": round(_model(a, b, c, d, pref, static, radial), 4)}
