"""Shared utilities for task scripts and grader."""
import importlib.util
import json
import math
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GOLDEN_PATH = os.path.join(ROOT, "golden", "test_data.json")
ORACLE_PATH = os.path.join(ROOT, "oracle", "implement.py")

INPUT_KEYS = ("value_a", "value_b", "value_c", "value_d")
INPUT_RANGES = {
    "value_a": (15.0, 80.0),
    "value_b": (0.0, 12.0),
    "value_c": (0.05, 1.0),
    "value_d": (1e13, 1e16),
}
TOLERANCE = 0.001  # absolute, output units (percentage points; output spans 0 .. 100)


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def load_oracle():
    return load_module(ORACLE_PATH, "oracle_impl")


def load_golden(path=GOLDEN_PATH):
    with open(path, "r", encoding="utf-8") as fh:
        return json.load(fh)


def sample_inputs(rng):
    """Random point in the input domain (value_d log-uniform)."""
    lo, hi = (math.log10(v) for v in INPUT_RANGES["value_d"])
    return (rng.uniform(*INPUT_RANGES["value_a"]), rng.uniform(*INPUT_RANGES["value_b"]),
            rng.uniform(*INPUT_RANGES["value_c"]), 10 ** rng.uniform(lo, hi))


def variant(a, sigma, q_d, j_nm, static_kappa=True, radial=True, distribution=True):
    """Oracle with individual twists frozen, for twist-integrity checks and baselines.

    static_kappa=False -> dynamic isotropic average: kappa^2 = 2/3 inside the Foerster kernel
                          (the textbook R0), instead of averaging E over the static p(kappa^2)
    radial=False       -> plain Gaussian P(r) ~ exp(-(r-a)^2 / 2 sigma^2) without the r^2 shell factor
    distribution=False -> single distance r = value_a (value_b ignored)
    """
    o = load_oracle()
    x = o.foerster_x(q_d, j_nm)

    def kernel(r):
        if static_kappa:
            return o.kappa_averaged_efficiency(r, x)
        return 1 / (1 + r ** 6 / (x * 2 / 3))

    if not distribution or sigma == 0:
        return round(100 * kernel(a), 4)
    lo = max(0.0, a - o.HALF_WIDTH_SIGMAS * sigma)
    hi = a + o.HALF_WIDTH_SIGMAS * sigma
    step = (hi - lo) / o.N_PANELS_R
    nodes, weights = o._GL_R
    num = den = 0.0
    for p in range(o.N_PANELS_R):
        mid = lo + (p + 0.5) * step
        for z, w in zip(nodes, weights):
            r = mid + 0.5 * step * z
            pr = w * (r * r if radial else 1.0) * math.exp(-((r - a) / sigma) ** 2 / 2)
            if pr > 0:
                num += pr * kernel(r)
                den += pr
    return round(100 * num / den, 4)
