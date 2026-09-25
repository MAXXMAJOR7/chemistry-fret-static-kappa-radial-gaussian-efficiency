"""Automated GATE 2 checks. Exits non-zero on any failure."""
import ast
import json
import math
import random
import sys

from helpers import (GOLDEN_PATH, INPUT_KEYS, INPUT_RANGES, ORACLE_PATH, load_golden,
                     load_oracle, sample_inputs, variant)

# Every numeric literal allowed in oracle/implement.py, with its justification.
ALLOWED_CONSTANTS = {
    6.02214076e23: "Avogadro constant (exact, SI 2019)",
    9000.0: "Foerster (1948) R0 prefactor numerator 9000 ln10 (1000 cm^3 L^-1 x 9)",
    128.0: "Foerster (1948) R0 prefactor denominator 128 pi^5",
    5: "pi^5 in the Foerster denominator 128 pi^5",
    4: "n^4 refractive-index dependence (Foerster 1948)",
    6: "r^6 distance law (Foerster 1948)",
    1.4: "refractive index for biomolecules in aqueous solution (Lakowicz; Wojcik et al. 2018)",
    1e-28: "nm^4 -> cm^4 (exact)",
    1e48: "cm^6 -> Angstrom^6 (exact)",
    3: "sqrt(3) of the isotropic kappa^2 density (Dale-Eisinger-Blumberg 1979; Loura 2012 eq 1)",
    2: "acosh(2): |kappa| max = 2 (collinear dipoles); r^2 radial shell; Gaussian 1/2; GL weights",
    10.0: "log10 base (ln 10) and +/-10 sigma integration window (tail < 1e-21)",
    100: "percent output; Newton iteration cap",
    48: "Gauss-Legendre order, orientation integral (fixed -> deterministic)",
    32: "Gauss-Legendre order per panel, distance integral",
    16: "number of distance panels",
    0.25: "Legendre root initial guess cos(pi (i - 1/4)/(n + 1/2)) (Abramowitz & Stegun 22.16.6)",
    0.5: "Legendre root initial guess; panel mid-points; half-interval maps",
    1e-16: "Newton convergence threshold (machine epsilon)",
    1: "structural",
    0: "structural (domain checks, r >= 0 clip)",
}

failures = []


def check(name, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {name}" + (f" -- {detail}" if detail else ""))
    if not cond:
        failures.append(name)


def numeric_literals(path):
    tree = ast.parse(open(path, encoding="utf-8").read())
    vals = []
    for node in ast.walk(tree):
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub) and isinstance(node.operand, ast.Constant):
            if isinstance(node.operand.value, (int, float)) and not isinstance(node.operand.value, bool):
                vals.append(-node.operand.value)
                node.operand.value = None
        elif isinstance(node, ast.Constant) and isinstance(node.value, (int, float)) and not isinstance(node.value, bool):
            vals.append(node.value)
    return vals


def kappa_mc(r, x, n=200000, seed=11):
    """Independent Monte Carlo of the static isotropic kappa^2 average (random D and A unit vectors)."""
    rng = random.Random(seed)

    def unit():
        z = rng.uniform(-1, 1)
        p = rng.uniform(0, 2 * math.pi)
        s = math.sqrt(1 - z * z)
        return s * math.cos(p), s * math.sin(p), z

    acc = 0.0
    for _ in range(n):
        d, a = unit(), unit()
        k2 = (d[0] * a[0] + d[1] * a[1] + d[2] * a[2] - 3 * d[2] * a[2]) ** 2
        acc += k2 * x / (k2 * x + r ** 6)
    return acc / n


def main():
    o = load_oracle()
    rng = random.Random(7)
    pts = [sample_inputs(rng) for _ in range(120)]
    outs = [o.compute(*p) for p in pts]

    # Determinism
    check("deterministic", all(o.compute(*p) == v for p, v in zip(pts, outs)))

    # Output schema
    out = o.oracle({"value_a": 40.0, "value_b": 5.0, "value_c": 0.5, "value_d": 1e15})
    check("output key is 'output' and scalar float",
          list(out) == ["output"] and isinstance(out["output"], float))
    check("4-decimal precision", all(round(v, 4) == v for v in outs))
    check("output within [0, 100]", all(0 <= v <= 100 for v in outs))

    # No invented constants
    unknown = [v for v in numeric_literals(ORACLE_PATH) if v not in ALLOWED_CONSTANTS]
    check("every oracle numeric literal is justified", not unknown, f"unjustified: {unknown}" if unknown else "")

    # Literature anchors
    pref = o.foerster_x(1.0, 1.0) * o.N_REFR ** 4
    check("Foerster prefactor reproduces 8.79e-5 (Lakowicz eq 13.5)", abs(pref - 8.79e-5) < 1e-7, f"{pref:.5e}")
    r0 = 0.211 * (2 / 3 * o.N_REFR ** -4 * 0.5 * 1e15) ** (1 / 6)
    r0_or = (2 / 3 * o.foerster_x(0.5, 1e15)) ** (1 / 6)
    check("R0 matches Lakowicz 0.211 (kappa^2 n^-4 Q J)^(1/6) to its 3 s.f.", abs(r0 / r0_or - 1) < 2e-3, f"{r0_or:.3f} vs {r0:.3f} A")
    x = 1.0
    mean_k2 = o.kappa_averaged_efficiency(10.0, x) * 1e6  # far limit: <E> -> <kappa^2> x / r^6 (1 + O(1e-6))
    check("static average far limit reproduces <kappa^2> = 2/3", abs(mean_k2 - 2 / 3) < 2e-6, f"{mean_k2:.8f}")
    check("static average near limit -> 1", abs(o.kappa_averaged_efficiency(1e-4, x) - 1) < 1e-5)
    mc = [(r, kappa_mc(r, x), o.kappa_averaged_efficiency(r, x)) for r in (0.7, 1.0, 1.3)]
    check("kappa^2 average agrees with Monte Carlo over random dipoles (3 sigma_MC ~ 3e-3)",
          all(abs(m - q) < 3e-3 for _, m, q in mc), "; ".join(f"r={r}: MC {m:.4f} vs {q:.4f}" for r, m, q in mc))
    # 100 - E ~ r^3 at short range (static) vs r^6 (dynamic)
    s1, s2 = (1 - o.kappa_averaged_efficiency(r, x) for r in (0.01, 0.02))
    slope = math.log(s2 / s1) / math.log(2)
    check("short-range tail exponent is 3 (static kappa^2 signature)", abs(slope - 3) < 0.01, f"{slope:.4f}")

    # Oracle agrees with the helper re-implementation
    check("helper variant(all twists) == oracle", all(variant(*p) == v for p, v in zip(pts, outs)))

    # Twist integrity
    t1 = max(abs(v - variant(*p, static_kappa=False)) for p, v in zip(pts, outs))
    t2 = max(abs(v - variant(*p, distribution=False)) for p, v in zip(pts, outs))
    t2b = max(abs(v - variant(*p, radial=False)) for p, v in zip(pts, outs))
    check("T1 static-kappa twist changes output (> 1000x tolerance)", t1 > 1.0, f"max diff {t1:.3f}")
    check("T2 distance-distribution twist changes output (> 1000x tolerance)", t2 > 1.0 and t2b > 1.0,
          f"vs single r {t2:.3f}, vs plain Gaussian {t2b:.3f}")
    check("T2 exactly neutral at value_b = 0",
          all(o.compute(p[0], 0.0, p[2], p[3]) == variant(p[0], 0.0, p[2], p[3], distribution=False) for p in pts))
    far = [o.compute(80.0, 0.0, 0.05, 1e13), variant(80.0, 0.0, 0.05, 1e13, static_kappa=False)]
    check("T1 neutral (within tolerance) in the far regime", abs(far[0] - far[1]) <= 0.001, f"{far}")
    both = o.compute(40.0, 8.0, 0.5, 1e15)
    only1 = variant(40.0, 8.0, 0.5, 1e15, distribution=False)
    only2 = variant(40.0, 8.0, 0.5, 1e15, static_kappa=False)
    none = variant(40.0, 8.0, 0.5, 1e15, static_kappa=False, distribution=False)
    check("twists interact (not additive)", abs((both - only1) - (only2 - none)) > 0.1,
          f"T2 effect with T1 {both - only1:.4f} vs without {only2 - none:.4f}")

    # Golden data
    g = load_golden()
    cats = [c["category"] for c in g["cases"]]
    check(">=2 discriminating edge cases", cats.count("discriminating") >= 2)
    check(">=2 control edge cases", cats.count("control") >= 2)
    check(">=1 boundary edge case", cats.count("boundary") >= 1)
    check("golden inputs neutral and in range", all(
        sorted(c["input"]) == sorted(INPUT_KEYS)
        and all(INPUT_RANGES[k][0] <= c["input"][k] <= INPUT_RANGES[k][1] for k in INPUT_KEYS)
        for c in g["cases"]))
    check("golden expected outputs match oracle",
          all(o.oracle(c["input"]) == c["expected"] for c in g["cases"]), GOLDEN_PATH)

    print(json.dumps({"failures": failures}))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
