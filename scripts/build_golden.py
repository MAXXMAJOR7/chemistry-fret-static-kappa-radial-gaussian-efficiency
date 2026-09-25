"""Regenerate golden/test_data.json from the oracle (deterministic, seeded)."""
import json
import math
import random

from helpers import GOLDEN_PATH, load_oracle

# (id, category, label, value_a, value_b, value_c, value_d)
EDGE_CASES = [
    ("edge_control_1", "control", "Single distance, r = 2 R0: weak-transfer regime, orientation twist nearly neutral (dynamic 2/3 gives 1.6076)",
     60.0, 0.0, 0.5, 1e14),
    ("edge_control_2", "control", "Single distance, r = 2.9 R0: far regime, E ~ <kappa^2> X / r^6, twist neutral to 4e-4",
     80.0, 0.0, 0.3, 1e14),
    ("edge_control_3", "control", "Single distance, largest r / R0 (5.7): both twists neutral, textbook r^-6 law",
     80.0, 0.0, 0.05, 1e13),
    ("edge_discrim_1", "discriminating", "Single distance, r = 0.9 R0: static kappa^2 average gives 48.6 vs textbook 65.0",
     40.0, 0.0, 0.5, 1e15),
    ("edge_discrim_2", "discriminating", "Single distance, r = 0.68 R0: static average 73.6 vs textbook 91.3",
     30.0, 0.0, 0.5, 1e15),
    ("edge_discrim_3", "discriminating", "sigma = 8 A at r ~ R0: r^2-weighted Gaussian LOWERS E (42.4) while a plain Gaussian raises it (49.3)",
     40.0, 8.0, 0.5, 1e15),
    ("edge_discrim_4", "discriminating", "Wide distribution centred at 1.1 R0: both twists fire together",
     80.0, 12.0, 1.0, 1e16),
    ("edge_discrim_5", "discriminating", "Small centre, wide sigma: r^2 Jacobian empties short distances (44.2 vs plain Gaussian 71.6, single r 85.6)",
     15.0, 12.0, 0.3, 1e14),
    ("edge_boundary_1", "boundary", "Shortest distance, largest R0 (r = 0.2 R0): 100 - E ~ r^3 tail from kappa^2 -> 0 orientations (99.16 vs textbook 99.99)",
     15.0, 0.0, 1.0, 1e16),
    ("edge_boundary_2", "boundary", "Smallest R0 with the widest distribution at the lowest centre: P(r) clipped at r = 0",
     15.0, 12.0, 0.05, 1e13),
    ("edge_boundary_3", "boundary", "sigma -> 0 limit (0.5 A): continuous with the single-distance value 48.6019",
     40.0, 0.5, 0.5, 1e15),
]

N_RANDOM = 40
SEED = 50


def main():
    o = load_oracle()
    cases = []
    for cid, cat, label, a, b, c, d in EDGE_CASES:
        inp = {"value_a": a, "value_b": b, "value_c": c, "value_d": d}
        cases.append({"id": cid, "category": cat, "label": label, "input": inp, "expected": o.oracle(inp)})
    rng = random.Random(SEED)
    for i in range(N_RANDOM):
        inp = {
            "value_a": round(rng.uniform(15.0, 80.0), 3),
            "value_b": round(rng.uniform(0.0, 12.0), 3),
            "value_c": round(rng.uniform(0.05, 1.0), 4),
            "value_d": float(f"{10 ** rng.uniform(13.0, 16.0):.4e}"),
        }
        cases.append({"id": f"rand_{i:02d}", "category": "random", "label": "uniform value_a..value_c, log-uniform value_d",
                      "input": inp, "expected": o.oracle(inp)})
    doc = {
        "schema": {
            "input": {
                "value_a": "float [15, 80]",
                "value_b": "float [0, 12]",
                "value_c": "float [0.05, 1]",
                "value_d": "float [1e13, 1e16]",
            },
            "output": {"output": "float, 4 decimals"},
        },
        "seed": SEED,
        "cases": cases,
    }
    assert all(math.isfinite(c["expected"]["output"]) for c in cases)
    with open(GOLDEN_PATH, "w", encoding="utf-8") as fh:
        json.dump(doc, fh, indent=2)
        fh.write("\n")
    print(f"wrote {len(cases)} cases to {GOLDEN_PATH}")


if __name__ == "__main__":
    main()
