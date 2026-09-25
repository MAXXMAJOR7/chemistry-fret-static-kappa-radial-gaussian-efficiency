import random, sys
sys.path.insert(0, "scripts")
from helpers import load_oracle, sample_inputs, variant
o = load_oracle(); rng = random.Random(50)
pts = [sample_inputs(rng) for _ in range(150)]
outs = [o.compute(*p) for p in pts]
print("frac in (1,99):", sum(1 < v < 99 for v in outs)/len(outs), "min", min(outs), "max", max(outs))
import statistics
for name, kw in [("dynamic k2", dict(static_kappa=False)), ("no r^2", dict(radial=False)),
                 ("single r", dict(distribution=False)), ("textbook", dict(static_kappa=False, distribution=False))]:
    errs = [abs(variant(*p, **kw) - v) for p, v in zip(pts, outs)]
    print(f"{name:12s} pass {sum(e<=0.001 for e in errs)}/{len(errs)} median {statistics.median(errs):.4f} max {max(errs):.3f}")
