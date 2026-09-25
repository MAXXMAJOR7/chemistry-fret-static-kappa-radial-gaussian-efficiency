import statistics, sys
sys.path.insert(0, "scripts")
from helpers import load_golden, variant
cases = load_golden()["cases"]
def ins(c): return tuple(c["input"][k] for k in ("value_a","value_b","value_c","value_d"))
for name, kw in [("Textbook: single r, kappa^2 = 2/3", dict(static_kappa=False, distribution=False)),
                 ("Dynamic kappa^2 = 2/3, full radial distribution (T1 off)", dict(static_kappa=False)),
                 ("Static kappa^2, single r = value_a (T2 off)", dict(distribution=False)),
                 ("Static kappa^2, plain Gaussian without r^2", dict(radial=False)),
                 ("Dynamic kappa^2, plain Gaussian", dict(static_kappa=False, radial=False))]:
    e = [abs(variant(*ins(c), **kw) - c["expected"]["output"]) for c in cases]
    print(f"| {name} | {sum(x<=0.001 for x in e)}/51 | {statistics.median(e):.3f} | {max(e):.2f} |")
