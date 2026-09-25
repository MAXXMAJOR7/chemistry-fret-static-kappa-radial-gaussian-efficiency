import sys
sys.path.insert(0, "scripts")
from helpers import load_oracle, variant
o = load_oracle()
cands = [
 (80, 0, 0.05, 1e13), (80, 0, 0.3, 1e14), (70, 0, 1.0, 1e16), (60, 0, 0.5, 1e14),
 (40, 0, 0.5, 1e15), (25, 0, 1.0, 1e16), (30, 0, 0.5, 1e15),
 (40, 8, 0.5, 1e15), (15, 12, 0.3, 1e14), (20, 12, 1.0, 1e16), (80, 12, 1.0, 1e16),
 (15, 0, 1.0, 1e16), (80, 0, 0.05, 1e13), (40, 0.5, 0.5, 1e15), (15,12,0.05,1e13),
]
x = lambda q,j: (o.foerster_x(q,j)*2/3)**(1/6)
for c in cands:
    print(c, "R0=%.1f"%x(c[2],c[3]), o.compute(*c), "dyn", variant(*c, static_kappa=False), "noR2", variant(*c, radial=False), "single", variant(*c, distribution=False))
