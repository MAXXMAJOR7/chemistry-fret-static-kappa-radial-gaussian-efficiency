import math, random, sys
sys.path.insert(0, "oracle")
import implement as o

# R0 prefactor vs Lakowicz 8.79e-5 (n^-4 included separately)
print("prefactor*n^4 =", o.foerster_x(1,1)*o.N_REFR**4)
x = o.foerster_x(1.0, 1e16); print("R0(2/3) max", (2/3*x)**(1/6))
x = o.foerster_x(0.05, 1e13); print("R0(2/3) min", (2/3*x)**(1/6))

# Monte Carlo check of kappa^2 average: random unit vectors for D and A, r along z
rng = random.Random(1)
def unit():
    z = rng.uniform(-1,1); p = rng.uniform(0, 2*math.pi); s = math.sqrt(1-z*z)
    return (s*math.cos(p), s*math.sin(p), z)
ks = []
for _ in range(400000):
    d, a = unit(), unit()
    k = d[0]*a[0]+d[1]*a[1]+d[2]*a[2] - 3*d[2]*a[2]
    ks.append(k*k)
print("MC <k2>", sum(ks)/len(ks))
X = 1.0
for r in (0.5, 0.8, 1.0, 1.2, 1.6):
    mc = sum(k*X/(k*X + r**6) for k in ks)/len(ks)
    dyn = 1/(1 + r**6/(2/3*X))
    print(r, "oracle", o.kappa_averaged_efficiency(r, X), "MC", mc, "dynamic", dyn)

# quadrature convergence: brute force midpoint in u with density
def brute(r, X, n=400000):
    c = r**6/X; A = math.log(2+math.sqrt(3)); s = 0; h = 2/n
    for i in range(n):
        u = (i+0.5)*h
        q = A if u < 1 else A - math.acosh(u)
        s += q*u*u/(u*u+c)*h
    return s/math.sqrt(3)
for r in (0.2, 1.0, 2.0):
    print("brute", r, brute(r,1.0), o.kappa_averaged_efficiency(r,1.0))
