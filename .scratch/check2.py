import math, random, sys, time
sys.path.insert(0, "oracle")
import implement as o

def ref(a, s, q, j, panels=200):
    x = o.foerster_x(q, j)
    if s == 0: return 100*o.kappa_averaged_efficiency(a, x)
    lo, hi = max(0.0, a-12*s), a+12*s
    n = panels*40; h=(hi-lo)/n; num=den=0
    for i in range(n+1):
        r = lo+i*h; w = (1 if i in (0,n) else (4 if i%2 else 2))*h/3
        p = w*r*r*math.exp(-((r-a)/s)**2/2)
        num += p*o.kappa_averaged_efficiency(r,x); den += p
    return 100*num/den
def dyn(a, s, q, j):
    x = o.foerster_x(q, j)*2/3
    return 100/(1+a**6/x)
rng = random.Random(3); worst = 0
t0=time.time()
for _ in range(60):
    a = rng.uniform(20,90); s = rng.uniform(0,12); q = rng.uniform(0.05,1); j = 10**rng.uniform(13,16)
    v = o.compute(a,s,q,j); e = abs(v-ref(a,s,q,j)); worst=max(worst,e)
print("worst quad err", worst, "time/call", (time.time()-t0)/60)
for args in [(20,12,1,1e16),(20,12,0.05,1e13),(90,12,0.05,1e13),(50,0,0.5,1e15),(50,8,0.5,1e15)]:
    print(args, o.compute(*args), "dyn single", round(dyn(*args),4), "ref", round(ref(*args),6))
