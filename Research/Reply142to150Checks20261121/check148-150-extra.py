# Contract-free additions (2026-11-21) to check148-150-quick.py. NOT claims of >>148-150 unless marked "caveat".
# A. >>150's two caveats (root 17 has other lower solutions; vertices 31 and 5591 succeed outside the table).
# B. p0 = 42698041: which of the 32 window vertices succeed (all admissible sources), and how the window sits in a closed set.
# C. 288 inputs: every Type I witness written as (k, a, h): 4e + 1 = k c, N_k = (k p + 1)/4, x = b N_k / h, e = a N_k / h,
#    a + h = k b, gcd(a, h) = 1, a h | N_k, h c = a p + b.  Distribution of k, split lower / upper.
# D. 288 inputs: window U = non-residue primes <= floor(cbrt p); closure of l inside U; does it contain a successful vertex?
# Units: source = one admissible (q, t); "succeeds" = some d | x^2 with c | d + p x (Type I) or c | d + x (Type II).
import sys, time, math
from collections import Counter
import sympy
sys.stdout.reconfigure(encoding='utf-8')
t0=time.time()
HARD={1,121,169,289,361,529}
P17=[30241,61681,106681,111721,148201,153409,163249,206641,216841,280009,288649,307129,309241,325081,347041,352081,374929]
def A(q): return q if q%4==3 else 3*q
def nr(q,p): return pow(q,(p-1)//2,p)==p-1
def dsq(x,f=None):
    out=[1]
    for r,e in (f or sympy.factorint(x)).items(): out=[a*r**k for a in out for k in range(2*e+1)]
    return out
def least_nr(p):
    l=2
    while not nr(l,p): l=sympy.nextprime(l)
    return l
def icbrt(n):
    r=round(n**(1/3))
    while r**3>n: r-=1
    while (r+1)**3<=n: r+=1
    return r
def vertex(p,q,lower_only=False):
    """returns (number of sources, list of successful (t, branch, d), set of non-residue prime factors of the x's)"""
    At=A(q); t=1; n=0; succ=[]; kids=set()
    while At*t*t<3*p:
        c=At*t*t
        if lower_only and 3*c>p: break
        x=(p+c)//4; f=sympy.factorint(x); n+=1
        kids|={r for r in f if nr(r,p)}
        for d in dsq(x,f):
            if (d+p*x)%c==0: succ.append((t,'I',d)); break
        for d in dsq(x,f):
            if (d+x)%c==0: succ.append((t,'II',d)); break
        t+=2
    return n,succ,kids

print('--- A. caveats of >>150')
for p in (289249,13217743009):
    n,succ,_=vertex(p,17,lower_only=True)
    print('   p=%d: root 17 lower sources %d, successful (t, branch, d) [first 4 of %d]: %s'%(p,n,len(succ),succ[:4]))
p=13217743009
table={31:(5,19),5591:(1,)}
for q in (31,5591):
    n,succ,_=vertex(p,q,lower_only=True)
    out=[s for s in succ if s[0] not in table[q]]
    print('   p=%d: vertex %d lower sources %d, successful outside the table rows [first 3 of %d]: %s'%(p,q,n,len(out),out[:3]))

print('--- B. p0 = 42698041, window = non-residue primes <= 349')
p0=42698041
U=[q for q in sympy.primerange(3,350) if nr(q,p0)]
S=[41,67,83,131,151,157,181,269,313,331]
res={}
for q in U:
    n,succ,kids=vertex(p0,q); res[q]=(n,succ,kids)
print('   sources per vertex (41):',res[41][0],'; successes at 41:',len(res[41][1]))
ok=[q for q in U if res[q][1]]
print('   window vertices with a successful source: %d of %d: %s'%(len(ok),len(U),ok))
for q in ok: print('      %d: %d sources, first success %s'%(q,res[q][0],res[q][1][0]))
print('   successes among S:',[q for q in S if res[q][1]])
out_ok=[q for q in U if q not in S and all(sympy.jacobi_symbol((-A(s))%q,q)==1 for s in S)]
print('   window vertices not in S that no member of S reaches (all ten symbols +1):',out_ok)
dist={41:0}; frontier=[41]
while frontier:
    nxt=[]
    for q in frontier:
        for r in res[q][2]:
            if r<=349 and r not in dist: dist[r]=dist[q]+1; nxt.append(r)
    frontier=nxt
print('   distance from the root inside the window:',sorted(Counter(dist.values()).items()),'; unreached:',[q for q in U if q not in dist])
print('   successful window vertices and their distance:',[(q,dist.get(q)) for q in ok])
print('   %.1f s'%(time.time()-t0))

print('--- C/D. 288 inputs')
inputs=sorted(set(p for p in sympy.primerange(1009,100000) if p%840 in HARD)|{2521,12889}|set(P17))
kall=Counter(); kmin=Counter(); nwit=Counter(); bad=0
kfam=Counter()
win=Counter(); winrows=[]
for p in inputs:
    l=least_nr(p); W=icbrt(p)
    info={}   # for window vertices: (kids in window, success)
    k3=any(r%3==2 for r in sympy.factorint((3*p+1)//4))
    kfam['inputs']+=1; kfam['N_3 has a prime factor = 2 mod 3']+=k3
    anyk3=False
    for q in sympy.primerange(3,p):
        if not nr(q,p): continue
        At=A(q); t=1; vs=False; vk=set()
        while At*t*t<3*p:
            c=At*t*t; x=(p+c)//4; f=sympy.factorint(x); side='upper' if 3*c>p else 'lower'
            ds=dsq(x,f)
            t1=[d for d in ds if (d+p*x)%c==0]
            t2=any((d+x)%c==0 for d in ds)
            if t1 or t2: vs=True
            if q<=W: vk|={r for r in f if r<=W and nr(r,p)}
            if t1:
                ks=[]
                for d in t1:
                    e=x*x//d; k=(4*e+1)//c
                    g=math.gcd(x,e); a=e//g; b=x//g; h=k*b-a; Nk=(k*p+1)//4
                    if not ((4*e+1)%c==0 and k%4==3 and (k*p+1)%4==0 and h>0 and math.gcd(a,h)==1 and Nk%(a*h)==0 and h*c==a*p+b and x*h==b*Nk and e*h==a*Nk): bad+=1
                    ks.append(k); kall[(side,k)]+=1; nwit[side]+=1
                    if k==3: anyk3=True
                kmin[(side,min(ks))]+=1
            t+=2
        if q<=W: info[q]=(vk,vs)
    kfam['inputs with a Type I witness of multiplier k = 3 at some square source']+=anyk3
    kfam['... of which N_3 has a prime factor = 2 mod 3']+=(anyk3 and k3)
    # window closure
    if l>W:
        win['root above the window (l > cbrt p)']+=1; winrows.append((p,l,W,'root outside')); continue
    seen={l}; fr=[l]
    while fr:
        nx=[]
        for q in fr:
            for r in info[q][0]:
                if r not in seen: seen.add(r); nx.append(r)
        fr=nx
    has=any(info[q][1] for q in seen)
    win['root in the window']+=1
    win['root in the window, closure inside the window has a successful vertex']+=has
    if not has: winrows.append((p,l,W,len(info),len(seen)))
print('   C: Type I witnesses checked against the (k,a,h) form: lower %d, upper %d; identities failing: %d'%(nwit['lower'],nwit['upper'],bad))
for side in ('lower','upper'):
    tot=sum(v for (s,k),v in kmin.items() if s==side)
    top=sorted(((k,v) for (s,k),v in kmin.items() if s==side),key=lambda z:-z[1])[:10]
    print('   C: %s sources with a Type I witness: %d; smallest multiplier k per source, most frequent: %s'%(side,tot,top))
    print('      k = 3: %d ; k <= 11: %d ; k <= 39: %d ; largest smallest-k: %d'%(kmin[(side,3)],sum(v for (s,k),v in kmin.items() if s==side and k<=11),sum(v for (s,k),v in kmin.items() if s==side and k<=39),max(k for (s,k) in kmin if s==side)))
for k in kfam: print('   C: %-78s %d'%(k,kfam[k]))
for k in win: print('   D: %-78s %d'%(k,win[k]))
print('   D: inputs with the root in the window but no successful vertex in the window closure [first 12 of %d] (p, l, W, window vertices, closure size):'%sum(1 for r in winrows if r[3]!='root outside'))
for r in [r for r in winrows if r[3]!='root outside'][:12]: print('      ',r)
print('%.1f s'%(time.time()-t0))
