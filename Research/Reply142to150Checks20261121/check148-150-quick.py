# Contract-free quick check (2026-11-21) of >>148, >>149, >>150 (no PR material accompanies these posts; PR #4 head is still 20b7ac7).
# Written after reading the three posts. Everything is recomputed from the primes.
# Units: vertex = non-residue prime q < p; source = admissible odd t with c = A(q) t^2 < 3p, A(q) = q (q = 3 mod 4) or 3q;
# x = (p + c)/4; lower = 3c < p; Type I witness: d | x^2, c | d + p x  (equivalently e = x^2/d with c | 4e + 1);
# Type II witness: d | x^2, c | d + x (equivalently c | 4d + p).
import sys, time, math, random
from fractions import Fraction
from collections import Counter
import sympy
sys.stdout.reconfigure(encoding='utf-8')
t0=time.time()
HARD={1,121,169,289,361,529}
nok=0; ndiff=0
def chk(name,got,exp):
    global nok,ndiff
    if got==exp: nok+=1; print('OK  ',name,got if len(str(got))<200 else '(long)')
    else: ndiff+=1; print('DIFF',name,'got',got,'expected',exp)
def A(q): return q if q%4==3 else 3*q
def leg(a,p): return sympy.jacobi_symbol(a%p,p)
def nr(q,p): return pow(q,(p-1)//2,p)==p-1
def dsq(x):
    out=[1]
    for r,e in sympy.factorint(x).items(): out=[a*r**k for a in out for k in range(2*e+1)]
    return out
def wit(p,c,x):
    ds=dsq(x)
    return sorted(d for d in ds if (d+p*x)%c==0),sorted(d for d in ds if (d+x)%c==0)
def es(p,c,x,d,typ):
    # exact three-unit-fraction identity from a factor pair
    N=p*x; t=d if typ=='I' else p*d
    assert (N+t)%c==0 and (N*N)%t==0 and (N+N*N//t)%c==0
    y=(N+t)//c; z=(N+N*N//t)//c
    assert Fraction(1,x)+Fraction(1,y)+Fraction(1,z)==Fraction(4,p)
    return (x,y,z)
def least_nr(p):
    l=2
    while not nr(l,p): l=sympy.nextprime(l)
    return l
def src(p,q,t):
    c=A(q)*t*t; assert t%2==1 and c<3*p and (p+c)%4==0
    return c,(p+c)//4
def edge(p,q,t,r):
    c,x=src(p,q,t); return x%r==0 and sympy.isprime(r) and nr(r,p)

# ================= >>148 =================
print('--- >>148')
# Thue step, concretely: some 1 <= |a|,|b| <= M with a*l = b (mod p)
def thue(p,l):
    M=math.isqrt(p)
    for a in range(1,M+1):
        b=(a*l)%p
        if b>p//2: b-=p
        if b!=0 and abs(b)<=M: return a,b,M
for p in (40609,2521,42698041):
    l=least_nr(p); a,b,M=thue(p,l)
    print('     p=%d l=%d M=%d: a=%d, b=%d, a*l-b = %d*p ; l<=M: %s'%(p,l,M,a,b,(a*l-b)//p,l<=M))
hard7=[p for p in sympy.primerange(1000,10**7) if p%840 in HARD]
chk('148: hard primes below 10^7 with least non-residue prime l > floor(sqrt p)',sum(least_nr(p)>math.isqrt(p) for p in hard7),0)
# boundary example
p=40609
chk('148: p=40609: 176*3461 - 15p, r = 3461 > p/12',(176*3461-15*p,12*3461>p),(1,True))
# p0
p0=42698041
chk('148: p0 prime, class, least non-residue prime',(sympy.isprime(p0),p0%840 in HARD,least_nr(p0)),(True,True,41))
U=[q for q in sympy.primerange(3,350) if nr(q,p0)]
chk('148: non-residue primes <= 349',len(U),32)
kids41=set(); t=1
while A(41)*t*t<3*p0:
    c,x=src(p0,41,t); kids41|={r for r in sympy.factorint(x) if nr(r,p0)}; t+=2
S=sorted({41}|{r for r in kids41 if r<=349})
chk('148: S = root and its direct children that lie in U',S,[41,67,83,131,151,157,181,269,313,331])
sq=0; allplus=[]
for b in U:
    s=sum(leg(-A(q),b) for q in S); sq+=s*s
    if all(leg(-A(q),b)==1 for q in S): allplus.append(b)
chk('148: sum over the 32 vertices of (sum of the 10 symbols)^2',sq,326)
print('     148: vertices b of U with all ten symbols +1 (the only ones that may lie outside C):',allplus)
# the criterion: for b in U, q in S, b != q: (-A(q)/b) = -1  <=>  some admissible odd t has b | x(q,t)
bad=0
for q in S:
    At=A(q); T=1; hitset=set()
    while At*T*T<3*p0:
        x=(p0+At*T*T)//4
        for b in U:
            if x%b==0: hitset.add(b)
        T+=2
    for b in U:
        if b!=q and ((leg(-At,b)==-1)!=(b in hitset)): bad+=1
chk('148: symbol criterion vs direct enumeration of every admissible t (10 x 31 pairs): disagreements',bad,0)
chk('148: A(q) t^2 < 3 p0 for q <= 349, t <= 347',3*349*347*347<3*p0,True)
c,x=src(p0,43,1)
w1,w2=wit(p0,c,x)
chk('148: vertex 43, t=1: c, x, 143 | x, x + 143 = 43 * 248248, lower source',(c,x,x%143,x+143,3*c<p0),(43,10674521,0,43*248248,True))
chk('148: vertex 43: 143 is a Type II witness; exact identity',(143 in w2,es(p0,c,x,143,'II')[0]),(True,x))
into43=[q for q in U if q!=43 and leg(-A(q),43)==-1]
chk('148: 61, 89, 157, 241, 277, 313 have an edge to 43',all(q in into43 for q in (61,89,157,241,277,313)),True)
print('     148: all vertices of U with an edge to 43:',into43,'; of these in S:',[q for q in into43 if q in S])
# explicit t for 157 -> 43 and 41 -> 157
for q,r in ((41,157),(157,43),(41,313),(313,43)):
    T=1; found=None
    while A(q)*T*T<3*p0:
        if ((p0+A(q)*T*T)//4)%r==0: found=T; break
        T+=2
    print('     148: smallest t with %d | x(%d, t): %s'%(r,q,found))
chk('148: 2521: non-residue primes <= floor(cbrt 2521) = 13',[q for q in sympy.primerange(3,14) if nr(q,2521)],[11])

# ================= >>149 =================
print('--- >>149')
chk('149: 4*461 + 1 = 15 * 123, 461 prime > 349',(4*461+1,sympy.isprime(461)),(15*123,True))
mods=[840]+list(sympy.primerange(11,350))+[461]
rems=[p0%840]+[p0%m for m in sympy.primerange(11,350)]+[(-123)%461]
from sympy.ntheory.modular import crt
base,M=crt(mods,rems); base=int(base); M=int(M)
random.seed(149); P=None; tries=0
k=1
while P is None:
    cand=base+k*M; tries+=1
    if sympy.isprime(cand): P=cand
    k+=1
print('     149: modulus has',len(str(M)),'digits; first prime in the class found at k =',k-1,'(',len(str(P)),'digits )')
U2=[q for q in sympy.primerange(3,350) if nr(q,P)]
chk('149: same 32 non-residue primes <= 349; p in a hard class',(U2==U,P%840 in HARD),(True,True))
same=all(leg(-A(q)*P,b)==leg(-A(q)*p0,b) for q in U for b in U if q!=b)
chk('149: the edge relation among the 32 vertices is unchanged (symbol (-A(q) p / b) for all ordered pairs)',same,True)
c,x=src(P,41,1)
chk('149: root 41, t=1 at the new p: c = 123, 461 | x, c | 4*461 + 1',(c,x%461,(4*461+1)%c),(123,0,0))
d=x*x//461
chk('149: Type I at the root of the new p (d = x^2/461): exact identity holds',es(P,c,x,d,'I')[0]==x,True)
# m(e): e | x^2  <=>  m(e) | x
def m_of(e):
    out=1
    for r,v in sympy.factorint(e).items(): out*=r**((v+1)//2)
    return out
bad=sum(((x*x)%e==0)!=(x%m_of(e)==0) for e,x in ((random.randrange(1,5000),random.randrange(1,5000)) for _ in range(20000)))
chk('149: e | x^2 <=> m(e) | x on 20,000 random pairs: disagreements',bad,0)

# ================= >>150 =================
print('--- >>150')
def lower_all_fail(p,q):
    n=0; t=1; At=A(q)
    while 3*At*t*t<p:
        c=At*t*t; x=(p+c)//4; w1,w2=wit(p,c,x)
        if w1 or w2: return n,False,(t,w1[:1],w2[:1])
        n+=1; t+=2
    return n,True,None
for p,name in ((289249,'ex1'),(13217743009,'ex2')):
    chk('150 %s: p prime, hard class, least non-residue prime'%name,(sympy.isprime(p),p%840 in HARD,least_nr(p)),(True,True,17))
    chk('150 %s: base congruences (mod 80, 111, 59, 17)'%name,(p%80,p%111,p%59,p%17),(49,94,31,11))
    chk('150 %s: 17 -(5)-> 37 and 17 -(47)-> 59 are edges'%name,(edge(p,17,5,37),edge(p,17,47,59)),(True,True))
    c,x=src(p,37,1); dd=49*200
    chk('150 %s: vertex 37, t=1: c=111, 20 | x, 111 | p + 4*9800, -c mod 7'%name,(c,x%20,(p+4*dd)%111,(-c)%7),(111,0,0,1))
    c,x=src(p,59,1)
    chk('150 %s: vertex 59, t=1: c=59, 59 | p + 28, -c mod 7'%name,(c,(p+28)%59,(-c)%7),(59,0,4))
    chk('150 %s: p mod 7, and 17 | N = (3p+1)/4'%name,(p%7,((3*p+1)//4)%17,(3*p+1)%4),(2,0,0))
p=289249
chk('150 ex1: lower sources of 37 and of 59 all fail (count, all fail)',(lower_all_fail(p,37)[:2],lower_all_fail(p,59)[:2]),((15,True),(20,True)))
N=(3*p+1)//4; x=8*N//17; e=7*N//17; c=4*x-p
chk('150 ex1: c = (7p+8)/17 = 119103 = 3*29*37^2, x, e',(c,(7*p+8)//17,sympy.factorint(c),x,e),(119103,119103,{3:1,29:1,37:2},102088,89327))
chk('150 ex1: e | x^2, 4e + 1 = 3c, upper source of vertex 29 with t = 37',((x*x)%e,4*e+1-3*c,src(p,29,37)==(c,x),3*c>p),(0,0,True,True))
chk('150 ex1: path 17 -(1)-> 263 -(19)-> 29',(edge(p,17,1,263),edge(p,263,19,29)),(True,True))
chk('150 ex1: every prime factor of x is <= p/17; exact identity',(max(sympy.factorint(x))*17<=p,es(p,c,x,x*x//e,'I')[0]),(True,x))
p=13217743009
t1=time.time()
chk('150 ex2: lower sources of 37 and of 59 all fail (count, all fail)',(lower_all_fail(p,37)[:2],lower_all_fail(p,59)[:2]),((3150,True),(4321,True)))
print('     (%.1f s for the 7,471 lower sources)'%(time.time()-t1))
c0=(7*p+8)//17
chk('150 ex2: (7p+8)/17 = 3*37*41*97*12329',sympy.factorint(c0),{3:1,37:1,41:1,97:1,12329:1})
rows=[(2,1847,11,4,26784,'I'),(9,31,5,4,944124556,'I'),(23,3407,1,4,2345778,'II'),(30,31,19,4,7425625,'II'),(37,135007,1,3,1113783248,'I'),(44,5591,1,4,283237470,'II' if False else 'I')]
okrows=[]
for res,q,t,j,f,typ in rows:
    c,x=src(p,q,t); d=7**j*f
    cong=((4*d+1)%c==0) if typ=='I' else ((4*d+p)%c==0)
    w1,w2=wit(p,c,x)
    okrows.append((nr(q,p) and sympy.isprime(q),3*c<p,(x*x)%f==0,(c*f)%7!=0,cong,(-c)%49==res,not w1 and not w2))
chk('150 ex2: six rows: (q is a vertex, lower, f | x^2, 7 does not divide c f, type congruence, -c mod 49 = stated residue, source fails)',okrows,[(True,)*7]*6)
chk('150 ex2: reachability 17 -(9)-> 31, 17 -(3113)-> 5591, 37 -(939)-> 1847, 37 -(1787)-> 3407, 37 -(12581)-> 135007',
    (edge(p,17,9,31),edge(p,17,3113,5591),edge(p,37,939,1847),edge(p,37,1787,3407),edge(p,37,12581,135007)),(True,)*5)
chk('150 ex2: p mod 49, p mod 101; residues 2 mod 7 not removed by the six rows',(p%49,p%101,sorted(set(range(2,49,7))-{r[0] for r in rows})),(16,67,[16]))
N=(3*p+1)//4; x=50*N//101; e=49*N//101; c=4*x-p
chk('150 ex2: c = (49p+50)/101 = 79167511 * 9^2, x, e',(c,(49*p+50)//101,c==79167511*81,x,e),(c,c,True,4907577850,4809426293))
chk('150 ex2: e | x^2, 4e + 1 = 3c, upper source of vertex 79167511 with t = 9',((x*x)%e,4*e+1-3*c,src(p,79167511,9)==(c,x),3*c>p),(0,0,True,True))
path=[(17,4911,70663),(70663,697,5943308047),(5943308047,1,171080813),(171080813,1,1716373181),(1716373181,1,79167511)]
chk('150 ex2: the five edges of the path',[edge(p,a,t,b) for a,t,b in path],[True]*5)
chk('150 ex2: every prime factor of x is <= p/17; exact identity',(max(sympy.factorint(x))*17<=p,es(p,c,x,x*x//e,'I')[0]),(True,x))

# ================= additions (NOT in the posts) =================
print('--- additions')
for p,(a,h) in ((289249,(7,17)),(13217743009,(49,101))):
    b=(a+h)//3; N=(3*p+1)//4
    print('     p=%d: (a,h,b)=(%d,%d,%d): h*c = a*p + b: %s ; (h/p) = %d, (b/p) = %d, (c/p) = %d'%(p,a,h,b,h*((a*p+b)//h)==a*p+b,leg(h,p),leg(b,p),leg((a*p+b)//h,p)))
# (i) ex2's first candidate c0 = (7p+8)/17 is not a square source, but x = 8N/17, e = 7N/17 is still a decomposition of 4/p
p=13217743009; N=(3*p+1)//4; x=8*N//17; e=7*N//17; c=4*x-p
chk('add: ex2, (a,h,b) = (7,17,8): c = (7p+8)/17, e | x^2, 4e+1 = 3c, exact 4/p identity (c is not of the form A(q) t^2)',
    (c==c0,(x*x)%e,4*e+1-3*c,es(p,c,x,x*x//e,'I')[0]==x),(True,0,0,True))
# (ii) the base condition p = 11 (mod 17) alone: (a,h,b) = (1,17,6)
for p in (289249,13217743009):
    N=(3*p+1)//4; x=6*N//17; e=N//17; c=4*x-p
    chk('add: p=%d, (a,h,b) = (1,17,6): c = (p+6)/17, e | x^2, 4e+1 = 3c, exact 4/p identity'%p,(17*c==p+6,(x*x)%e,4*e+1-3*c,es(p,c,x,x*x//e,'I')[0]==x),(True,0,0,True))
# (iii) how often the k = 3 family applies: N = (3p+1)/4 has a prime factor = 2 (mod 3)  [then (a,h) = (1,h) works]
cnt=Counter()
for p in hard7:
    N=(3*p+1)//4
    cnt['hard primes < 10^7']+=1
    cnt['N has a prime factor = 2 mod 3']+=any(r%3==2 for r in sympy.factorint(N))
print('     add:',dict(cnt))
print(nok,'OK /',ndiff,'DIFF,',round(time.time()-t0,1),'s')
