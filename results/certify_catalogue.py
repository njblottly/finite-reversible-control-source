#!/usr/bin/env python3
"""Exact catalogue geometry and shared-controller circuit search.
Python 3.10+, standard library only. Field Q(sqrt(2),sqrt(3),sqrt(5),sqrt(7),i).
All decisions in geometry, probability bands and matching search are exact.
"""
from fractions import Fraction as F
from math import isqrt
from functools import lru_cache
from pathlib import Path
import json, time, hashlib, argparse
ROOT=Path(__file__).resolve().parent
PRIMES=(2,3,5,7)
RADS=[1]*16
for m in range(16):
    for j,p in enumerate(PRIMES):
        if m>>j&1:RADS[m]*=p
MUL=[]
for a in range(32):
    row=[]
    for b in range(32):
        c=RADS[(a&15)&(b&15)]*(-1 if a&16 and b&16 else 1)
        row.append((a^b,c))
    MUL.append(row)
class K:
    def __init__(self,x=0):self.d=x if isinstance(x,dict) else ({0:F(x)} if x else {})
    def __add__(a,b):
        b=b if isinstance(b,K) else K(b);d=a.d.copy()
        for m,v in b.d.items():
            d[m]=d.get(m,0)+v
            if not d[m]:del d[m]
        return K(d)
    __radd__=__add__
    def __neg__(a):return K({m:-v for m,v in a.d.items()})
    def __sub__(a,b):return a+-asK(b)
    def __rsub__(a,b):return asK(b)+-a
    def __mul__(a,b):
        b=asK(b);d={}
        for m,x in a.d.items():
            for n,y in b.d.items():
                k,c=MUL[m][n];d[k]=d.get(k,0)+x*y*c
        return K({m:v for m,v in d.items() if v})
    __rmul__=__mul__
    def conj(a):return K({m:(-v if m&16 else v) for m,v in a.d.items()})
    def sign(a):
        assert not any(m&16 for m in a.d)
        if not a.d:return 0
        for bits in (48,96,192,384,768):
            lo,hi=bounds(a,bits)
            if lo>0:return 1
            if hi<0:return -1
        raise ArithmeticError('Unresolved sign')
    def approx(a):
        return sum(float(v)*(RADS[m&15]**.5)*(1j if m&16 else 1) for m,v in a.d.items())
def asK(x):return x if isinstance(x,K) else K(x)
@lru_cache(None)
def root_interval(m,bits):
    den=1<<bits;z=isqrt(RADS[m]*den*den)
    return F(z,den), F(z if z*z==RADS[m]*den*den else z+1,den)
def bounds(a,bits):
    lo=hi=F(0)
    for m,v in a.d.items():
        l,u=root_interval(m,bits)
        lo+=v*(l if v>=0 else u);hi+=v*(u if v>=0 else l)
    return lo,hi
def sqrt_int(n):
    if n==0:return K()
    coeff=1;mask=0
    for j,p in enumerate(PRIMES):
        k=0
        while n%p==0:n//=p;k+=1
        coeff*=p**(k//2)
        if k%2:mask|=1<<j
    assert n==1
    return K({mask:F(coeff)})
ONE=K(1);ZERO=K();II=K({16:F(1)});SQ2=sqrt_int(2)
z=(ONE+II)*SQ2*F(1,2)
zs=[ONE]
for j in range(1,8):zs.append(zs[-1]*z)
I=(ONE,ZERO,ZERO,ONE)
H=tuple(x*SQ2*F(1,2) for x in (ONE,ONE,ONE,-ONE))
T=(ONE,ZERO,ZERO,z)
def mm(a,b):return tuple(sum((a[2*i+k]*b[2*k+j] for k in range(2)),K()) for i in range(2) for j in range(2))
def overlap(a,b):
    tr=sum((x.conj()*y for x,y in zip(a,b)),K());return tr.conj()*tr
def make_catalogue():
    out=[]
    for g in range(8):out.append((ZERO,ONE,zs[g],ZERO))
    for g in range(8):out.append((ONE,ZERO,ZERO,-zs[g]))
    for n in range(1,8):
        a=sqrt_int(2*n)*F(1,4);b=sqrt_int(2*(8-n))*F(1,4)
        for beta in range(8):
            for gamma in range(8):out.append((a,b*zs[beta],b*zs[gamma],-a*zs[(beta+gamma)%8]))
    return out
REPS=[0,8]+[16+64*u+g for u in range(7) for g in range(8)]
def decode_label(j):
    if j<8:return 0,(-j)%8
    if j<16:return 1,j-8
    u,r=divmod(j-16,64);b,g=divmod(r,8);return 2+8*u+g,b
READOUT=[F(1),F(0)]+[F(7-u,8) for u in range(7) for g in range(8)]
TPERM=[0,1]+[2+8*u+(g+1)%8 for u in range(7) for g in range(8)]
TINVERSE=[0,1]+[2+8*u+(g-1)%8 for u in range(7) for g in range(8)]
ANCHORS={1:26,33:1}
def build_geometry():
    cat=make_catalogue();threshold=F(2)-F(13,50)**2;threshold=threshold**2
    edges=[];graph=[set() for _ in REPS]
    started=time.time()
    for i,idx in enumerate(REPS):
        target=mm(H,mm(T,cat[idx]))
        for j,a in enumerate(cat):
            if (overlap(a,target)-threshold).sign()>=0:
                q,s=decode_label(j)
                if i in ANCHORS and q!=ANCHORS[i]:continue
                if q in ANCHORS.values() and ANCHORS.get(i)!=q:continue
                edges.append([i,q,s]);graph[i].add(q)
        if i%10==0:print('geometry rows',i+1,'of 58',flush=True)
    assert len(edges)==174,(len(edges),sum(map(len,graph)))
    assert [1,26,5] in edges and [33,1,4] in edges
    # Calibration fixes shifts as well as target orbits.
    edges=[e for e in edges if e[0] not in ANCHORS or e in [[1,26,5],[33,1,4]]]
    payload={'edges':edges,'graph':[sorted(x) for x in graph], 'readout':[str(x) for x in READOUT], 'geometry_seconds':time.time()-started}
    (ROOT/'data/catalogue_geometry.json').write_text(json.dumps(payload,indent=2))
    return payload

def complete(graph,assigned,certificate=False):
    used=set(assigned.values());match={v:k for k,v in assigned.items()}
    locked=set(assigned)
    def aug(i,seen):
        for j in graph[i]:
            if j in seen:continue
            seen.add(j)
            if j not in match or (match[j] not in locked and aug(match[j],seen)):
                match[j]=i;return True
        return False
    for i in range(58):
        if i not in locked and not aug(i,set()):return None
    return {i:j for j,i in match.items()}

@lru_cache(None)
def ideal(word):
    v=(ONE,ZERO)
    for g in word:
        a=H if g=='H' else T
        v=(a[0]*v[0]+a[1]*v[1],a[2]*v[0]+a[3]*v[1])
    return v[1].conj()*v[1]
def acceptable(word,delta):
    p=ideal(word)
    return {i for i,r in enumerate(READOUT) if (K(r)-p-delta).sign()<=0 and (p-K(r)-delta).sign()<=0}

def solve_suite(graph,words,delta=F(1,3),max_nodes=0):
    """Exact shared-table search. Every leaf must have a full matching extension."""
    allowed=[acceptable(w,delta) for w in words]
    assigned=ANCHORS.copy();used=set(assigned.values());nodes=0;halls=0
    def run(k,pos,state):
        nonlocal nodes,halls
        nodes+=1
        if max_nodes and nodes>max_nodes:raise TimeoutError('Search node limit; no impossibility claim')
        if k==len(words):return complete(graph,assigned)
        w=words[k]
        if pos==len(w):
            if state not in allowed[k]:return None
            return run(k+1,0,1)
        if w[pos]=='T':return run(k,pos+1,TPERM[state])
        src=TINVERSE[state]
        if src in assigned:return run(k,pos+1,assigned[src])
        for dst in graph[src]:
            if dst in used:continue
            if pos==len(w)-1 and dst not in allowed[k]:continue
            assigned[src]=dst;used.add(dst)
            if complete(graph,assigned) is not None:
                result=run(k,pos+1,dst)
                if result is not None:return result
            else:halls+=1
            used.remove(dst);del assigned[src]
        return None
    sol=run(0,0,1)
    return sol,{'nodes':nodes,'hall_rejections':halls}
def output(sigma,word):
    state=1
    for c in word:state=TPERM[state] if c=='T' else sigma[TINVERSE[state]]
    return READOUT[state]

def run_search(geometry,max_length):
    from itertools import product
    graph=geometry['graph'];candidates=[]
    # Final T does not change the computational readout, and initial T fixes |0>.
    for n in range(1,max_length+1):
        for chars in product('HT',repeat=n):
            w=''.join(chars)
            if w[0]!='H' or w[-1]!='H':continue
            candidates.append(w)
    words=[];records=[]
    for iteration in range(100):
        try:
            solution,counts=solve_suite(graph,words,max_nodes=2000000)
        except TimeoutError as exc:
            result={'status':'unresolved_node_limit','max_length':max_length,'circuits':words,'records':records,'message':str(exc)}
            (ROOT/'data/circuit_search.json').write_text(json.dumps(result,indent=2))
            print(str(exc),flush=True)
            return result
        print('iteration',iteration,'tests',len(words),counts,'feasible',solution is not None,flush=True)
        records.append({'words':list(words),'feasible':solution is not None,**counts})
        if solution is None:break
        violating=[]
        for w in candidates:
            diff=K(output(solution,w))-ideal(w)
            if (diff-F(1,3)).sign()>0 or (-diff-F(1,3)).sign()>0:
                violating.append((len(w),-abs(float(diff.approx().real)),w))
        if not violating:
            print('all candidate circuits feasible under one table',flush=True);break
        words.append(min(violating)[2]);print('add',words[-1],flush=True)
    result={'max_length':max_length,'circuits':words,'records':records,'feasible':solution is not None,'solution':solution}
    if solution is None:
        # Greedy deletion gives an inclusion-minimal suite, not minimum cardinality.
        for word in words[:]:
            trial=[x for x in words if x!=word]
            witness,_=solve_suite(graph,trial,max_nodes=2000000)
            if witness is None:words=trial
        result['reduced_circuits']=words
        witnesses={}
        for w in words:
            sol,ct=solve_suite(graph,[x for x in words if x!=w],max_nodes=2000000)
            assert sol is not None
            witnesses[w]={'sigma':sol,'counts':ct}
        result['deletion_witnesses']=witnesses
        _,result['final_counts']=solve_suite(graph,words,max_nodes=2000000)
    (ROOT/'data/circuit_search.json').write_text(json.dumps(result,indent=2))
    return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--rebuild',action='store_true');p.add_argument('--search',type=int,default=0);args=p.parse_args()
    f=ROOT/'data/catalogue_geometry.json'
    geo=build_geometry() if args.rebuild or not f.exists() else json.loads(f.read_text())
    print('exact geometric edges',len(geo['edges']),'quotient pairs',sum(map(len,geo['graph'])),flush=True)
    if args.search:run_search(geo,args.search)
