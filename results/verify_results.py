#!/usr/bin/env python3
"""Reproduce manuscript certificates; figures are generated separately.
Runs with Python's standard library. --rebuild-geometry certifies every edge anew.
"""
import argparse, json, hashlib, platform, time
from fractions import Fraction as F
from pathlib import Path
import certify_catalogue as c
ROOT=Path(__file__).resolve().parent

def repeated_certificate(geo):
    graph=geo['graph'];cycle=[1,26,56,35,40,52,17,53,18,48,43,41,50,42,49,51,33]
    assigned=dict(zip(cycle,cycle[1:]+cycle[:1]));witness=c.complete(graph,assigned)
    assert witness is not None
    cs=[c.K(1),-(2+c.SQ2)*F(1,4)]
    for n in range(2,194):cs.append(2*cs[1]*cs[-1]-cs[-2])
    probs=[(6-c.SQ2)*F(1,17)*(1-x) for x in cs]
    allowed=[];errors=[];seq=[c.READOUT[i] for i in cycle]
    for n,p in enumerate(probs):
        allowed.append({i for i,r in enumerate(c.READOUT) if (c.K(r)-p-F(1,3)).sign()<=0 and (p-c.K(r)-F(1,3)).sign()<=0})
        e=c.K(seq[n%17])-p;errors.append(abs(e.approx()))
        if n<=192:assert (e-F(13,40)).sign()<0 and (-e-F(13,40)).sign()<0
        else:assert (e-F(17,50)).sign()>0 or (-e-F(17,50)).sign()>0
    table=c.ANCHORS.copy();used=set(table.values());counts={'nodes':0,'hall_rejections':0,'band_rejections':0,'injection_rejections':0,'max_depth':0}
    def dfs(n,state):
        counts['nodes']+=1
        if state not in allowed[n]:counts['band_rejections']+=1;return False
        counts['max_depth']=max(counts['max_depth'],n)
        if n==193:return True
        if state in table:return dfs(n+1,table[state])
        for target in graph[state]:
            if target in used:counts['injection_rejections']+=1;continue
            table[state]=target;used.add(target)
            if c.complete(graph,table) is not None:
                if dfs(n+1,target):return True
            else:counts['hall_rejections']+=1
            used.remove(target);del table[state]
        return False
    assert not dfs(0,1)
    assert counts['max_depth']==192
    # Each edge's unique shift constructs a full 464-frame permutation.
    shifts={(i,j):s for i,j,s in geo['edges']}
    assert len(shifts)==len(geo['edges'])
    lift=[8*witness[i]+(k+shifts[i,witness[i]])%8 for i in range(58) for k in range(8)]
    assert len(set(lift))==464
    # Check every lifted transition against physical matrices, including polar phases.
    cat=c.make_catalogue()
    def original(i,k):
        if i==0:return (-k)%8
        if i==1:return 8+k
        return c.REPS[i]+8*k
    for i in range(58):
        for k in range(8):
            out_i,out_k=divmod(lift[8*i+k],8)
            src=cat[original(i,k)];dst=cat[original(out_i,out_k)]
            target=c.mm(c.H,c.mm(c.T,src))
            assert (c.overlap(dst,target)-(F(2)-F(13,50)**2)**2).sign()>=0
    # Full calibrated transitions: P(T^{-1})=H and P(T^{-1}H)=I.
    assert lift[8*1+3]==8*26
    assert lift[8*33]==8*1+4
    assert sum(shifts[i,assigned[i]] for i in cycle)%8==0
    out={'quotient_cycle':cycle,'sigma':witness,'sequence':[str(x) for x in seq], 'max_error_through_192':max(errors[:193]),'error_193':errors[193], 'exhaustive_search':counts,'lift_in_orbit_order':lift}
    (ROOT/'data/repeated_witness.json').write_text(json.dumps(out,indent=2))
    return out

def verify_length12(geo):
    from itertools import product
    item=json.loads((ROOT/'data/length12_feasible.json').read_text())
    sigma={int(i):j for i,j in item['solution'].items()}
    assert set(sigma)==set(range(58)) and len(set(sigma.values()))==58
    assert all(j in geo['graph'][i] for i,j in sigma.items())
    assert all(sigma[i]==j for i,j in c.ANCHORS.items())
    checked=0;maxerr=0
    for n in range(1,13):
        for chars in product('HT',repeat=n):
            w=''.join(chars)
            p=c.ideal(w);e=c.K(c.output(sigma,w))-p
            assert (e-F(1,3)).sign()<=0 and (-e-F(1,3)).sign()<=0
            maxerr=max(maxerr,abs(e.approx()));checked+=1
    return {'words_checked':checked,'max_error':maxerr}

def four_frame():
    # Enumerate the entire four-frame matching class and certify its radii.
    from itertools import permutations
    tables=[p for p in permutations(range(4)) if all(p[i] in (i,(i+1)%4) for i in range(4))]
    assert tables==[(0,1,2,3),(1,2,3,0)]
    # Rational comparisons of exact Q(sqrt(2)) expressions.
    a=(2-c.SQ2)*F(1,4);b=c.SQ2*F(1,4);mix=(c.SQ2-1)*F(1,4)
    assert (a-F(1,4)).sign()<0 and (b-F(1,4)).sign()>0
    assert (mix).sign()>0
    return {'allowed_tables':[list(p) for p in tables],'individual_radius':a.approx().real,'shared_radius':b.approx().real,'mixed_table_radius':mix.approx().real}

def main():
    p=argparse.ArgumentParser();p.add_argument('--rebuild-geometry',action='store_true');args=p.parse_args()
    start=time.time();path=ROOT/'data/catalogue_geometry.json'
    geo=c.build_geometry() if args.rebuild_geometry or not path.exists() else json.loads(path.read_text())
    repeated=repeated_certificate(geo);print('192/193 horizon verified',repeated['exhaustive_search'],flush=True)
    allwords=verify_length12(geo);print('all-word finite-depth witness verified',allwords,flush=True)
    out={'python':platform.python_version(),'platform':platform.platform(),'four_frame':four_frame(),'repeated':repeated['exhaustive_search'],'all_words_12':allwords,'elapsed_seconds':time.time()-start}
    out['source_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),ROOT/'certify_catalogue.py']}
    (ROOT/'data/verification_report.json').write_text(json.dumps(out,indent=2));print('All certificate checks passed.',flush=True)
if __name__=='__main__':main()
