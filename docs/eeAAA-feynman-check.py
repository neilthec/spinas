#!/usr/bin/env python3
"""Exact fixed-reference checks, independent of C++ (requires SymPy).
Run: python3 docs/eeAAA-feynman-check.py
"""
from itertools import permutations
from sympy import Matrix as M, I, Rational as R, sqrt, expand as simplify
m=R(3)
p={1:(5,0,0,4),2:(5,0,0,-4),
   3:(R(5,2),R(5,2),0,0),4:(R(10,3),0,R(10,3),0),
   5:(R(25,6),-R(5,2),-R(10,3),0)}
orders=list(permutations((3,4,5)))
def lower(p):
    e,x,y,z=p
    return M([[e+z,x-I*y],[x+I*y,e-z]])
def adj(a): return M([[a[1,1],-a[0,1]],[-a[1,0],a[0,0]]])
def row(a): return M([[a[1],-a[0]]])
def val(a): return simplify(a[0])
lo={i:lower(v) for i,v in p.items()}; up={i:adj(a) for i,a in lo.items()}
a={}; s={}; sr={}
for i in (3,4,5):
    e,x,y,_=p[i]
    a[i]=M([1,(x+I*y)/e])
    s[i]=M([[e,x-I*y]])
    sr[i]=M([s[i][1],-s[i][0]])
    assert simplify(a[i]*s[i]-lo[i])==M.zeros(2)
def ext(i,spin):
    c,t=(1,0) if i==1 else (0,1)
    big,small=R(3),R(1)
    if spin==-1:
        return (M([big*c,big*t]),M([-small*c,-small*t]),
                M([[big*t,-big*c]]),M([[small*t,-small*c]]))
    return (M([-small*t,small*c]),M([big*t,-big*c]),
            M([[small*c,small*t]]),M([[big*c,big*t]]))
cyclic={i:(a[j],row(a[j])) for i,j in ((3,4),(4,5),(5,3))}
physical={i:(-lo[j]*sr[i],s[i]*up[j]) for i,j in ((3,2),(4,1),(5,1))}
mixed={3:physical[3],4:cyclic[4],5:physical[5]}
reverse={i:(a[j],row(a[j])) for i,j in ((3,5),(4,3),(5,4))}
def eps(q,ward=0):
    result={}
    for i in (3,4,5):
        col,bra=q[i]; den=val(bra*a[i])
        e,f=sqrt(2)*col*s[i]/den,sqrt(2)*sr[i]*bra/den
        if i==ward: e,f=lo[i],up[i]
        result[i]=M.vstack(M.hstack(M.zeros(2),e),M.hstack(f,M.zeros(2)))
    return result
def diagram(q,o,s1,s2,ward=0):
    i,j,k=o; qa,qc=lo[1]-lo[i],lo[2]-lo[k]; eye=M.eye(2)
    pa=M.vstack(M.hstack(m*eye,-qa),M.hstack(-adj(qa),m*eye))
    pc=M.vstack(M.hstack(m*eye,qc),M.hstack(adj(qc),m*eye))
    ar,br,_,_=ext(1,s1); _,_,al,bl=ext(2,s2)
    v=eps(q,ward); den=(qa.det()-m*m)*(qc.det()-m*m)
    # Expand after each multiplication to keep exact rational arithmetic small.
    x=M.vstack(ar,br)
    for factor in (v[i],pa,v[j],pc,v[k]):
        x=(factor*x).applyfunc(simplify)
    return val(-M.hstack(al,bl)*x/den)

def expanded(q,o,s1,s2):
    i,j,k=o; qa,qc=lo[1]-lo[i],lo[2]-lo[k]
    ar,br,_,_=ext(1,s1); _,_,al,bl=ext(2,s2)
    ac,ab=q[i]; bc,bb=q[j]; cc,cb=q[k]
    L,S=val(al*cc),val(bl*sr[k]); T,C=val(s[i]*br),val(ab*ar)
    U,V=val(s[k]*adj(qc)*bc),val(s[j]*adj(qa)*ac)
    W,Z=val(bb*qa*sr[i]),val(cb*qc*sr[j])
    kj,ji=val(s[k]*sr[j]),val(s[j]*sr[i]); rr,rc=val(bb*ac),val(cb*bc)
    n=(-L*U*V*T+m*L*U*ji*C+m*m*L*kj*rr*T-m*L*kj*W*C
       -m*S*rc*V*T+m*m*S*rc*ji*C+m*S*Z*rr*T-S*Z*W*C)
    den=(qa.det()-m*m)*(qc.det()-m*m)
    for t in (3,4,5): den*=val(q[t][1]*a[t])
    return simplify(-2*sqrt(2)*n/den)
def h(s1,s2):
    _,_,al,_=ext(1,s1); ar,_,_,_=ext(2,s2)
    den=1; total=0
    for i in (3,4,5): den*=((lo[1]-lo[i]).det()-m*m)
    for k in (3,4,5):
        i,j=[t for t in (3,4,5) if t!=k]
        total+=val(s[i]*sr[j])**2*val(s[k]*up[1]*lo[2]*sr[k])/((lo[2]-lo[k]).det()-m*m)
    return simplify(2*sqrt(2)*m*val(al*ar)*total/den)
def physical_factored(o,s1,s2):
    # The physical X,Y,A,B expressions in the written derivation, with no qi.
    i,j,k=o; r={3:2,4:1,5:1}
    qa,qc=lo[1]-lo[i],lo[2]-lo[k]
    ar,br,_,_=ext(1,s1); _,_,al,bl=ext(2,s2)
    T=val(s[i]*br); C=val(s[i]*up[r[i]]*ar)
    X=val(s[j]*adj(qa)*lo[r[i]]*sr[i])*T+m*val(s[j]*sr[i])*C
    Y=-m*val(s[j]*up[r[j]]*lo[r[i]]*sr[i])*T-val(s[j]*up[r[j]]*qa*sr[i])*C
    L=val(s[k]*up[r[k]]*ext(2,s2)[0]); S=val(bl*sr[k])
    A=L*val(s[k]*adj(qc)*lo[r[j]]*sr[j])-m*S*val(s[k]*up[r[k]]*lo[r[j]]*sr[j])
    B=-m*L*val(s[k]*sr[j])+S*val(s[k]*up[r[k]]*qc*sr[j])
    den=(qa.det()-m*m)*(qc.det()-m*m)
    for t in (3,4,5): den*=val(s[t]*up[r[t]]*a[t])
    return simplify(-2*sqrt(2)*(A*X+B*Y)/den)

assert lo[1]+lo[2]==lo[3]+lo[4]+lo[5]
assert lo[1].det()==lo[2].det()==m*m
checks=0
for s1 in (-1,1):
    for s2 in (-1,1):
        for name,q in [('cyclic',cyclic),('physical',physical),('mixed',mixed),('reverse',reverse)]:
            terms=[]
            for o in orders:
                t=diagram(q,o,s1,s2)
                assert simplify(t-expanded(q,o,s1,s2))==0,(name,o,s1,s2)
                checks+=1; terms.append(t)
                if name=='physical':
                    assert simplify(t-physical_factored(o,s1,s2))==0,(o,s1,s2)
                    checks+=1
            assert simplify(sum(terms)-h(s1,s2))==0,(name,s1,s2)
            checks+=1
            for leg in (3,4,5):
                assert simplify(sum(diagram(q,o,s1,s2,leg) for o in orders))==0,(name,leg,s1,s2)
                checks+=1
        print(f'spins {s1:+d},{s2:+d}: exact matrix/expanded, gauge and H checks PASS',flush=True)
print(f'PASS: {checks} exact assertions; six diagrams, four reference assignments, all fermion spins.')

# Negative control: repeat the old mistake of moving references with the order.
bad=0
for o in orders:
    moved={i:(a[j],row(a[j])) for i,j in zip(o,(o[1],o[2],o[0]))}
    bad+=diagram(moved,o,-1,-1)
assert simplify(bad-h(-1,-1))!=0
print('PASS: moving the cyclic references with each diagram is rejected by the exact negative control.')
