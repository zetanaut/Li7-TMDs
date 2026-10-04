"""Legacy exact spin, transverse-tensor, and selected SIDIS calculations.

The finite algebra and check names are preserved from the original validator.
Importing this module performs no calculations or writes.
"""
from __future__ import annotations
import itertools as it
import sympy as s
from sympy.physics.matrices import mgamma


def run_checks(*, emit=None) -> tuple[dict[str, object], dict[str, int]]:
    """Return the legacy report and computed tensor-map ranks."""
    I, R = s.I, s.Rational
    checks: dict[str, object] = {}
    computed_ranks: dict[str, int] = {}

    def zero(x):
        if isinstance(x, s.MatrixBase):
            return all(s.simplify(e) == 0 for e in x)
        return s.simplify(x) == 0

    def check(name: str, condition: bool) -> None:
        if not condition:
            raise AssertionError(name)
        checks[name] = 'PASS'
        if emit is not None:
            emit('PASS: ' + name)

    # Spin moments and inverse density-matrix map.
    Jz = s.diag(R(3,2),R(1,2),-R(1,2),-R(3,2))
    Jp = s.Matrix([[0,s.sqrt(3),0,0],[0,0,2,0],[0,0,0,s.sqrt(3)],[0,0,0,0]])
    J = [(Jp+Jp.T)/2,(Jp-Jp.T)/(2*I),Jz]
    Q = {(i,j):s.simplify((J[i]*J[j]+J[j]*J[i])/2-R(5,4)*int(i==j)*s.eye(4))
         for i,j in it.product(range(3),repeat=2)}
    O={}
    for i,j,k in it.product(range(3),repeat=3):
        sym=sum((J[a]*J[b]*J[c] for a,b,c in it.permutations((i,j,k))),s.zeros(4))/6
        sub=R(41,60)*(int(i==j)*J[k]+int(i==k)*J[j]+int(j==k)*J[i])
        O[i,j,k]=s.simplify(sym-sub)
    check('three-dimensional quadrupole trace',zero(sum((Q[i,i] for i in range(3)),s.zeros(4))))
    check('three-dimensional octupole traces',all(zero(sum((O[i,i,k] for i in range(3)),s.zeros(4))) for k in range(3)))
    herm=[]
    for a in range(4):
        m=s.zeros(4);m[a,a]=1;herm.append(m)
    for a in range(4):
        for b in range(a+1,4):
            m=s.zeros(4);m[a,b]=m[b,a]=1;herm.append(m)
            m=s.zeros(4);m[a,b]=I;m[b,a]=-I;herm.append(m)
    for rho in herm:
        rec=s.trace(rho)*s.eye(4)/4
        rec+=sum((s.trace(rho*J[i])*J[i] for i in range(3)),s.zeros(4))/5
        rec+=sum((s.trace(rho*op)*op for op in Q.values()),s.zeros(4))/6
        rec+=R(2,9)*sum((s.trace(rho*op)*op for op in O.values()),s.zeros(4))
        if not zero(rho-rec):
            raise AssertionError('spin-density inverse map')
    check('spin-density inverse map on all 16 Hermitian basis matrices',True)
    check('octupole longitudinal eigenvalues',O[2,2,2]==s.diag(R(3,10),-R(9,10),R(9,10),-R(3,10)))
    check('octupole helicity-one matrix',zero(O[2,2,0]+I*O[2,2,1]-s.Matrix([[0,2*s.sqrt(3)/5,0,0],[0,0,-R(6,5),0],[0,0,0,2*s.sqrt(3)/5],[0,0,0,0]])))
    check('quadrupole helicity-two matrix',zero(Q[0,0]-Q[1,1]+2*I*Q[0,1]-Jp**2))

    # STF_2 tensor component reconstruction; u+i v is its independent complex component.
    def tensor(n:int,u,v=0):
        if n==0:return {():s.sympify(u)}
        def component(idx):
            b=sum(idx)
            return (-1)**(b//2)*(u if b%2==0 else v)
        return {idx:s.sympify(component(idx)) for idx in it.product(range(2),repeat=n)}

    def ktensor(n:int,x,y):
        if n==0:return {():s.Integer(1)}
        z=s.expand((x+I*y)**n)
        return tensor(n,s.re(z)/2**(n-1),s.im(z)/2**(n-1))

    def dual(a):
        out={}
        for idx in a:
            j=(1-idx[0],)+idx[1:]
            out[idx]=(1 if idx[0]==0 else -1)*a[j]
        return out

    def contract(a,b):return s.expand(sum(a[i]*b[i] for i in a))

    def vec_low(a,m,k):
        return s.Matrix([sum(a[(i,)+idx]*k[m-1][idx] for idx in it.product(range(2),repeat=m-1)) for i in range(2)])

    def vec_high(a,m,k):
        return s.Matrix([sum(k[m+1][(i,)+idx]*a[idx] for idx in it.product(range(2),repeat=m)) for i in range(2)])

    def stf2(mat):return s.simplify((mat+mat.T)/2-s.eye(2)*s.trace(mat)/2)

    def mat(a):return s.Matrix(2,2,lambda i,j:a[i,j])

    def glu_low(a,m,k):
        if m==1:
            v=s.Matrix([a[(0,)],a[(1,)]])
            q=s.Matrix([k[1][(0,)],k[1][(1,)]])
            return stf2(v*q.T)
        if m==2:return mat(a)
        if m==3:
            return s.Matrix(2,2,lambda i,j:sum(a[i,j,b]*k[1][(b,)] for b in range(2)))
        raise ValueError(m)

    def glu_high(a,m,k):
        return s.Matrix(2,2,lambda i,j:sum(k[m+2][(i,j)+idx]*a[idx] for idx in it.product(range(2),repeat=m)))

    x,y,a,b,c,d=s.symbols('x y a b c d',real=True)
    k={n:ktensor(n,x,y) for n in range(6)}
    E=s.Matrix([[0,1],[-1,0]])
    check('STF momentum tensors through rank five',all(all(zero(sum(t[(j,j)+tail] for j in range(2))) for tail in it.product(range(2),repeat=n-2)) for n,t in k.items() if n>=2))
    A=mat(tensor(2,a,b));B=mat(tensor(2,c,d))
    check('STF of rank-two product vanishes',zero(stf2(A*B)))
    check('dual STF of rank-two product vanishes',zero(stf2(E*A*B)))
    A3,B3=tensor(3,a,b),tensor(3,c,d)
    C=s.Matrix(2,2,lambda i,j:sum(A3[(i,u,v)]*B3[(j,u,v)] for u,v in it.product(range(2),repeat=2)))
    check('dual STF of rank-three double contraction vanishes',zero(stf2(E*C)))
    for m in (1,2,3):
        am=tensor(m,a,b);zp=x+I*y;ap=a+I*b
        vl=vec_low(am,m,k);vh=vec_high(am,m,k)
        check(f'quark complex low branch m={m}',zero(vl[0]+I*vl[1]-ap*s.conjugate(zp)**(m-1)))
        check(f'quark complex high branch m={m}',zero(vh[0]+I*vh[1]-s.conjugate(ap)*zp**(m+1)/2))
        gl=glu_low(am,m,k);gh=glu_high(am,m,k)
        ll={1:ap*zp/2,2:ap,3:ap*s.conjugate(zp)}[m]
        check(f'gluon complex low branch m={m}',zero(gl[0,0]+I*gl[0,1]-ll))
        check(f'gluon complex high branch m={m}',zero(gh[0,0]+I*gh[0,1]-s.conjugate(ap)*zp**(m+2)/4))
        check(f'sine invariant m={m}',zero(-contract(k[m],dual(am))-s.im(s.conjugate(ap)*zp**m)))

    # Independent real rank at a nonzero transverse momentum for 16 target coordinates.
    coords=[(K,m,part) for K in range(4) for m in range(K+1) for part in range(1 if m==0 else 2)]
    ks={n:ktensor(n,s.Integer(2),s.Integer(1)) for n in range(6)}

    def catalogue(r:int):
        out=[]
        for K in range(4):
            for m in range(K+1):
                if m>0 or K%2==0:out.append(('f',K,m,m))
                if m>0 or K%2==1:out.append(('g',K,m,m))
                for n in ([r] if m==0 else [abs(m-r),m+r]):out.append(('h',K,m,n))
        return out

    def basis_column(label,r):
        ch,K,m,n=label;out=[]
        for Kc,mc,part in coords:
            value=s.zeros(4,1)
            if (Kc,mc)==(K,m):
                am=tensor(m,1 if part==0 else 0,1 if part==1 else 0)
                if ch in ('f','g'):
                    ordinary=(K%2==0 if ch=='f' else K%2==1)
                    scalar=contract(ks[m],am) if ordinary else -contract(ks[m],dual(am))
                    value[0 if ch=='f' else 1]=scalar
                elif r==1:
                    if m==0:v=s.Matrix([ks[1][(i,)] for i in range(2)])
                    else:v=vec_low(am,m,ks) if n==m-1 else vec_high(am,m,ks)
                    if K%2==0:v=E*v
                    value[2:4,0]=v
                else:
                    if m==0:v=mat(ks[2])
                    else:v=glu_low(am,m,ks) if n==abs(m-2) else glu_high(am,m,ks)
                    if K%2==1:v=E*v
                    value[2]=v[0,0];value[3]=v[0,1]
            out.extend(value)
        return s.Matrix(out)
    for r,kind in ((1,'quark'),(2,'gluon')):
        cat=catalogue(r)
        rank=s.Matrix.hstack(*(basis_column(l,r) for l in cat)).rank()
        computed_ranks[kind] = int(rank)
        check(f'{kind} catalogue has 32 entries',len(cat)==32)
        check(f'{kind} 64-by-32 tensor map has rank 32',rank==32)
        counts=[sum(l[1]==K for l in cat) for K in range(4)]
        check(f'{kind} multipole counts 2,6,10,14',counts==[2,6,10,14])
        checks[kind+'_catalogue']=[{'channel':ch,'K':K,'m':m,'orbital_rank':n} for ch,K,m,n in cat]

    # Gamma algebra: projections and leading-power electromagnetic hard tensor.
    g=[mgamma(i) for i in range(4)];g5=mgamma(5)
    gp=(g[0]+g[3])/s.sqrt(2);gm=(g[0]-g[3])/s.sqrt(2)
    def isig(A,B):return -(A*B-B*A)/2
    pq=[gp,gp*g5]+[isig(g[i],gp)*g5 for i in (1,2)]
    pd=[gm,gm*g5]+[isig(g[i],gm)*g5 for i in (1,2)]
    bq=[gm/2,g5*gm/2]+[isig(g[i],gm)*g5/2 for i in (1,2)]
    bd=[gp/2,g5*gp/2]+[isig(g[i],gp)*g5/2 for i in (1,2)]
    check('quark projection Gram matrix',zero(s.Matrix([[s.trace(B*P)/2 for P in pq] for B in bq])-s.eye(4)))
    check('fragmentation projection Gram matrix',zero(s.Matrix([[s.trace(B*P)/2 for P in pd] for B in bd])-s.eye(4)))
    expected={(0,0):s.eye(2),(1,0):I*E,(2,2):s.diag(-1,1),(2,3):s.Matrix([[0,-1],[-1,0]]),(3,2):s.Matrix([[0,-1],[-1,0]]),(3,3):s.diag(1,-1)}
    for a0 in range(4):
        for b0 in (0,2,3):
            W=s.Matrix(2,2,lambda i,j:s.simplify(s.trace(bq[a0]*g[i+1]*bd[b0]*g[j+1])))
            check(f'electromagnetic trace channel {a0},{b0}',zero(W-expected.get((a0,b0),s.zeros(2))))
    F,G,Tx,Ty,D,H,px,py,eps,lam,dep=s.symbols('F G Tx Ty D H px py eps lam dep',real=True)
    Tv=s.Matrix([Tx,Ty]);Dv=-E*s.Matrix([px,py])*H
    W=s.eye(2)*F*D+I*E*G*D-(Tv*Dv.T+Dv*Tv.T-s.eye(2)*(Tv.dot(Dv)))
    L=(s.eye(2)+eps*s.diag(1,-1)-I*lam*dep*E)/2
    reduced=sum(L[i,j]*W[i,j] for i,j in it.product(range(2),repeat=2))
    check('SIDIS scalar-helicity-Collins master contraction',zero(reduced-(F*D+lam*dep*G*D+eps*H*(Tx*py+Ty*px))))

    # Population and ideal resolved-transition identities.
    p0,p1,p2,p3=s.symbols('p0 p1 p2 p3',real=True)
    d1,d2,d3=p0-p1,p1-p2,p2-p3
    SL=R(3,2)*p0+R(1,2)*p1-R(1,2)*p2-R(3,2)*p3
    LL=p0+p3-p1-p2;LLL=R(3,10)*(p0-p3-3*(p1-p2))
    check('ideal transition vector moment',zero(SL-(3*d1+4*d2+3*d3)/2))
    check('ideal transition quadrupole moment',zero(LL-(d1-d3)))
    check('ideal transition octupole moment',zero(LLL-R(3,10)*(d1-2*d2+d3)))
    return checks, computed_ranks
