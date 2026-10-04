"""Exact spin rotations, octupole preparations, response tomography, NMR example."""
from __future__ import annotations
from itertools import product
import sympy as S
from sympy.physics.wigner import wigner_d_small, clebsch_gordan
from spin_foundations import cartesian, spherical, zero, HELICITIES, R, I4

DIRECTIONS=(
    (0,0,1),(1,0,0),(0,1,0),(R(2,3),R(2,3),-R(1,3)),
    (R(1,3),R(2,3),-R(2,3)),(R(2,3),R(1,3),-R(2,3)),
    (R(2,11),R(6,11),-R(9,11)),
)


def x3(n,J):
    axis=sum((n[a]*J[a] for a in range(3)),S.zeros(4))
    return S.sqrt(5)/3*(axis**3-R(41,20)*axis)


def real_octupole_basis(T):
    basis=[T[3,0]]
    for m in range(1,4):
        basis += [(T[3,m]+T[3,m].H)/S.sqrt(2),
                  (T[3,m]-T[3,m].H)/(S.I*S.sqrt(2))]
    return basis


def response_matrix(J,T):
    basis=real_octupole_basis(T)
    return S.Matrix([[S.simplify(S.trace(x3(n,J)*B)) for B in basis]
                     for n in DIRECTIONS])


def run_checks():
    J,Q,O,_,_=cartesian();T=spherical();checks={}
    def check(name,condition,payload=None):
        if condition is not True:raise AssertionError(name)
        checks[name]=payload or {'method':'exact matrix/polynomial identity'}
    phi=S.pi/2;theta=2*S.atan(R(1,2))
    uy=wigner_d_small(R(3,2),-theta)
    uz=S.diag(*(S.exp(-S.I*phi*m) for m in HELICITIES))
    U=uz*uy
    n=(0,R(4,5),R(3,5))
    check('spin.rotation_calibration',zero(U*J[2]*U.H-sum((n[a]*J[a] for a in range(3)),S.zeros(4))))
    check('spin.rotation_coherence',all(zero(uz*T[K,m]*uz.H-S.exp(-S.I*m*phi)*T[K,m])
          for K,m in T) and not zero(U*T[3,1]*U.H-uz.H*(uy*T[3,1]*uy.T)*uz))
    check('spin.rotation_octupole',zero(U*T[3,0]*U.H-x3(n,J)))
    Oamp=S.symbols('O',real=True)
    rotated_component_cases=0
    for n in DIRECTIONS:
        X=x3(n,J)
        check('spin.direction_unit_'+str(DIRECTIONS.index(n)),sum(v*v for v in n)==1)
        for a,b,c in product(range(3),repeat=3):
            derived=S.simplify(S.trace((I4/4+Oamp*S.sqrt(5)/3*X)*O[a,b,c]))
            expected=R(5,2)*Oamp*(n[a]*n[b]*n[c] -
                (n[a]*int(b==c)+n[b]*int(a==c)+n[c]*int(a==b))/5)
            if not zero(derived-expected):raise AssertionError('spin.rotated_cartesian_tensor')
        moments={(a,b,c):S.simplify(S.trace((I4/4+Oamp*S.sqrt(5)/3*X)*O[a,b,c]))
                 for a,b,c in product(range(3),repeat=3)}
        plus=n[0]+S.I*n[1]
        observed=(moments[2,2,2],
                  moments[2,2,0]+S.I*moments[2,2,1],
                  2*(moments[2,0,0]-moments[2,1,1])+4*S.I*moments[2,0,1],
                  moments[0,0,0]-3*moments[0,1,1]+S.I*(3*moments[0,0,1]-moments[1,1,1]))
        predicted=(Oamp*(5*n[2]**3-3*n[2])/2,
                   Oamp*(5*n[2]**2-1)*plus/2,
                   5*Oamp*n[2]*plus**2,
                   R(5,2)*Oamp*plus**3)
        if any(not zero(a-b) for a,b in zip(observed,predicted)):
            raise AssertionError('spin.rotated_components')
        rotated_component_cases+=4
    check('spin.rotated_cartesian_tensor',True,{'directions':len(DIRECTIONS),'components':27})
    check('spin.rotated_components',True,{'exact_cases':rotated_component_cases,'signed_amplitudes':True})
    # At a transverse direction, the LLT and TTT components are both nonzero.
    trans=(1,0,0);X=x3(trans,J);beta=S.sqrt(5)/12
    rho=I4/4+beta*X
    Rmom={(a,b,c):S.simplify(S.trace(rho*O[a,b,c])) for a,b,c in product(range(3),repeat=3)}
    llt=Rmom[2,2,0]+S.I*Rmom[2,2,1]
    ttt=Rmom[0,0,0]-3*Rmom[0,1,1]
    check('spin.transverse_llt_ttt',llt!=0 and ttt!=0 and all(zero(S.trace(rho*op)) for op in J))
    basis=real_octupole_basis(T)
    check('spin.real_octupole_basis',all(zero(S.trace(A*B)-int(a==b))
          for a,A in enumerate(basis) for b,B in enumerate(basis)))
    response=response_matrix(J,T);det=S.simplify(response.det())
    check('spin.seven_direction_tomography',det!=0,
          {'determinant':str(det),'directions':[list(map(str,n)) for n in DIRECTIONS]})
    coeff=S.symbols('w0:7',real=True)
    check('spin.response_recovery',all(zero(v-w) for v,w in zip(response.inv()*(response*S.Matrix(coeff)),coeff)))
    lower=I4/4+J[0]+Q[0,1]+T[1,1]+T[1,1].H
    check('spin.lower_rank_cancellation',all(zero(S.trace(x3(n,J)*lower)) for n in DIRECTIONS))
    p=S.symbols('p0:4',real=True)
    sl=sum(p[a]*HELICITIES[a] for a in range(4))
    ll=p[0]+p[3]-p[1]-p[2]
    lll=R(3,10)*(p[0]-p[3]-3*(p[1]-p[2]))
    sL,sLL,sLLL=S.symbols('S_L S_LL S_LLL',real=True)
    inverse=[(1+R(6,5)*sL+sLL+R(2,3)*sLLL)/4,
             (1+R(2,5)*sL-sLL-2*sLLL)/4,
             (1-R(2,5)*sL-sLL+2*sLLL)/4,
             (1-R(6,5)*sL+sLL-R(2,3)*sLLL)/4]
    solved=S.solve([sum(p)-1,sl-sL,ll-sLL,lll-sLLL],p)
    check('spin.population_inverse',all(zero(inverse[i]-solved[p[i]]) for i in range(4)))
    jp=J[0]+S.I*J[1]
    check('spin.nmr_strengths',[S.simplify(jp[a,a+1]**2) for a in range(3)]==[3,4,3])
    d1,d2,d3,C=S.symbols('d_1 d_2 d_3 C',real=True,nonzero=True)
    low=(1-d1-2*d2-3*d3)/4
    nmr=[low+d3+d2+d1,low+d3+d2,low+d3,low]
    replaced={d1:p[0]-p[1],d2:p[1]-p[2],d3:p[2]-p[3]}
    check('spin.nmr_inverse',all(zero((nmr[i]-p[i]).subs(replaced)-(1-sum(p))/4)
          for i in range(4)) and zero(sum(nmr)-1) and
          [S.simplify(nmr[i]-nmr[i+1]) for i in range(3)]==[d1,d2,d3] and
          [3*C*d1,4*C*d2,3*C*d3]==[3*C*(nmr[0]-nmr[1]),4*C*(nmr[1]-nmr[2]),3*C*(nmr[2]-nmr[3])])
    # Independent CG isometry L=1 ⊗ s=1/2, basis ordered (Lz,sz).
    lvals=(1,0,-1);svals=(R(1,2),-R(1,2))
    V=S.Matrix(6,4,lambda row,col:clebsch_gordan(1,R(1,2),R(3,2),
          lvals[row//2],svals[row%2],HELICITIES[col]))
    spinz=S.diag(*(sz for _l in lvals for sz in svals))
    orbital=S.diag(*((3*l*l-2)*sz for l in lvals for sz in svals))
    check('spin.ls_isometry',zero(V.H*V-I4))
    check('spin.ls_projection',zero(V.H*spinz*V-J[2]/3) and
          zero(V.H*orbital*V-O[2,2,2]-R(2,15)*J[2]))
    return checks
