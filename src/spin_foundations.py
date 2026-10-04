"""Independent Cartesian and CG constructions for a spin-3/2 target.

Conventions: descending J_z basis, active U_z(phi)=exp(-i phi J_z).
All arithmetic is exact. Source: eq:J,Q,QO,trace_metrics,rho,spherical.
"""
from __future__ import annotations

from itertools import permutations, product
import sympy as S
from sympy.physics.wigner import clebsch_gordan

R = S.Rational
I4 = S.eye(4)
HELICITIES = (R(3, 2), R(1, 2), -R(1, 2), -R(3, 2))


def zero(value):
    entries = value if isinstance(value, S.MatrixBase) else (value,)
    return all(S.simplify(x) == 0 for x in entries)


def cartesian():
    jp = S.zeros(4)
    for a in range(3):
        m = HELICITIES[a + 1]
        jp[a, a + 1] = S.sqrt(R(3, 2) * R(5, 2) - m * (m + 1))
    jm = jp.T
    J = ((jp + jm) / 2, (jp - jm) / (2 * S.I), S.diag(*HELICITIES))
    cq, co = S.symbols('c_q c_o', real=True)
    q_raw = {(a, b): (J[a] * J[b] + J[b] * J[a]) / 2
             for a, b in product(range(3), repeat=2)}
    o_raw = {(a, b, c): sum((J[i] * J[j] * J[k]
                             for i, j, k in permutations((a, b, c))), S.zeros(4)) / 6
             for a, b, c in product(range(3), repeat=3)}
    q_trace = sum((q_raw[a, a] - cq * I4 for a in range(3)), S.zeros(4))
    o_trace = sum((o_raw[a, a, 2] - co * (J[2] + 2 * int(a == 2) * J[a])
                   for a in range(3)), S.zeros(4))
    q_const = S.solve(S.trace(q_trace) / 4, cq)[0]
    o_const = S.solve(S.trace(o_trace * J[2]) / S.trace(J[2] ** 2), co)[0]
    Q = {(a, b): q_raw[a, b] - q_const * int(a == b) * I4
         for a, b in q_raw}
    O = {(a, b, c): o_raw[a, b, c] - o_const *
         (int(a == b) * J[c] + int(a == c) * J[b] + int(b == c) * J[a])
         for a, b, c in o_raw}
    return J, Q, O, q_const, o_const


def spherical():
    """CG coefficients constructed without using Cartesian spin operators."""
    j = R(3, 2)
    return {(K, m): S.sqrt(R(2 * K + 1, 4)) * S.Matrix(4, 4, lambda row, col:
            clebsch_gordan(j, K, j, HELICITIES[col], m, HELICITIES[row]))
            for K in range(4) for m in range(-K, K + 1)}


def hermitian_basis():
    out = []
    for a in range(4):
        x = S.zeros(4); x[a, a] = 1; out.append(x)
    for a in range(4):
        for b in range(a + 1, 4):
            x = S.zeros(4); x[a, b] = x[b, a] = 1; out.append(x)
            x = S.zeros(4); x[a, b] = S.I; x[b, a] = -S.I; out.append(x)
    return out


def reconstruct(X, J, Q, O, *, q_weight=R(1, 6), o_weight=R(2, 9)):
    rec = S.trace(X) * I4 / 4
    rec += sum((S.trace(X * op) * op for op in J), S.zeros(4)) / 5
    rec += q_weight * sum((S.trace(X * op) * op for op in Q.values()), S.zeros(4))
    rec += o_weight * sum((S.trace(X * op) * op for op in O.values()), S.zeros(4))
    return rec


def run_checks():
    J, Q, O, cq, co = cartesian()
    T = spherical()
    checks = {}
    def check(name, condition):
        if condition is not True:
            raise AssertionError(name)
        checks[name] = {'method': 'exact SymPy matrix identity'}
    check('spin.ladder_commutator', all(zero(J[a] * J[b] - J[b] * J[a] -
          S.I * sum((S.LeviCivita(a, b, c) * J[c] for c in range(3)), S.zeros(4)))
          for a, b in product(range(3), repeat=2)))
    check('spin.casimir', zero(sum((x*x for x in J), S.zeros(4)) - R(15, 4)*I4))
    check('spin.derived_subtractions', cq == R(5,4) and co == R(41,60))
    check('spin.traces', all(zero(sum((Q[a,a] for a in range(3)), S.zeros(4))) for _ in (0,))
          and all(zero(sum((O[a,a,c] for a in range(3)), S.zeros(4))) for c in range(3)))
    check('spin.spherical_orthonormal', all(zero(S.trace(A.H*B) - (key==other))
          for key,A in T.items() for other,B in T.items()))
    check('spin.spherical_conjugation', all(zero(T[K,m].H-(-1)**m*T[K,-m])
          for K,m in T))
    check('spin.anchor_conversion', zero(T[0,0]-I4/2) and
          zero(T[1,0]-J[2]/S.sqrt(5)) and zero(T[2,0]-Q[2,2]/2) and
          zero(T[3,0]-S.sqrt(5)*O[2,2,2]/3))
    check('spin.helicity_selection', all(all(T[K,m][a,b]==0 for a,b in product(range(4),repeat=2)
          if HELICITIES[a]-HELICITIES[b] != m) for K,m in T))
    check('spin.spherical_generators', all(zero(J[2]*A-A*J[2]-m*A) for (K,m),A in T.items()))
    jp=J[0]+S.I*J[1];jm=J[0]-S.I*J[1]
    check('spin.spherical_ladder',all(zero(jp*T[K,m]-T[K,m]*jp-
          S.sqrt((K-m)*(K+m+1))*T[K,m+1]) for K in range(4) for m in range(-K,K)) and
          all(zero(jm*T[K,m]-T[K,m]*jm-
          S.sqrt((K+m)*(K-m+1))*T[K,m-1]) for K in range(4) for m in range(-K+1,K+1)))
    check('spin.cartesian_density_inverse', all(zero(reconstruct(X,J,Q,O)-X)
          for X in hermitian_basis()))
    check('spin.spherical_density_inverse', all(zero(sum((S.trace(A.H*X)*A for A in T.values()),S.zeros(4))-X)
          for X in hermitian_basis()))
    mats=[S.Matrix([[1,S.I,0,0],[0,2,1,0],[1,0,S.I,1],[0,0,0,1]]),
          S.Matrix([[1,0,0,0],[0,1,S.I,0],[0,0,1,0],[0,0,0,2]])]
    check('spin.positive_complex_states', all(zero(reconstruct((A*A.H)/S.trace(A*A.H),J,Q,O)
          -(A*A.H)/S.trace(A*A.H)) for A in mats))
    # Full Cartesian contractions, including repeated-index multiplicities.
    check('spin.trace_metrics',
          all(zero(S.trace(J[a]*J[b])-5*int(a==b)) for a,b in product(range(3),repeat=2)) and
          all(zero(S.trace(Q[a,b]*Q[c,d])-3*(int(a==c and b==d)+int(a==d and b==c)-R(2,3)*int(a==b and c==d)))
              for a,b,c,d in product(range(3),repeat=4)) and
          zero(S.trace(O[2,2,2]**2)-R(9,5)))
    beta=S.symbols('beta', nonnegative=True, real=True)
    X3=S.sqrt(5)/3*(J[2]**3-R(41,20)*J[2])
    rho=I4/4+beta*X3
    def pi3(a,b,c,d,e,f):
        left=(a,b,c);right=(d,e,f)
        leading=R(sum(S.prod(int(left[i]==right[perm[i]]) for i in range(3))
                    for perm in set(permutations(range(3)))),6)
        subtract=sum(int(left[i]==left[j])*sum(int(left[k]==right[r])*int(right[t]==right[u])
                     for r,t,u in ((0,1,2),(1,0,2),(2,0,1)))
                     for i,j,k in ((0,1,2),(0,2,1),(1,2,0)))/S.Integer(15)
        return leading-subtract
    check('spin.full_octupole_metric',all(zero(S.trace(O[a,b,c]*O[d,e,f])-R(9,2)*pi3(a,b,c,d,e,f))
          for a,b,c,d,e,f in product(range(3),repeat=6)))
    check('spin.cross_rank_orthogonality',all(zero(S.trace(A*B))
          for collection_a,collection_b in ((J,tuple(Q.values())),(J,tuple(O.values())),
                                          (tuple(Q.values()),tuple(O.values())))
          for A in collection_a for B in collection_b))
    check('spin.cartesian_permutation',all(Q[a,b]==Q[b,a] for a,b in Q) and
          all(O[a,b,c]==O[perm] for (a,b,c) in O for perm in set(permutations((a,b,c)))))
    check('spin.metric_inverse_weights',R(1,5)==1/S.trace(J[2]**2) and
          R(1,6)==1/(S.trace(Q[2,2]**2)/(R(2,3))) and
          R(2,9)==1/(S.trace(O[2,2,2]**2)/(R(2,5))))
    check('spin.octupole_spectrum', sorted([S.simplify(v) for v in X3.eigenvals()],key=str)==
          sorted([S.sqrt(5)*v/10 for v in (1,-3,3,-1)],key=str))
    edge=S.sqrt(5)/6
    check('spin.preparation_positivity', all(v>=0 for v in (rho.subs(beta,edge)).eigenvals()) and
          all(v>0 for v in (rho.subs(beta,edge/2)).eigenvals()) and
          all(v==R(1,4) for v in rho.subs(beta,0).eigenvals()) and
          any(v<0 for v in (rho.subs(beta,edge+R(1,10))).eigenvals()))
    partner=I4/4-beta*X3
    check('spin.pure_octupole', all(zero(S.trace(state*op)) for state in (rho,partner)
          for op in (*J,*Q.values())) and
          zero(S.trace(rho*O[2,2,2])-3*beta/S.sqrt(5)) and
          zero(S.trace(partner*O[2,2,2])+3*beta/S.sqrt(5)) and
          zero(S.trace(rho-partner)))
    return checks
