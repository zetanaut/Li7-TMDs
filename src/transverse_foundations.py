"""Two-dimensional STF tensors from harmonic projection and helicity.

Source: eq:Kn, Kcomplex, dual, K45, dyadic_not_STF, massconversion.
The Cartesian route uses a harmonic polynomial in auxiliary variables;
the helicity route uses (kx+i ky)^n independently.
"""
from __future__ import annotations
from itertools import product
import sympy as S

R=S.Rational
X,Y=S.symbols('kx ky', real=True)
M=S.symbols('M_A', positive=True)
U,V=S.symbols('u v', real=True)
E=S.Matrix([[0,1],[-1,0]])


def harmonic_polynomial(n,kx=X,ky=Y,mass=M):
    if n==0: return S.Integer(1)
    dot=kx*U+ky*V
    radius=U*U+V*V
    k2=kx*kx+ky*ky
    coefficient=S.Integer(1)
    out=dot**n
    for j in range(1,n//2+1):
        coefficient *= -R((n-2*j+2)*(n-2*j+1),4*j*(n-j))
        out += coefficient*k2**j*radius**j*dot**(n-2*j)
    return S.expand(out/mass**n)


def cartesian_tensor(n,kx=X,ky=Y,mass=M):
    if n==0:return {():S.Integer(1)}
    p=S.Poly(harmonic_polynomial(n,kx,ky,mass),U,V)
    return {idx:S.expand(p.coeff_monomial(U**(n-sum(idx))*V**sum(idx))/S.binomial(n,sum(idx)))
            for idx in product(range(2),repeat=n)}


def helicity_tensor(n,kx=X,ky=Y,mass=M):
    if n==0:return {():S.Integer(1)}
    plus=S.expand((kx+S.I*ky)**n)/(2**(n-1)*mass**n)
    real=S.expand_complex(plus).as_real_imag()[0]
    imag=S.expand_complex(plus).as_real_imag()[1]
    return {idx:S.expand((-1)**(sum(idx)//2)*(real if sum(idx)%2==0 else imag))
            for idx in product(range(2),repeat=n)}


def dual(t):
    if () in t: raise ValueError('rank-zero tensor has no dual')
    return {idx:sum(E[idx[0],j]*t[(j,)+idx[1:]] for j in range(2)) for idx in t}


def stf2(A):return S.expand((A+A.T)/2-S.eye(2)*S.trace(A)/2)


def dyadic_regression(kx=X,ky=Y):
    a,b=S.symbols('S_xx S_xy',real=True)
    target=S.Matrix([[a,b],[b,-a]])
    k=S.Matrix([kx,ky]);k2=k.dot(k)
    rank2=k*k.T-S.eye(2)*k2/2
    # One Lorentz index of S is lowered with g_T=-delta. Braces are
    # the unnormalized sum of the two free-index permutations.
    contraction=-(target*rank2+(target*rank2).T)
    dyadic=-(target*(k*k.T)+(target*(k*k.T)).T)
    return target,contraction,dyadic,stf2(contraction),stf2(dyadic)


def run_checks():
    checks={}
    def check(name,condition,payload=None):
        if condition is not True:raise AssertionError(name)
        checks[name]=payload or {'domain':'exact polynomial over Q(kx,ky,M_A); M_A>0'}
    for n in range(6):
        c=cartesian_tensor(n);h=helicity_tensor(n)
        check(f'stf.rank_{n}.independent',all(S.expand(c[key]-h[key])==0 for key in c),
              {'rank':n,'components':len(c),'method':'harmonic projection versus complex helicity'})
        if n>=2:
            check(f'stf.rank_{n}.trace',all(S.expand(c[(0,0)+tail]+c[(1,1)+tail])==0
                  for tail in product(range(2),repeat=n-2)))
        if n>=1:
            d=dual(c)
            check(f'stf.rank_{n}.dual',all(S.expand(dual(d)[key]+c[key])==0 for key in c)
                  and (n<2 or all(S.expand(d[(0,0)+tail]+d[(1,1)+tail])==0
                  for tail in product(range(2),repeat=n-2))))
            normalized=cartesian_tensor(n,X/M,Y/M,S.Integer(1))
            check(f'stf.rank_{n}.legacy_mass',all(S.expand(c[key]-normalized[key])==0 for key in c))
            m0=S.symbols('M_0',positive=True)
            alt=cartesian_tensor(n,X,Y,m0)
            check(f'stf.rank_{n}.mass_conversion',all(S.cancel(c[key]-(m0/M)**n*alt[key])==0 for key in c),
                  {'rank':n,'relation':'c(M_0)=(M_0/M_A)^n c(M_A)'})
    tscale=S.symbols('lambda',real=True)
    for n in range(1,6):
        c=cartesian_tensor(n)
        rotated=cartesian_tensor(n,-Y,X,M)
        quarter=all(S.expand(rotated[idx]-sum(S.prod((0 if idx[j]==old[j] else (-1 if idx[j]==0 else 1))
                   for j in range(n))*c[old] for old in product(range(2),repeat=n)))==0
                   for idx in c)
        # Active +pi/2 rotation matrix [[0,-1],[1,0]].
        check(f'stf.rank_{n}.rotation',quarter)
        check(f'stf.rank_{n}.homogeneity',all(S.expand(cartesian_tensor(n,tscale*X,tscale*Y,M)[idx]-tscale**n*c[idx])==0 for idx in c))
        norm=S.expand(sum(value*value for value in c.values()))
        check(f'stf.rank_{n}.norm',S.cancel(norm-(X*X+Y*Y)**n/(2**(n-1)*M**(2*n)))==0)
        check(f'stf.rank_{n}.zero_momentum',all(value.subs({X:0,Y:0})==0 for value in c.values()))
    aa,bb,cc,dd=S.symbols('aa bb cc dd',real=True)
    def target_tensor(n,a,b):
        return {idx:(-1)**(sum(idx)//2)*(a if sum(idx)%2==0 else b)
                for idx in product(range(2),repeat=n)}
    for n in range(1,4):
        A=target_tensor(n,aa,bb);c=cartesian_tensor(n)
        direct=S.expand(sum(A[idx]*c[idx] for idx in A))
        sine=S.expand(-sum(dual(A)[idx]*c[idx] for idx in A))
        wave=S.expand((aa-S.I*bb)*(X+S.I*Y)**n/M**n)
        check(f'stf.rank_{n}.scalar_sine',S.expand(direct-S.re(wave))==0 and S.expand(sine-S.im(wave))==0)
    A2=S.Matrix([[aa,bb],[bb,-aa]]);B2=S.Matrix([[cc,dd],[dd,-cc]])
    check('stf.rank_two_products',stf2(A2*B2)==S.zeros(2) and stf2(E*A2*B2)==S.zeros(2))
    A3=target_tensor(3,aa,bb);B3=target_tensor(3,cc,dd)
    product3=S.Matrix(2,2,lambda i,j:sum(A3[i,u,v]*B3[j,u,v] for u,v in product(range(2),repeat=2)))
    check('stf.rank_three_product_dual',stf2(product3)==S.zeros(2) and stf2(E*product3)==S.zeros(2))
    target,contract,dyadic,cstf,dstf=dyadic_regression()
    residual=-(X*X+Y*Y)*target
    check('stf.dyadic_null',cstf==S.zeros(2))
    check('stf.dyadic_residual',all(S.expand(dstf[i,j]-residual[i,j])==0 for i,j in product(range(2),repeat=2)),
          {'residual':[[str(residual[i,j]) for j in range(2)] for i in range(2)]})
    check('stf.dyadic_pure_trace',all(S.expand(contract[i,j]-(S.trace(contract)/2 if i==j else 0))==0
          for i,j in product(range(2),repeat=2)))
    return checks
