"""Rank-zero to rank-five Fourier/STF identities on convergent fixtures."""
from __future__ import annotations
import math
import numpy as np
import mpmath as mp
import sympy as s
from transverse_foundations import cartesian_tensor

R=s.Rational
bx,by=s.symbols('b_x b_y',real=True)
lam,mass=s.symbols('Lambda M_A',positive=True)
kx,ky=s.symbols('k_x k_y',real=True)

def normalized_component(rank,branch=1):
    if rank==0:return s.Integer(1)
    if rank<0 or rank>5 or branch not in (-1,1):raise ValueError('tensor rank/branch')
    return 2**(1-rank)*(kx+branch*s.I*ky)**rank/mass**rank

def derivative_component(rank,indices):
    """Apply (-i)^n STF(derivatives)/M^n to the Gaussian scalar."""
    scalar=s.exp(-lam**2*(bx**2+by**2)/4)
    tensor=cartesian_tensor(rank,kx,ky,mass)[indices]
    polynomial=s.Poly(tensor,kx,ky)
    return s.expand(sum(co*(-s.I)**rank*s.diff(scalar,bx,a,by,b)
                        for (a,b),co in polynomial.terms()))

def exact_gaussian():
    scalar=s.exp(-lam**2*(bx**2+by**2)/4)
    checked=0
    for rank in range(6):
        source=cartesian_tensor(rank,kx,ky,mass)
        target=cartesian_tensor(rank,bx,by,mass)
        for index in source:
            got=derivative_component(rank,index)
            expected=(s.I*lam**2/2)**rank*target[index]*scalar
            if s.simplify(got-expected)!=0:
                raise AssertionError(f'Fourier Cartesian derivative n={rank} index={index}')
            checked+=1
            if rank>0 and s.simplify(got.subs({bx:0,by:0}))!=0:
                raise AssertionError(f'Fourier b=0 Cartesian limit n={rank}')
        if rank:
            plus=s.expand(sum(derivative_component(rank,(0,)*rank) for _ in (0,))+
                          s.I*derivative_component(rank,(1,)+(0,)*(rank-1)))
            # The named complex STF component is T_x...x + i T_yx...x.
            expected=2**(1-rank)*(s.I*lam**2/2)**rank*(bx+s.I*by)**rank*scalar/mass**rank
            if s.simplify(plus-expected)!=0:
                raise AssertionError(f'Fourier complex STF normalization n={rank}')
    # Inverse transform of the scalar Gaussian.  The two 1D integrations
    # supply (2*pi)^-2 exactly; b=0 yields the original normalization.
    norm=s.simplify((s.pi/(lam**2/4))/(2*s.pi)**2)
    if norm!=1/(s.pi*lam**2):raise AssertionError('inverse Fourier normalization')
    m0=s.symbols('M_0',positive=True)
    for rank in range(1,6):
        a=cartesian_tensor(rank,bx,by,mass)
        b=cartesian_tensor(rank,bx,by,m0)
        if any(s.simplify(a[key]-(m0/mass)**rank*b[key])!=0 for key in a):
            raise AssertionError(f'Fourier mass conversion n={rank}')
    radial=s.symbols('k',positive=True)
    quartic=s.integrate(2*s.pi*radial*s.exp(-(radial/lam)**4),
                        (radial,0,s.oo))
    if s.simplify(quartic-s.pi**s.Rational(3,2)*lam**2/2)!=0:
        raise AssertionError('quartic radial normalization')
    return {'Cartesian_components':checked,'ranks':list(range(6)),
            'scalar_transform':'exp(-Lambda**2*b_T**2/4)',
            'inverse_normalization':str(norm),
            'b_zero':'scalar=1; tensor ranks 1..5 vanish without dividing by b',
            'mass_rescaling':'c(M_A)=(M_0/M_A)^n c(M_0)',
            'quartic_integral':str(quartic),
            'domain':'Lambda,M_A>0; all polynomial-weighted Gaussian integrals finite',
            'source_labels':['eq:Kn','eq:Kcomplex','eq:Fourierderivative','eq:Besseltransform']}

def _radial(k,scale,kind):
    if kind=='gaussian':return mp.exp(-(k/scale)**2)/(mp.pi*scale**2)
    if kind=='quartic':return 2*mp.exp(-(k/scale)**4)/(mp.pi**1.5*scale**2)
    raise ValueError(kind)

def cartesian_quadrature(rank,branch,bvec,scale,M,kind,order,cutoff):
    points,weights=np.polynomial.legendre.leggauss(order)
    grid=cutoff*points; weights=cutoff*weights
    xx,yy=np.meshgrid(grid,grid,indexing='ij')
    kk=(xx*xx+yy*yy)/scale**2
    radial=(np.exp(-kk)/(math.pi*scale**2) if kind=='gaussian' else
            2*np.exp(-kk*kk)/(math.pi**1.5*scale**2))
    value=np.exp(1j*(xx*bvec[0]+yy*bvec[1]))*((xx+branch*1j*yy)/M)**rank*radial
    return complex(np.sum(weights[:,None]*weights[None,:]*value))

def bessel_quadrature(rank,branch,bvec,scale,M,kind,dps=40):
    with mp.workdps(dps):
        bx0,by0=map(mp.mpf,map(str,bvec)); radius=mp.sqrt(bx0**2+by0**2)
        phase=((bx0+branch*1j*by0)/radius)**rank if radius else 1
        # Explicit k dk measure.  The finite upper bound has a checked tail.
        upper=mp.mpf('8')*mp.mpf(str(scale)) if kind=='gaussian' else 4*mp.mpf(str(scale))
        integrand=lambda k:k*(k/mp.mpf(str(M)))**rank*mp.besselj(rank,k*radius)*_radial(k,mp.mpf(str(scale)),kind)
        integral=mp.quad(integrand,[0,upper/2,upper])
        return complex(2*mp.pi*(1j)**rank*phase*integral)

def inverse_gaussian_harmonic(rank,k,scale,M,dps=45):
    """Independent inverse Hankel integral with (2*pi)^-2 normalization."""
    with mp.workdps(dps):
        scale=mp.mpf(str(scale));mass0=mp.mpf(str(M));k=mp.mpf(str(k))
        a=scale**2/4
        radial=lambda b:(1j*scale**2/(2*mass0))**rank*b**rank*mp.exp(-a*b*b)
        upper=16/scale
        got=(-1j)**rank/(2*mp.pi)*mp.quad(lambda b:b*mp.besselj(rank,k*b)*radial(b),
                                          [0,upper/2,upper])
        expected=(k/mass0)**rank*mp.exp(-(k/scale)**2)/(mp.pi*scale**2)
        return complex(got),float(expected)

def numeric_comparison():
    cases=[]
    for kind,rank_set,scale,M,bvec,cut in (
       ('gaussian',range(6),1.7,2.3,(0.45,0.63),10.2),
       ('quartic',(0,1,3,5),1.35,2.2,(0.32,0.47),4.05)):
        for rank in rank_set:
            for branch in (-1,1):
                radial=bessel_quadrature(rank,branch,bvec,scale,M,kind)
                coarse=cartesian_quadrature(rank,branch,bvec,scale,M,kind,40,cut)
                fine=cartesian_quadrature(rank,branch,bvec,scale,M,kind,72,cut)
                scale_for_error=max(abs(radial),1e-8)
                error=abs(fine-radial)
                if error>5e-8*scale_for_error:
                    raise AssertionError(f'Fourier routes {kind} n={rank} branch={branch}: {error}')
                if abs(fine-radial)>abs(coarse-radial)+2e-10*scale_for_error:
                    raise AssertionError('Fourier quadrature did not converge')
                analytic=None
                if kind=='gaussian':
                    z=bvec[0]+branch*1j*bvec[1]
                    analytic=(1j*scale**2/(2*M))**rank*z**rank*math.exp(-scale**2*sum(q*q for q in bvec)/4)
                    if abs(radial-analytic)>2e-11*scale_for_error:
                        raise AssertionError(f'Gaussian Bessel expectation n={rank}')
                cases.append({'kind':kind,'rank':rank,'branch':branch,
                              'b':list(bvec),'Lambda':scale,'M_A':M,
                              'cartesian_orders':[40,72],
                              'cartesian_values':[[v.real,v.imag] for v in (coarse,fine)],
                              'bessel':[radial.real,radial.imag],
                              'analytic':None if analytic is None else [analytic.real,analytic.imag],
                              'coarse_absolute_error':abs(coarse-radial),
                              'absolute_error':error,'relative_scale':scale_for_error,
                              'cutoff':cut})
    # For an odd real momentum component the transform is imaginary and obeys
    # F*(b)=F(-b), rather than being pointwise real.
    odd=(bessel_quadrature(1,1,(.45,.63),1.7,2.3,'gaussian')+
         bessel_quadrature(1,-1,(.45,.63),1.7,2.3,'gaussian'))/2
    opposite=(bessel_quadrature(1,1,(-.45,-.63),1.7,2.3,'gaussian')+
              bessel_quadrature(1,-1,(-.45,-.63),1.7,2.3,'gaussian'))/2
    if abs(odd.conjugate()-opposite)>2e-12 or abs(odd.imag)<1e-3:
        raise AssertionError('odd Fourier conjugation')
    inverse=[]
    for rank in (0,1,5):
        for k in (.65,1.1):
            got,expected=inverse_gaussian_harmonic(rank,k,1.7,2.3)
            if abs(got-expected)>2e-12*max(abs(expected),1e-8):
                raise AssertionError(f'inverse Gaussian Hankel n={rank} k={k}')
            inverse.append({'rank':rank,'k':k,'value':[got.real,got.imag],
                            'expected':expected,'absolute_error':abs(got-expected)})
    return {'cases':cases,'case_count':len(cases),
            'inverse_cases':inverse,
            'max_absolute_error':max(c['absolute_error'] for c in cases),
            'max_scaled_error':max(c['absolute_error']/c['relative_scale'] for c in cases),
            'quartic_normalization':'2/(pi**(3/2)*Lambda**2)',
            'measure':'d2k and k dk; positive Fourier exponential',
            'dimensions':'c: momentum^-2; transform: dimensionless; b: momentum^-1; M_A: momentum'}
