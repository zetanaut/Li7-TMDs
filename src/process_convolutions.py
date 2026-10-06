"""Synthetic transverse Gaussian integrals; no nuclear TMD prediction.

Three routes: direct Cartesian hard-kernel quadrature, curated harmonic
quadrature, and exact conditional Gaussian moments. Widths have momentum²
units. sigma=-1 for SIDIS p, +1 for DY k_B; this is kinematic only.
"""
from __future__ import annotations
import math
import numpy as np
import sympy as s
from numpy.polynomial.hermite import hermgauss
from process_responses import (a,b,c,d,H,T,E,derived_expression,
                               expected_expression)

R=s.Rational
FIXTURE={'width_a':R(7,5),'width_b':R(11,10),'recoil':R(9,10),
         'mass_a':R(14,5),'mass_b':R(7,5),'z':R(3,5),
         'H':(3+4*s.I)/5,'T':(5+12*s.I)/13,'E':(8+15*s.I)/17}


def sigma_p(process):
    if process=='SIDIS':return -1
    if process=='DY':return 1
    raise ValueError(process)


def recoil_argument(process,fixture=FIXTURE):
    return fixture['recoil']/fixture['z'] if process=='SIDIS' else fixture['recoil']


def gaussian_moment_1d(degree,mu,var):
    # Generating function exp(mu*t+var*t²/2), independently of delta solver.
    return sum(s.factorial(degree)/(s.factorial(degree-2*j)*s.factorial(j)*2**j)
               *mu**(degree-2*j)*var**j for j in range(degree//2+1))


def exact_integral(process,row,fixture=FIXTURE):
    sig=sigma_p(process);wa=fixture['width_a'];wb=fixture['width_b']
    r=recoil_argument(process,fixture);ma=fixture['mass_a'];mb=fixture['mass_b']
    # Complete the square directly in the two original radial exponents.
    W=wa+wb;mu=wa*r/W;var=wa*wb/(2*W)
    assert s.simplify(s.diff((s.Symbol('x')**2/wa+(r-s.Symbol('x'))**2/wb),s.Symbol('x')).subs(s.Symbol('x'),mu))==0
    values={H:fixture['H'],T:fixture['T'],E:fixture['E']}
    expr=expected_expression(process,row).subs(values)
    x,y=s.symbols('x y',real=True)
    poly=s.Poly(s.expand(expr.subs({a:x/ma,b:y/ma,
                                    c:s.Integer(sig)*(r-x)/mb,d:-s.Integer(sig)*y/mb})),x,y)
    expectation=sum(coefficient*gaussian_moment_1d(px,mu,var)*
                    gaussian_moment_1d(py,0,var)
                    for (px,py),coefficient in poly.terms())
    norm=s.pi*wa*wb/W*s.exp(-r*r/W)
    return s.simplify(norm*expectation)


def quadrature(process,row,order=8,route='cartesian',fixture=FIXTURE):
    sig=sigma_p(process);wa=float(fixture['width_a']);wb=float(fixture['width_b'])
    r=float(recoil_argument(process,fixture));ma=float(fixture['mass_a']);mb=float(fixture['mass_b'])
    mu=wa*r/(wa+wb);scale=math.sqrt(wa*wb/(wa+wb))
    nodes,weights=hermgauss(order)
    # Cartesian route retains reflection-odd terms in the actual integrand.
    expr=(derived_expression(process,row,reflect=False) if route=='cartesian'
          else expected_expression(process,row))
    fn=s.lambdify((a,b,c,d,H,T,E),expr,'numpy')
    h=complex(fixture['H']);t=complex(fixture['T']);e=complex(fixture['E'])
    total=0j
    for ix,x in enumerate(nodes):
        kx=mu+scale*x;px=sig*(r-kx)
        for iy,y in enumerate(nodes):
            ky=scale*y;py=-sig*ky
            total+=weights[ix]*weights[iy]*complex(fn(kx/ma,ky/ma,px/mb,py/mb,h,t,e))
    norm=wa*wb/(wa+wb)*math.exp(-r*r/(wa+wb))
    return float((norm*total).real)


def check_row_integral(process,row,orders=(2,4,8,12)):
    exact=float(s.N(exact_integral(process,row),16))
    traces={route:[quadrature(process,row,n,route) for n in orders]
            for route in ('cartesian','harmonic')}
    scale=max(1,abs(exact))
    residuals={route:[abs(value-exact) for value in vals] for route,vals in traces.items()}
    if any(residuals[route][-1]>2e-11*scale for route in traces):
        raise AssertionError(f'{process} convolution mismatch {row}: {exact} {traces}')
    if row[3]==4 and any(residuals[route][0]<100*residuals[route][-1]
                         for route in traces):
        raise AssertionError(f'{process} highest-rank quadrature did not converge')
    return {'exact_gaussian':exact,'orders':list(orders),'quadrature':traces,
            'absolute_residuals':residuals,'relative_scale':scale,
            'mass_a':float(FIXTURE['mass_a']),'mass_b':float(FIXTURE['mass_b']),
            'width_a':float(FIXTURE['width_a']),'width_b':float(FIXTURE['width_b']),
            'z':float(FIXTURE['z']),'sigma_p':sigma_p(process)}


def momentum_sign_control():
    wa=FIXTURE['width_a'];wb=FIXTURE['width_b'];r=FIXTURE['recoil']
    W=wa+wb
    mean_p_plus=wb*r/W
    mean_p_minus=-wb*r/W
    unweighted=s.pi*wa*wb/W*s.exp(-r*r/W)
    assert mean_p_plus==-mean_p_minus and mean_p_plus!=0
    assert unweighted==s.pi*wa*wb/W*s.exp(-r*r/W)
    return {'unweighted_same':True,'odd_p_plus':float(mean_p_plus),
            'odd_p_minus':float(mean_p_minus),'odd_residual':float(2*mean_p_plus)}
