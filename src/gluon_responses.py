"""All-term gluon hard responses by Cartesian trace and helicity routes."""
from __future__ import annotations
import sympy as s
from correlator_foundations import catalogue,cartesian_value,helicity_value
from gluon_stokes import density,reconstruct
from gluon_response_fixture import ROWS


def _factor(row):
    K,m,ch,n,c,p,q,parity=row
    phi,psi,beta=s.symbols('phi psi beta',real=True)
    t,M=s.symbols('t M',positive=True)
    if ch=='h': angle=p*phi+q*psi-2*beta
    else: angle=p*phi+q*psi
    weight=s.Rational(c.numerator,c.denominator) if hasattr(c,'denominator') else s.Integer(c)
    return weight*t**n*({'one':1,'sin':s.sin(angle),'cos':s.cos(angle)}[parity])


def response(label,route='cartesian'):
    phi,psi,beta=s.symbols('phi psi beta',real=True)
    t,M=s.symbols('t M',positive=True)
    fun=cartesian_value if route=='cartesian' else helicity_value
    vals=[fun(label,part,M*t*s.cos(phi),M*t*s.sin(phi),M)
          for part in range(1 if label.m==0 else 2)]
    weight=[s.Integer(1)] if label.m==0 else [s.cos(label.m*psi),s.sin(label.m*psi)]
    F,G,C,S=(sum(weight[i]*vals[i][j] for i in range(len(weight))) for j in range(4))
    if route=='cartesian':
        bu,bg,bc,bs=(int(label.channel==x) for x in ('f','g','h','none'))
        # For the linear-polarization channel, the analyzer direction is
        # an independent angle; for scalar channels only one Stokes slot acts.
        if label.channel=='h':bc,bs=s.cos(2*beta),s.sin(2*beta)
        return s.expand(s.trace(reconstruct(bu,bg,bc,bs)*density(F,G,C,S)))
    if label.channel=='f':return s.expand(F)
    if label.channel=='g':return s.expand(G)
    return s.expand(C*s.cos(2*beta)+S*s.sin(2*beta))


def verify_rows():
    labels={(z.K,z.m,z.channel,z.n):z for z in catalogue('gluon')}
    assert len(labels)==len(ROWS)==32
    out={}
    for row in ROWS:
        key=row[:4]
        label=labels[key]
        expected=_factor(row)
        for route in ('cartesian','helicity'):
            actual=response(label,route)
            # Exponential Laurent form retains the complete high harmonics;
            # generic trigsimp can leave equivalent rank-five forms unequal.
            residual=s.simplify(s.expand((actual-expected).rewrite(s.exp),power_exp=True))
            if residual!=0:
                raise AssertionError(f'{label.id()} {route}: {residual}')
        out[label.id()]={'K':label.K,'m':label.m,'n':label.n,
                         'channel':label.channel,'routes':['cartesian.trace','helicity.scalar'],
                         'mass_domain':'M_A>0; t=|k_T|/M_A',
                         'source_labels':['eq:Ffull','eq:Gfull','tab:gkernels','eq:gkernelrule']}
    return out
