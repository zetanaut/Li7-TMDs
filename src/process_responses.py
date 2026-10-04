"""Response polynomials derived from Cartesian target covariants.

H,T,E are formal unit phases for recoil, target and negative-lepton
azimuths. Reflection in the recoil axis projects radial convolutions.
"""
from __future__ import annotations
import sympy as s
from correlator_foundations import Label, cartesian_value, catalogue
from response_fixtures import SIDIS_ROWS,DY_ROWS,row_id

a,b,c,d=s.symbols('a b c d',real=True)
H,T,E=s.symbols('H T E',nonzero=True)
M_A,M_B=s.symbols('M_A M_B',positive=True)
I=s.I


def conjugate_phase(expr):
    return s.conjugate(expr).subs({s.conjugate(H):1/H,
                                    s.conjugate(T):1/T,s.conjugate(E):1/E})


def even_reflection(expr):
    return s.expand((expr+expr.subs({b:-b,d:-d}))/2)


def target_parts(m):
    if m==0:return ((1,0),)
    z=(T/H)**m
    return (((z+1/z)/2,(z-1/z)/(2*I)),)


def cartesian_projection(label):
    re,im=target_parts(label.m)[0]
    zero=cartesian_value(label,0,M_A*a,M_A*b,M_A)
    if label.m:
        one=cartesian_value(label,1,M_A*a,M_A*b,M_A)
    else:one=(0,0,0,0)
    return tuple(s.expand(re*x+im*y) for x,y in zip(zero,one))


def derived_expression(process,row, *, reflect=True):
    K,m,channel,n,*_=row
    label=Label('quark',channel,K,m,n)
    f,g,tx,ty=cartesian_projection(label)
    project=even_reflection if reflect else s.expand
    if channel=='f':return project(f)
    if channel=='g':return project(g)*(1 if process=='SIDIS' else -1)
    tplus=tx+I*ty
    if process=='SIDIS':
        # Delta_T=-epsilon p H/M_h; contraction with normalized lepton tensor.
        complex_term=H**2*tplus*(c+I*d)
        expression=(complex_term-conjugate_phase(complex_term))/(2*I)
    else:
        # Minus-beam annihilation antiquark: bar T_+=+i v_B h_1^perp.
        complex_term=E**-2*tplus*I*(c+I*d)
        expression=(complex_term+conjugate_phase(complex_term))/2
    return project(expression)


def reviewed_weight(name):
    # Deliberately separate polynomial transcription of eq:weightpolys.
    return {
      'w0':s.Integer(1),'w1':a,'w2':a*a-b*b,'w3':a**3-3*a*b*b,
      'b0':a*c-b*d,'l1':c,'l2':a*c+b*d,
      'l3':(a*a-b*b)*c+2*a*b*d,
      'r1':((a*a-b*b)*c-2*a*b*d)/2,
      'r2':((a**3-3*a*b*b)*c-(3*a*a*b-b**3)*d)/2,
      'r3':((a**4-6*a*a*b*b+b**4)*c-(4*a**3*b-4*a*b**3)*d)/2,
    }[name]


def expected_expression(process,row):
    K,m,channel,n,weight,phase,sign=row
    z=H**phase[0]*T**phase[1]*E**phase[2]
    if process=='SIDIS':
        trig=('sin' if (channel=='f' and K%2) or (channel=='g' and K%2==0) or
              (channel=='h' and K%2) else 'cos')
        if channel in ('f','g') and m==0:trig='one'
    else:
        trig='sin' if channel in ('f','h') else 'cos'
        if channel=='g' and m==0:trig='one'
    angular={'one':s.Integer(1),'cos':(z+1/z)/2,
             'sin':(z-1/z)/(2*I)}[trig]
    return s.expand(sign*angular*reviewed_weight(weight))


def compare_row(process,row):
    actual=derived_expression(process,row)
    expected=expected_expression(process,row)
    residual=s.cancel(actual-expected)
    if residual!=0:
        raise AssertionError(f'{row_id(process,row)}: residual={residual}')
    return {'source_label':'tab:SF' if process=='SIDIS' else
            ('eq:DYnumberSF' if row[2]=='f' else 'eq:DYhelicitySF' if row[2]=='g' else 'eq:DYchiralSF'),
            'label':{'K':row[0],'m':row[1],'channel':row[2],'orbital_rank':row[3]},
            'weight':row[4],'phase_coefficients':list(row[5]),'signed_prefactor':row[6],
            'polynomial_coefficients':len(s.Poly(s.expand(reviewed_weight(row[4])),a,b,c,d).terms()),
            'routes':['Cartesian target covariant and hard contraction',
                      'independently curated Fourier polynomial fixture']}


def compare_all():
    known=set(catalogue('quark'))
    out={}
    for process,rows in (('SIDIS',SIDIS_ROWS),('DY',DY_ROWS)):
        for row in rows:
            label=Label('quark',row[2],row[0],row[1],row[3])
            assert label in known
            out[row_id(process,row)]=compare_row(process,row)
    return out
