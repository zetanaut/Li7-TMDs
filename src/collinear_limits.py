"""Finite angular selection and conditional straight-link collinear algebra.

The angular map is normalized by 1/(2*pi).  No radial TMD integral or
small-b matching coefficient is evaluated here.
"""
from __future__ import annotations
from math import factorial
import sympy as s
from correlator_foundations import (catalogue, cartesian_value, helicity_value,
                                    target_component_operators)
from process_reversal import derived_coefficient_sign
from spin_foundations import cartesian

R=s.Rational
X,Y=s.symbols('kx ky',real=True)
r,M=s.symbols('r M_A',positive=True)
EXPECTED_CANDIDATES={
 'quark':('quark.f.00[0]','quark.f.20[0]','quark.g.10[0]',
          'quark.g.30[0]','quark.h.11[0]','quark.h.21[0]','quark.h.31[0]'),
 'gluon':('gluon.f.00[0]','gluon.f.20[0]','gluon.g.10[0]',
          'gluon.g.30[0]','gluon.h.22[0]','gluon.h.32[0]'),
}
EXPECTED_COLLINEAR={
 'quark':tuple(x for x in EXPECTED_CANDIDATES['quark'] if x!='quark.h.21[0]'),
 'gluon':tuple(x for x in EXPECTED_CANDIDATES['gluon'] if x!='gluon.h.32[0]'),
}

def circle_average(expr, radius=r):
    """Cartesian monomial integral over a circle, divided by 2*pi."""
    poly=s.Poly(s.expand(expr),X,Y)
    out=0
    for (a,b),coefficient in poly.terms():
        if a%2 or b%2:continue
        u,v=a//2,b//2
        out+=coefficient*radius**(a+b)*R(factorial(2*u)*factorial(2*v),
                                         4**(u+v)*factorial(u)*factorial(v)*factorial(u+v))
    return s.cancel(out)

def selection(species):
    rows=[]
    for label in catalogue(species):
        parts=[]
        for part in range(1 if label.m==0 else 2):
            cart=cartesian_value(label,part,X,Y,M)
            harmonic=helicity_value(label,part,X,Y,M)
            lhs=tuple(circle_average(v) for v in cart)
            rhs=tuple(circle_average(v) for v in harmonic)
            if any(s.simplify(a-b)!=0 for a,b in zip(lhs,rhs)):
                raise AssertionError(f'angular route mismatch {label.id()}')
            parts.append(lhs)
        candidate=any(any(v!=0 for v in component) for component in parts)
        if candidate!=(label.n==0):
            raise AssertionError(f'orbital metadata mismatch {label.id()}')
        link_sign=derived_coefficient_sign(label)
        rows.append({'id':label.id(),'rank':label.n,'candidate':candidate,
                     'straight_link':candidate and link_sign==1,
                     'link_sign':link_sign,
                     'angular_components':[[str(v) for v in p] for p in parts]})
    candidates=tuple(x['id'] for x in rows if x['candidate'])
    survivors=tuple(x['id'] for x in rows if x['straight_link'])
    if set(candidates)!=set(EXPECTED_CANDIDATES[species]) or set(survivors)!=set(EXPECTED_COLLINEAR[species]):
        raise AssertionError(f'{species} collinear selection/source fixture')
    return {'species':species,'rows':rows,'candidates':list(candidates),
            'survivors':list(survivors),'angular_measure':'dphi/(2*pi)',
            'straight_link_assumption':'source operator PT symmetry at common scale/scheme'}

def operator_projection():
    ops=target_component_operators();J,Q,O,_,_=cartesian()
    fixture={'f.00':s.eye(4),'f.20':Q[2,2],
             'g.10':J[2],'g.30':O[2,2,2],
             'q.h.11':J[0],'q.h.31':O[2,2,0],
             'g.h.22':Q[0,0]-Q[1,1]}
    actual={'f.00':ops[0,0,0],'f.20':ops[2,0,0],
            'g.10':ops[1,0,0],'g.30':ops[3,0,0],
            'q.h.11':ops[1,1,0],'q.h.31':ops[3,1,0],
            'g.h.22':ops[2,2,0]}
    if any(s.simplify(actual[k]-v)!=s.zeros(4) for k,v in fixture.items()):
        raise AssertionError('collinear polarization normalization')
    f,fq,g,go,h,ho=s.symbols('f f_Q g g_O h h_O',real=True)
    Uo,Ui=f+fq,f-fq
    Go,Gi=R(3,2)*g+R(3,10)*go,g/2-R(9,10)*go
    if any(s.simplify(v)!=0 for v in (
       f-(Uo+Ui)/2,fq-(Uo-Ui)/2,
       g-(3*Go+Gi)/5,go-(Go-3*Gi)/3)):
        raise AssertionError('population inverse')
    jp=J[0]+s.I*J[1];op=O[2,2,0]+s.I*O[2,2,1]
    raising=s.expand(h*jp+ho*op)
    a,b=raising[0,1],raising[1,2]
    if s.simplify(a-s.sqrt(3)*(h+R(2,5)*ho)) or s.simplify(b-(2*h-R(6,5)*ho)):
        raise AssertionError('outer/middle coherence')
    return {'population':{'U_outer':str(Uo),'U_inner':str(Ui),
                          'G_outer':str(Go),'G_inner':str(Gi)},
            'coherence':{'outer':str(a),'middle':str(b)},
            'source_labels':['eq:qcolproj','eq:gcolproj','eq:populationresponses',
                             'eq:populationinverse','eq:coherenceblock']}

def operations_and_tail():
    radius,lam=s.symbols('R Lambda',positive=True)
    k2=s.symbols('k2',nonnegative=True)
    cutoff=s.integrate(2*s.pi*s.symbols('k',positive=True)/
                (s.pi*(s.symbols('k',positive=True)**2+lam**2)),
                (s.symbols('k',positive=True),0,radius))
    if s.simplify(cutoff-s.log(1+radius**2/lam**2))!=0:
        raise AssertionError('UV cutoff dependence')
    at_zero={species:[label.id() for label in catalogue(species) if label.n==0]
             for species in ('quark','gluon')}
    if any(set(at_zero[x])!=set(EXPECTED_CANDIDATES[x]) for x in at_zero):
        raise AssertionError('regular fixture k=0 selection')
    return {'cutoff_integral':str(cutoff),'unbounded_limit':'logarithmic',
            'regular_k_zero':at_zero,'distinct_operations':
            ['forward k_T=0 evaluation','fixed-radius angular average',
             'renormalized straight-link operator']}

def inclusive_born_bookkeeping():
    from process_dirac import pair_table
    table=pair_table('SIDIS')
    identity=s.eye(2); epsilon=s.Matrix([[0,1],[-1,0]])
    if table['F','D']!=identity or table['G','D']!=s.I*epsilon:
        raise AssertionError('inclusive Born transverse number/helicity current')
    x,e1,e2=s.symbols('x e_1 e_2',positive=True)
    q1,a1,q2,a2=s.symbols('q_1 qbar_1 q_2 qbar_2',real=True)
    F1=(e1**2*(q1+a1)+e2**2*(q2+a2))/2
    F2=2*x*F1
    if s.simplify(F2-x*(e1**2*(q1+a1)+e2**2*(q2+a2)))!=0:
        raise AssertionError('inclusive charge or target normalization')
    return {'number_current':str(table['F','D']),
            'helicity_current':str(table['G','D']),
            'F1_two_flavor':str(F1),'F2_over_2xF1':'1',
            'applies_to_target_ranks':[0,2],
            'helicity_target_ranks':[1,3],
            'SIDIS_to_inclusive_assumptions':['energy-weighted sum over hadrons',
                'fragmentation momentum sum rule','fragmentation variable Jacobians',
                'no unrestricted integration of low-P_hT approximation'],
            'source_label':'eq:inclusiveSF'}
