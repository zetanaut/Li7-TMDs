"""Analytic orbital projectors versus independent Cartesian linear-map inverse.

Eight-point Fourier sums exactly equal the scalar angular integrals for
harmonics |m|<=3. The two-branch inverses implement eq:qinverse/ginverse.
"""
from __future__ import annotations
from itertools import product
import sympy as S
from correlator_foundations import catalogue,helicity_value,map_matrix

R=S.Rational
T=R(2,3)


def _complex_h(species,label,part,kx,ky):
    v=helicity_value(label,part,kx,ky,S.Integer(1))
    return v[2]+S.I*v[3]


def analytic_recover(species, coefficients, *, radial=T):
    radial=S.sympify(radial)
    if radial.is_zero is True:
        raise ValueError('analytic inverse undefined at k_T=0; use a nonzero radial point or a defined limit')
    if radial.is_zero is not False:
        raise ValueError('analytic inverse requires a proven nonzero radial point')
    labels=catalogue(species)
    if len(coefficients)!=len(labels):raise ValueError('coefficient count mismatch')
    recovered=[None]*len(labels)
    for i,label in enumerate(labels):
        if label.channel not in ('f','g'):continue
        matching=[j for j,other in enumerate(labels) if other.channel==label.channel and
                  (other.K,other.m)==(label.K,label.m)]
        assert matching==[i]
        m=label.m
        def response(phi):
            kx=radial*S.cos(phi);ky=radial*S.sin(phi)
            return sum(coefficients[j]*helicity_value(labels[j],0,kx,ky,S.Integer(1))[0 if label.channel=='f' else 1]
                       for j in matching)
        if m==0:
            recovered[i]=S.simplify(sum(response(S.pi*j/4) for j in range(8))/8)
        else:
            ordinary=(label.K%2==0)==(label.channel=='f')
            wave=S.cos if ordinary else S.sin
            recovered[i]=S.simplify(sum(wave(m*S.pi*j/4)*response(S.pi*j/4)
                              for j in range(8))*R(1,4)/radial**m)
    for K in range(4):
        for m in range(K+1):
            slots=[i for i,label in enumerate(labels) if label.channel=='h' and
                   (label.K,label.m)==(K,m)]
            if not slots:continue
            phase=S.I if (species=='quark' and K%2==0) or (species=='gluon' and K%2==1) else 1
            # physical plus = -i B for dual sectors, so B = i physical.
            def response(part):
                return S.expand(phase*sum(coefficients[i]*_complex_h(species,labels[i],part,radial,0)
                                          for i in slots))
            if m==0:
                pref=radial if species=='quark' else radial*radial/2
                recovered[slots[0]]=S.simplify(response(0)/pref)
            else:
                low,high=slots
                v1=response(0);vi=response(1)/S.I
                nlow=labels[low].n
                c=1 if species=='quark' else {1:R(1,2),2:S.Integer(1),3:S.Integer(1)}[m]
                recovered[low]=S.simplify((v1+vi)/(2*c*radial**nlow))
                recovered[high]=S.simplify((v1-vi)/(radial**labels[high].n*(1 if species=='quark' else R(1,2))))
    return recovered


def linear_recover(species,coefficients):
    labels,M=map_matrix(species)
    rows=M.T.rref()[1]
    selected=M.extract(rows,list(range(32)))
    response=M*S.Matrix(coefficients)
    return list(selected.inv()*response.extract(rows,[0]))


def run_checks():
    out={}
    for species in ('quark','gluon'):
        labels=catalogue(species)
        for j in range(32):
            unit=[S.Integer(int(i==j)) for i in range(32)]
            recovered=analytic_recover(species,unit)
            if any(S.simplify(x-y)!=0 for x,y in zip(unit,recovered)):
                raise AssertionError(f'{species}.analytic_projector unit={labels[j].id()}')
        out[f'{species}.analytic_projectors']={'unit_cases':32,'point':'k=(2/3,0), M_A=1',
           'scalar_method':'eight-point exact Fourier quadrature','helicity_flip':'published two-branch inverses'}
        mixed=[S.Rational(j+1,j+2) for j in range(32)]
        analytic=analytic_recover(species,mixed)
        linear=linear_recover(species,mixed)
        if any(S.simplify(a-b)!=0 or S.simplify(a-c)!=0 for a,b,c in zip(mixed,analytic,linear)):
            raise AssertionError(f'{species}.independent_linear_recovery')
        out[f'{species}.independent_linear_recovery']={'coefficients':32,'method':'analytic Fourier/branch versus exact selected-row inverse'}
    return out
