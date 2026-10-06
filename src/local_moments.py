"""Finite rotational content and positive-x antiquark bookkeeping.

This does not identify an ordinary integral with a renormalized local
operator matrix element without the source's convergence assumptions.
"""
from __future__ import annotations
from collections import Counter
import sympy as s


def tensor_weights(n, bilinear):
    if n<1 or bilinear not in ('vector','axial','tensor'):
        raise ValueError('local moment domain')
    # SO(3) restriction of a Lorentz vector is 0 + 1; the six tensor
    # components are two copies of 1 (electric and magnetic vectors).
    weights=Counter({-1:1,0:2,1:1}) if bilinear!='tensor' else Counter({-1:2,0:2,1:2})
    for _ in range(n-1):
        following=Counter()
        for m0,multiplicity in weights.items():
            for delta,count in ((-1,1),(0,2),(1,1)):
                following[m0+delta]+=multiplicity*count
        weights=following
    return {m:c for m,c in weights.items() if c}

def spin_multiplicities(weights):
    top=max(weights)
    return {j:weights.get(j,0)-weights.get(j+1,0) for j in range(top+1)
            if weights.get(j,0)-weights.get(j+1,0)}

def coupled_representations(n,bilinear):
    # Independent Clebsch-Gordan spin coupling of 0+1 derivative slots.
    contents={0:1,1:1} if bilinear!='tensor' else {1:2}
    for _ in range(n-1):
        following=Counter()
        for j,multiplicity in contents.items():
            following[j]+=multiplicity      # time derivative
            for k in range(abs(j-1),j+2):
                following[k]+=multiplicity  # spatial derivative
        contents=dict(following)
    return contents

def moment_sign(channel,n):
    if n<1 or channel not in ('vector','axial','transversity'):
        raise ValueError('moment channel/domain')
    return (-1)**(n-(channel=='axial'))

def bilinear_charge_signs():
    from process_dirac import basis
    g,g5,_,_,sigma=basis()
    C=s.I*g[2]*g[0]
    structures={'vector':g[1],'axial':g[1]*g5,
                'transversity':s.I*sigma(g[1],g[2])*g5}
    out={}
    for channel,G in structures.items():
        transformed=s.simplify(C.inv()*G.T*C)
        eta=1 if transformed==G else -1 if transformed==-G else None
        if eta is None:raise AssertionError(f'charge conjugation {channel}')
        out[channel]=eta
    if out!={'vector':-1,'axial':1,'transversity':-1}:
        raise AssertionError('source bilinear C parity')
    return out

def check_rank_and_charge():
    rows=[]
    charge=bilinear_charge_signs()
    for n in (1,2,3,4):
        for kind in ('vector','axial','tensor'):
            weights=tensor_weights(n,kind)
            reps=spin_multiplicities(weights)
            other=coupled_representations(n,kind)
            if reps!=other or max(reps)!=n:
                raise AssertionError(f'rotation restriction {kind} N={n}')
            rows.append({'N':n,'bilinear':kind,'multiplicities':reps,
                         'target_rank_sectors':[k for k in range(4) if k in reps],
                         'ambient_max_rank':max(reps)})
    q,anti=s.symbols('q qbar',real=True)
    parity=[]
    for n in (1,2,3,4):
        for channel in ('vector','axial','transversity'):
            sign=moment_sign(channel,n)
            parity.append({'N':n,'channel':channel,'antiquark_sign':sign,
                           'combination':str(q+sign*anti)})
            if sign!=charge[channel]*(-1)**(n-1):
                raise AssertionError('bilinear C parity and derivative count')
    expected={('vector',1):-1,('axial',1):1,('transversity',1):-1,
              ('axial',2):-1,('transversity',2):1}
    if any(moment_sign(ch,n)!=sign for (ch,n),sign in expected.items()):
        raise AssertionError('source local antiquark combinations')
    induction=[]
    for n in range(1,8):
        for kind in ('vector','axial','tensor'):
            top=max(spin_multiplicities(tensor_weights(n,kind)))
            if top!=n:raise AssertionError(f'rotation induction N={n} {kind}')
        induction.append(n)
    return {'rotational_rows':rows,'parity_rows':parity,
            'bilinear_charge_signs':charge,'induction_checked_N':induction,
            'induction_rule':'(0+1) tensor power adds at most one unit of rotational rank per derivative',
            'first_moment_absences':['f_1LL vector K2','g_1LLL axial K3',
                                     'h_1LLT tensor K3'],
            'second_moment_absences':['g_1LLL axial K3','h_1LLT tensor K3'],
            'first_allowed_reduced_octupole_N':3,
            'projection_statement':'trace and Young projection cannot increase ambient maximum rank',
            'source_labels':['eq:localoperators','eq:momentselection',
                             'eq:localmoments','eq:secondmoments']}

def nuclear_bookkeeping():
    x,y,k,p=s.symbols('x y k p',positive=True)
    z=s.symbols('x_N',positive=True)
    n=s.symbols('N',integer=True,positive=True)
    shift=k-x*p/y
    if s.simplify(shift-(k-z*p)).subs(x,y*z)!=0:
        raise AssertionError('nuclear transverse shift')
    factor=s.simplify((y*z)**(n-1)*y/y)
    if s.simplify(factor-y**(n-1)*z**(n-1))!=0:
        raise AssertionError('nuclear Mellin Jacobian')
    # Counterexample: counts fix integrals of nucleon numbers, not their
    # first y moments.  Delta support at y=1/14 gives 3+4 nucleons but
    # only half of the nuclear plus momentum in the nucleon component.
    number=3+4; momentum=3*s.Rational(1,14)+4*s.Rational(1,14)
    if number!=7 or momentum==1:raise AssertionError('nuclear count inference')
    if nuclear_kinematics(s.Rational(1,4),s.Rational(1,2),2,1)!=(
       s.Rational(1,2),s.Rational(3,2)):
        raise AssertionError('nuclear fraction/shift fixture')
    for bad in ((1,0),(s.Rational(3,4),s.Rational(1,2)),(-1,s.Rational(1,2))):
        try:nuclear_kinematics(*bad,2,1)
        except ValueError:pass
        else:raise AssertionError('nuclear support rejection')
    return {'x_N':'x/y','k_relative':'k_T-(x/y)*p_NT',
            'positive_x_support':'0 < x <= y <= 1; y=0 excluded',
            'Mellin_factor':'y**(N-1)',
            'nucleon_counts':[3,4],
            'count_only_momentum_counterexample':str(momentum),
            'conditional_total_tensor_witness':{'quark_flavor_1':'1',
                                                'gluon':'-1','total':'0'},
            'total_tensor_constraint':'sum of all partonic contributions conditional on total QCD momentum generator; no flavorwise zero',
            'source_labels':['eq:nuclearshift','eq:nuclearconv','eq:nuclearmomentum']}

def nuclear_kinematics(x,y,k,p):
    x,y,k,p=map(s.sympify,(x,y,k,p))
    if not all(v.is_number for v in (x,y,k,p)) or not (0<x<=y<=1):
        raise ValueError('positive-x impulse support requires 0<x<=y<=1')
    return s.simplify(x/y),s.simplify(k-x*p/y)
