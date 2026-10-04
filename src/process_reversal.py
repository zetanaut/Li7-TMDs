"""Finite PT/link algebra, conditional on the source's field/link rule.

No Wilson-line matrix element, evolution kernel, or factorization theorem is
computed here. Momentum labels are fixed by the combined PT comparison.
"""
from __future__ import annotations
from dataclasses import dataclass,replace
import sympy as s
from correlator_foundations import catalogue,target_component_operators
from spin_foundations import HELICITIES


@dataclass(frozen=True)
class Link:
    species:str
    paths:tuple[str,...]
    color_word:str='single_trace'
    return_path_dagger:bool=True
    endpoint:str='standard'
    def __post_init__(self):
        if self.species not in ('quark','gluon'):raise ValueError('parton species')
        if len(self.paths)!=(1 if self.species=='quark' else 2):raise ValueError('path count')
        if any(path not in ('+','-') for path in self.paths):raise ValueError('path direction')
        if self.color_word!='single_trace':raise NotImplementedError('extra color loop requires reviewed definition')
        if not self.return_path_dagger or self.endpoint!='standard':
            raise NotImplementedError('unreviewed return path or endpoint')
    def reverse(self):
        return replace(self,paths=tuple('-' if p=='+' else '+' for p in self.paths))
    @property
    def class_name(self):
        return 'same_direction' if len(self.paths)==1 or self.paths[0]==self.paths[1] else 'mixed_direction'


@dataclass(frozen=True)
class Coefficient:
    species:str
    channel:str
    K:int
    m:int
    n:int
    parent:str
    flavor:str
    antiquark:bool
    link:Link
    mu:str
    zeta:str
    scheme:str
    def reverse(self):
        return replace(self,link=self.link.reverse())
    def comparable(self,other):
        same=replace(self,link=other.link)
        if same!=other or self.link.reverse()!=other.link:
            raise ValueError('mismatched species, flavor, parent, scale, scheme, or link')
        return True


def time_reversal_unitary():
    J=s.Rational(3,2)
    return s.Matrix(4,4,lambda r,c:(-1)**(J-HELICITIES[c])
                    if HELICITIES[r]==-HELICITIES[c] else 0)


def reverse_density(rho):
    U=time_reversal_unitary()
    return U*rho.conjugate()*U.H


def field_eta(species,channel):
    """Analytic input eq:PTprojection, including conjugation of explicit i."""
    if channel=='f':return 1
    if channel=='g':return -1
    if channel=='h':return -1 if species=='quark' else 1
    raise ValueError(channel)


def derived_coefficient_sign(label):
    ops=target_component_operators()
    U=time_reversal_unitary()
    for part in range(1 if label.m==0 else 2):
        op=ops[label.K,label.m,part]
        transformed=U*op.conjugate()*U.H
        target_sign=(-1)**label.K
        if s.simplify(transformed-target_sign*op)!=s.zeros(4):
            raise AssertionError(f'target PT: {label.id()} part={part}')
    return field_eta(label.species,label.channel)*target_sign


def compare_sign_tables():
    # Reviewed expectation is separately encoded as the two source tables.
    expected={'quark':{'f':(1,-1,1,-1),'g':(-1,1,-1,1),'h':(-1,1,-1,1)},
              'gluon':{'f':(1,-1,1,-1),'g':(-1,1,-1,1),'h':(1,-1,1,-1)}}
    out={}
    for species in ('quark','gluon'):
        rows=[]
        for label in catalogue(species):
            value=derived_coefficient_sign(label)
            if value!=expected[species][label.channel][label.K]:
                raise AssertionError(label.id())
            rows.append({'id':label.id(),'sign':int(value),'K':label.K,'m':label.m,
                         'orbital_rank':label.n,'channel':label.channel})
        octupole=[row for row in rows if row['K']==3]
        out[species]={'rows':rows,'octupole_odd':sum(x['sign']==-1 for x in octupole),
                      'octupole_even':sum(x['sign']==1 for x in octupole)}
    return out


def check_reversal():
    U=time_reversal_unitary()
    assert U*U.conjugate()==-s.eye(4)
    A=s.Matrix([[1,s.I,0,1],[s.I,2,1,0],[0,1,1,s.I],[1,0,s.I,1]])
    rho=A*A.H/s.trace(A*A.H)
    assert any(s.im(rho[i,j])!=0 for i in range(4) for j in range(4))
    assert s.simplify(reverse_density(reverse_density(rho))-rho)==s.zeros(4)
    assert s.simplify(reverse_density(rho)-U*rho*U.H)!=s.zeros(4)
    for species in ('quark','gluon'):
        paths=(('+',),('-',)) if species=='quark' else (('+','+'),('-','-'),('+','-'),('-','+'))
        for pair in paths:
            link=Link(species,pair)
            assert link.reverse().reverse()==link
            assert link.reverse().class_name==link.class_name
    tables=compare_sign_tables()
    assert (tables['quark']['octupole_odd'],tables['quark']['octupole_even'])==(3,11)
    assert (tables['gluon']['octupole_odd'],tables['gluon']['octupole_even'])==(10,4)
    return {'spin_state_square':'-I','density_involution':True,
            'complex_coherences':True,'link_pairs':6,'tables':tables,
            'analytic_input':'eq:PTprojection field/link transformation',
            'momentum_labels':'unchanged under combined PT'}


def conditional_evolution_check():
    # A common scalar operator commutes with the sign action. Generic channel
    # mixing does not: D M != M D when channels have opposite signs.
    t=s.symbols('t',real=True)
    D=s.diag(-1,1)
    common=s.exp(-t)*s.eye(2)
    mixing=s.Matrix([[1,1],[0,1]])
    assert D*common==common*D and D*mixing!=mixing*D
    return {'commuting_common_operator':True,'generic_mixing_rejected':True,
            'assumptions':['matched mu/zeta/scheme','common spin-independent operator',
                           'operator intertwines link prescriptions'],
            'qcd_kernel_computed':False}
