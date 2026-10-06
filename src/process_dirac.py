"""Four-dimensional leading-power Dirac algebra for the two photon processes.

The Dirac and chiral matrices are built independently from Pauli blocks.
The latter is a representation check, not an independent QCD derivation.
Metric (+---), epsilon^0123=+1, gamma5=i gamma0 gamma1 gamma2 gamma3.
"""
from __future__ import annotations

import sympy as s
from functools import lru_cache
from sympy.physics.matrices import msigma

I=s.I
Z=s.zeros(2)
E=s.eye(2)


def gamma(representation='dirac'):
    if representation == 'dirac':
        g0=s.diag(1,1,-1,-1)
        spatial=[s.BlockMatrix([[Z,p],[-p,Z]]).as_explicit() for p in (msigma(1),msigma(2),msigma(3))]
    elif representation == 'chiral':
        g0=s.BlockMatrix([[Z,E],[E,Z]]).as_explicit()
        spatial=[s.BlockMatrix([[Z,p],[-p,Z]]).as_explicit() for p in (msigma(1),msigma(2),msigma(3))]
    else: raise ValueError('unknown gamma representation')
    return (g0,*spatial)


def basis(representation='dirac'):
    g=gamma(representation)
    g5=I*g[0]*g[1]*g[2]*g[3]
    gp=(g[0]+g[3])/s.sqrt(2)
    gm=(g[0]-g[3])/s.sqrt(2)
    def sigma(a,b):return I*(a*b-b*a)/2
    return g,g5,gp,gm,sigma


def reconstruction(side, channel, representation='dirac'):
    """Unit projection: side=A (plus mover), B (minus mover), or FF."""
    g,g5,gp,gm,sigma=basis(representation)
    large=gm if side=='A' else gp
    if channel=='F' or channel=='D':return large/2
    if channel=='G':return g5*large/2
    if channel.startswith('T'):
        axis=int(channel[-1])
        return I*sigma(g[axis],large)*g5/2
    raise ValueError(channel)


def projection(side, channel, representation='dirac'):
    g,g5,gp,gm,sigma=basis(representation)
    opposite=gp if side=='A' else gm
    if channel in ('F','D'):return opposite
    if channel=='G':return opposite*g5
    if channel.startswith('T'):
        return I*sigma(g[int(channel[-1])],opposite)*g5
    raise ValueError(channel)


def hard_trace(left, right, i, j, representation='dirac'):
    g=gamma(representation)
    return s.simplify(s.trace(reconstruction('A',left,representation)*g[i]*
                              reconstruction('B',right,representation)*g[j]))


@lru_cache(maxsize=None)
def pair_table(process, representation='dirac'):
    left=('F','G','T1','T2')
    right=('D','T1','T2') if process=='SIDIS' else left
    return {(a,b):s.Matrix(2,2,lambda i,j:hard_trace(a,b,i+1,j+1,representation))
            for a in left for b in right}


def check_algebra():
    out={}
    metric=(1,-1,-1,-1)
    for representation in ('dirac','chiral'):
        g,g5,gp,gm,sigma=basis(representation)
        for a in range(4):
            for b in range(4):
                assert g[a]*g[b]+g[b]*g[a]==2*metric[a]*int(a==b)*s.eye(4)
                assert sigma(g[a],g[b])==I*(g[a]*g[b]-g[b]*g[a])/2
            assert g[a].H==g[0]*g[a]*g[0]
            assert g5*g[a]+g[a]*g5==s.zeros(4)
        assert g5**2==s.eye(4) and g5.H==g5
        assert gp*gp==s.zeros(4) and gm*gm==s.zeros(4)
        for side in ('A','B'):
            for channel in ('F','G','T1','T2'):
                for other in ('F','G','T1','T2'):
                    value=s.trace(reconstruction(side,channel,representation)*
                                  projection(side,other,representation))/2
                    assert s.simplify(value-int(channel==other))==0,(representation,side,channel,other,value)
        out[representation]={'clifford_pairs':16,'projection_pairs':32,'gamma5':'four-dimensional'}
    conversion=s.BlockMatrix([[E,-E],[E,E]]).as_explicit()/s.sqrt(2)
    assert conversion*conversion.H==s.eye(4)
    assert all(conversion*gamma('dirac')[i]*conversion.H==gamma('chiral')[i]
               for i in range(4))
    out['basis_conversion']={'explicit_unitary':True,
                             'construction':'independent Pauli-block Dirac and chiral matrices'}
    for process in ('SIDIS','DY'):
        a=pair_table(process,'dirac');b=pair_table(process,'chiral')
        assert a==b
        delta=s.eye(2)
        epsilon=s.Matrix([[0,1],[-1,0]])
        for pair,matrix in a.items():
            left,right=pair
            if left in ('F','G') and right in ('D','F','G'):
                target=I*epsilon if (left=='G')^(right=='G') else delta
            elif left.startswith('T') and right.startswith('T'):
                p=int(left[-1])-1;q=int(right[-1])-1
                target=s.Matrix(2,2,lambda i,j:-int(i==p)*int(j==q)-
                                int(j==p)*int(i==q)+int(i==j)*int(p==q))
            else:target=s.zeros(2)
            assert matrix==target,(process,pair,matrix,target)
        out[process]={'pairs':len(a),'transverse_components':4*len(a),
                      'independent_representation_equal':True,
                      'curated_tensor_structure_equal':True}
    return out
