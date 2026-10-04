"""Literal source-index joint correlator and conditional spectral examples.

The pre-existing foundation map used the opposite sigma_2 convention.
Here the conversion is checked explicitly; partial transpose is never used
as a positivity-preserving operation on a test matrix.
"""
from __future__ import annotations
import sympy as s
from sympy.physics.matrices import msigma
from sympy.physics.wigner import wigner_d_small
from correlator_foundations import (catalogue,coordinates,helicity_value,
                                    target_component_operators,joint_covariant)

def source_covariant(species,column,kx=s.Integer(2),ky=s.Integer(1),mass=s.Integer(3)):
    labels=catalogue(species);label=labels[column]
    ops=target_component_operators()
    pauli=([s.eye(2),msigma(3),msigma(1),-msigma(2)] if species=='quark' else
           [s.eye(2),-msigma(2),msigma(3),msigma(1)])
    out=s.zeros(8)
    for K,m,part in coordinates():
        if (K,m)!=(label.K,label.m):continue
        vals=helicity_value(label,part,kx,ky,mass)
        out+=sum((s.kronecker_product(vals[a]*pauli[a],ops[K,m,part])/2
                  for a in range(4)),s.zeros(8))
    return s.simplify(out)

def partial_parton_transpose(M):
    return s.BlockMatrix([[M[4*b:4*b+4,4*a:4*a+4]
                           for b in range(2)] for a in range(2)]).as_explicit()

def source_mapping_check():
    counts={};examples=[]
    for species in ('quark','gluon'):
        for i in range(32):
            source=source_covariant(species,i)
            old=joint_covariant(species,i,s.Integer(2),s.Integer(1),s.Integer(3))
            if s.simplify(source-partial_parton_transpose(old)/2)!=s.zeros(8):
                raise AssertionError(f'{species} parton transpose convention col={i}')
            if source!=source.H:
                raise AssertionError('literal source matrix Hermiticity')
            if catalogue(species)[i].id() in (
                'quark.h.11[0]','gluon.g.10[0]'):
                source_unscaled=2*source
                changed=[(a,b) for a in range(8) for b in range(8)
                         if s.simplify(source_unscaled[a,b]-old[a,b])!=0]
                if not changed:raise AssertionError('imaginary source convention not probed')
                a,b=changed[0]
                examples.append({'label':catalogue(species)[i].id(),'indices':[a,b],
                                 'source_twice':str(source_unscaled[a,b]),
                                 'old_auxiliary':str(old[a,b])})
        counts[species]=32
    bell=s.zeros(8,1);bell[0]=bell[5]=1/s.sqrt(2)
    gram=bell*bell.H
    pt=partial_parton_transpose(gram)
    if set(pt.eigenvals())!={s.Rational(1,2),-s.Rational(1,2),s.Integer(0)}:
        raise AssertionError('partial transpose Bell witness')
    return {'source_indexed_columns':counts,
            'conversion':'source = partial-parton-transpose(foundation auxiliary map)/2',
            'partial_transpose_psd_preserving':False,
            'Bell_partial_transpose_negative_eigenvalue':'-1/2',
            'minimal_entries':examples,
            'source_labels':['eq:qD','eq:gD','eq:jointM']}

def parity_unitary(species):
    uy=wigner_d_small(s.Rational(3,2),-s.pi)
    up=-s.I*msigma(2) if species=='quark' else s.I*msigma(3)
    return s.kronecker_product(up,uy)

def spectral_recovery(species):
    W=[source_covariant(species,i,s.Integer(2),s.Integer(0),s.Integer(3)) for i in range(32)]
    U=parity_unitary(species)
    for col in W:
        if s.simplify(U*col*U.H-col)!=s.zeros(8):
            raise AssertionError('source parity space disagrees with foundation map')
    A=s.Matrix(3,8,lambda r,c:s.Rational((r+2)*(c+1)%7-3,5)+
               s.I*s.Rational((2*r+c)%5-2,7))
    raw=A.H*A
    M=s.simplify((raw+U*raw*U.H)/2)
    if M!=M.H or s.simplify(U*M*U.H-M)!=s.zeros(8):
        raise AssertionError('parity averaged spectral Gram')
    columns=s.Matrix.hstack(*(s.Matrix(64,1,lambda j,_:W[i][j//8,j%8]) for i in range(32)))
    target=s.Matrix(64,1,lambda j,_:M[j//8,j%8])
    pivots=columns.T.rref()[1]
    if len(pivots)!=32:raise AssertionError('source map rank')
    recovered=columns.extract(pivots,range(32)).inv()*target.extract(pivots,[0])
    if s.simplify(columns*recovered-target)!=s.zeros(64,1):
        raise AssertionError('positive source Gram not recovered by covariant basis')
    if any(s.im(s.simplify(c))!=0 for c in recovered):
        raise AssertionError('recovered coefficients nonreal')
    return {'species':species,'matrix_rank':M.rank(),'source_map_rank':32,
            'recovered_coefficients':32,'parity_invariant':True,
            'raw_Gram_rank':raw.rank(),'spectral_rows':3,
            'kinematics':['2','0','3'],
            'factorization':'(A.H*A + U*A.H*A*U.H)/2; both summands PSD',
            'partial_transpose_applied_to_Gram':False}
