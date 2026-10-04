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
            'basis_order':'parton Cartesian/helicity outer index, target m=3/2,1/2,-1/2,-3/2 inner index',
            'coefficient_rank_preserved':True,
            'matrix_rank_or_psd_preserved':False,
            'partial_transpose_psd_preserving':False,
            'Bell_partial_transpose_negative_eigenvalue':'-1/2',
            'minimal_entries':examples,
            'source_labels':['eq:qD','eq:gD','eq:jointM']}

def spectral_index_check():
    """Derive source target contraction and hard rate from amplitude labels.

    The source joint indices are (i,Lambda;j,Lambda').  A row is a
    removal amplitude A_X,i,Lambda, so M=A.H*A.  The target density is
    rho_Lambda,Lambda'=c_Lambda*c_Lambda'^* in the same basis.
    """
    rho_seed=s.Matrix([[1,s.I,0,1],[s.I,2,1,0],[0,1,1,s.I],[1,0,s.I,1]])
    rho=rho_seed*rho_seed.H/s.trace(rho_seed*rho_seed.H)
    if rho==rho.T:raise AssertionError('mixed target must probe imaginary coherence')
    witness={}
    for species in ('quark','gluon'):
        for col in range(32):
            M=source_covariant(species,col)
            D=s.Matrix(2,2,lambda i,j:s.trace(rho*M[4*i:4*i+4,4*j:4*j+4]))
            indexed=s.Matrix(2,2,lambda i,j:sum(
                rho[b,a]*M[4*i+a,4*j+b] for a in range(4) for b in range(4)))
            if s.simplify(D-indexed)!=s.zeros(2):
                raise AssertionError(f'{species} spectral mixed-index contraction {col}')
            if col==0:witness[species]=str(s.simplify(D[0,0]))
    # Explicit complex amplitude labels and a circular hard analyzer.
    A=s.Matrix(3,8,lambda x,y:s.Rational((x+2)*(y+1)%7-3,5)+
               s.I*s.Rational((2*x+y)%5-2,7))
    M=A.H*A
    c=s.Matrix([1,s.I,2,-s.I])/s.sqrt(7)
    e=s.Matrix([1,s.I])/s.sqrt(2)
    H=s.Matrix([1+2*s.I,2-s.I])
    D=s.Matrix(2,2,lambda i,j:sum(s.conjugate(c[a])*M[4*i+a,4*j+b]*c[b]
                                 for a in range(4) for b in range(4)))
    amplitude=s.Matrix([sum(H[i]*e[i]*A[x,4*i+a]*c[a]
                            for i in range(2) for a in range(4)) for x in range(3)])
    rate=s.expand(sum(s.conjugate(amplitude[x])*amplitude[x] for x in range(3)))
    source_rate=sum(s.conjugate(H[i]*e[i]*c[a])*M[4*i+a,4*j+b]*
                    H[j]*e[j]*c[b]
                    for i in range(2) for j in range(2)
                    for a in range(4) for b in range(4))
    if s.simplify(rate-source_rate)!=0:raise AssertionError('source spectral rate')
    return {'columns_checked':{'quark':32,'gluon':32},
            'mixed_target_imaginary_coherence':True,
            'mixed_contraction':'D_ij = sum rho_(Lambda prime,Lambda) M_(i Lambda,j Lambda prime)',
            'spectral_matrix':'M_(i Lambda,j Lambda prime) = sum_X A*_(X,i,Lambda) A_(X,j,Lambda prime)',
            'hard_rate':'sum_X |sum_(i,Lambda) H_i e_i A_(X,i,Lambda) c_Lambda|^2',
            'source_operator_witness':witness,
            'amplitude_rate':str(s.simplify(rate))}

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
    # A nontrivial source-consistent finite-kT Gram need not remain PSD
    # under the coordinate-map partial transpose.  This is a negative
    # control on treating the auxiliary matrix as a physical density.
    import numpy as np
    pt=partial_parton_transpose(M)
    pt_min=float(np.linalg.eigvalsh(np.array(pt.evalf(18),dtype=complex)).min())
    if pt_min>=-1e-3:raise AssertionError('finite-kT partial-transpose control lost')
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
            'finite_k_partial_transpose_min_eigenvalue':pt_min,
            'kinematics':['2','0','3'],
            'factorization':'(A.H*A + U*A.H*A*U.H)/2; both summands PSD',
            'partial_transpose_applied_to_Gram':False}
