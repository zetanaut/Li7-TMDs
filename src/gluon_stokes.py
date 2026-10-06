"""Index-safe gluon Gram, Stokes reconstruction, and basis covariance."""
from __future__ import annotations
import sympy as s
import numpy as np


def density(F,G,C,S):
    return s.Matrix([[F+C,S+s.I*G],[S-s.I*G,F-C]])/2


def stokes(B):
    return ((B[0,0]+B[1,1])/2,(B[0,1]-B[1,0])/(2*s.I),
            (B[0,0]-B[1,1])/2,(B[0,1]+B[1,0])/2)


def reconstruct(bU,bG,bC,bS):
    return s.Matrix([[bU+bC,bS+s.I*bG],[bS-s.I*bG,bU-bC]])


def exact_checks():
    u,g,c,z,F,G,C,S=s.symbols('u g c z F G C S',real=True)
    B=reconstruct(u,g,c,z);D=density(F,G,C,S)
    assert stokes(B)==(u,g,c,z)
    assert s.expand(s.trace(B*D)-(u*F+g*G+c*C+z*S))==0
    # The characteristic polynomial establishes the eigenvalues without
    # asking the symbolic eigensolver to choose square-root branches.
    lam=s.symbols('lambda')
    assert s.expand((lam*s.eye(2)-B).det()-((lam-u)**2-g*g-c*c-z*z))==0
    return {'source_labels':['eq:gD','eq:hardcontraction','eq:Bhard','eq:bproject'],
            'symbolic_entries':4,'characteristic_polynomial':True}


def numerical_checks():
    # These are amplitudes, not arbitrary Hermitian entries; B is positive.
    amplitudes=np.array([[1+2j,.3-.4j],[.7+0.2j,-.6+1.1j]],complex)
    B=amplitudes@amplitudes.conj().T
    examples={
        'circular_plus':np.array([1,1j])/np.sqrt(2),
        'circular_minus':np.array([1,-1j])/np.sqrt(2),
        'elliptic':np.array([2,1+2j])/3,
        'real_x':np.array([1,0],complex),
    }
    rates={}
    for name,e in examples.items():
        e=e/np.linalg.norm(e)
        D=np.outer(e.conj(),e)
        lhs=np.trace(B@D)
        rhs=np.sum(np.abs(e@amplitudes)**2)
        assert abs(lhs-rhs)<2e-14,(name,lhs,rhs)
        if name.startswith('circular'):
            wrong=np.trace(B@D.T)
            assert abs(lhs-wrong)>1e-3
        rates[name]=float(lhs.real)
    mixed=(.37*np.outer(examples['circular_plus'].conj(),examples['circular_plus'])
           +.63*np.outer(examples['circular_minus'].conj(),examples['circular_minus']))
    assert np.trace(B@mixed).real>0
    U=np.array([[1,1j],[1,-1j]],complex)/np.sqrt(2)
    assert np.allclose(U.conj().T@U,np.eye(2))
    Anew=U.T@amplitudes
    Bnew=U.T@B@U.conj()
    for e in examples.values():
        e=e/np.linalg.norm(e)
        ep=U.conj().T@e
        D=np.outer(e.conj(),e)
        Dp=U.T@D@U.conj()
        assert np.allclose(Anew.T@ep,amplitudes.T@e)
        assert abs(np.trace(Bnew@Dp)-np.trace(B@D))<2e-14
    eig=np.linalg.eigvalsh(B)
    bu=(B[0,0].real+B[1,1].real)/2
    bg=B[0,1].imag;bc=(B[0,0].real-B[1,1].real)/2;bs=B[0,1].real
    assert np.allclose(eig,[bu-np.linalg.norm([bg,bc,bs]),bu+np.linalg.norm([bg,bc,bs])])
    assert eig.min()>0 and bu>0
    return {'rates':rates,'complex_polarizations':3,'mixed_state':True,
            'basis_covariance':True,'positive_eigenvalues':eig.tolist(),
            'stokes':[float(x) for x in (bu,bg,bc,bs)]}
