#!/usr/bin/env python3
"""Born hard Stokes response for l g -> l' Q Qbar.

Requires Python >=3.10 and NumPy. No network access or TMD model is used.
Returned B is coupling-stripped: multiply by e^4 e_Q^2 g_s^2 T_F / Q^4,
then apply the appropriate phase-space and flux factors for the observable. The initial gluon spin
average is supplied by D_in=I/2, not by another factor in B.

Example:
    python src/gluon_born_response.py --output results/gluon_born_report.json
Conventions: g=(+---), epsilon^{0123}=+1, E_xy=+1, i at the
conjugate field in the gluon Gram correlator. Rate = trace(B D_g).
"""
from __future__ import annotations
import argparse
import itertools
import json
from pathlib import Path
import numpy as np
from numpy.typing import NDArray

CMatrix = NDArray[np.complex128]
RVector = NDArray[np.float64]
METRIC = np.diag([1., -1., -1., -1.])
I4 = np.eye(4, dtype=complex)

def gamma_matrices() -> list[CMatrix]:
    """Dirac representation of gamma^mu, with {gamma^mu,gamma^nu}=2g^mu nu."""
    pauli = [np.array([[0,1],[1,0]],complex),
             np.array([[0,-1j],[1j,0]],complex),
             np.array([[1,0],[0,-1]],complex)]
    z = np.zeros((2,2),complex)
    return [np.diag([1,1,-1,-1]).astype(complex)] + [
        np.block([[z,p],[-p,z]]) for p in pauli]

GAMMA = gamma_matrices()

def dot(a: RVector,b: RVector) -> float:
    return float(a @ METRIC @ b)

def slash(p: RVector) -> CMatrix:
    return sum((GAMMA[mu]*(METRIC @ p)[mu] for mu in range(4)),np.zeros((4,4),complex))

def bar(a: CMatrix) -> CMatrix:
    return GAMMA[0] @ a.conj().T @ GAMMA[0]

def epsilon_lower(mu: int,nu: int,rho: int,sigma: int) -> int:
    """epsilon_0123=-1, consistent with epsilon^0123=+1."""
    idx=(mu,nu,rho,sigma)
    if len(set(idx)) < 4:
        return 0
    inversions=sum(idx[a] > idx[b] for a in range(4) for b in range(a+1,4))
    return -(-1)**inversions

def leptonic(l: RVector,lp: RVector,helicity: float) -> CMatrix:
    """Historical source L_mu,nu, with conjugate-current-first ordering.

    ``helicity`` retains the published source label.  For the physical
    incoming electron helicity h in an amplitude-first Born trace, the
    required current tensor is this matrix transposed at lambda=h.
    """
    lo,lpo=METRIC @ l,METRIC @ lp
    result=2*(np.outer(lo,lpo)+np.outer(lpo,lo)-METRIC*dot(l,lp)).astype(complex)
    for mu,nu in itertools.product(range(4),repeat=2):
        result[mu,nu] += 2j*helicity*sum(
            epsilon_lower(mu,nu,r,s)*l[r]*lp[s]
            for r,s in itertools.product(range(4),repeat=2))
    return result

def physical_amplitude_first_leptonic(l: RVector,lp: RVector,h: float) -> CMatrix:
    """Candidate J_mu,nu=sum j_mu j_nu*, derived by index interchange.

    This is diagnostic only; ``evaluate`` continues to implement the
    approved source formula and historical reference values.
    """
    return leptonic(l,lp,h).T

def kinematics(sqrt_s: float,Q2: float,mass: float,theta: float,
               phi: float,lepton_energy: float) -> dict[str,RVector]:
    if not (np.isfinite([sqrt_s,Q2,mass,theta,phi,lepton_energy]).all()
            and sqrt_s > 2*mass > 0 and Q2 > 0 and lepton_energy > 0):
        raise ValueError('Require finite inputs, Q2>0, mass>0, sqrt_s>2*mass and lepton_energy>0.')
    s=sqrt_s**2
    eg=(s+Q2)/(2*sqrt_s)
    eq=(s-Q2)/(2*sqrt_s)
    k=np.array([eg,0.,0.,eg])
    q=np.array([eq,0.,0.,-eg])
    p=np.sqrt(s/4-mass**2)
    p1=np.array([sqrt_s/2,p*np.sin(theta)*np.cos(phi),
                 p*np.sin(theta)*np.sin(phi),p*np.cos(theta)])
    p2=np.array([sqrt_s/2,-p1[1],-p1[2],-p1[3]])
    lz=(-Q2/2-eq*lepton_energy)/eg
    lx2=lepton_energy**2-lz**2
    if lx2 <= 0 or lepton_energy-eq <= 0:
        raise ValueError('This lepton energy cannot realize the specified virtual-photon kinematics.')
    l=np.array([lepton_energy,np.sqrt(lx2),0.,lz])
    return dict(l=l,lp=l-q,q=q,k=k,p1=p1,p2=p2)

def evaluate(*,sqrt_s: float=5.,Q2: float=4.,mass: float=1.5,
             theta: float=.8,phi: float=.4,lepton_energy: float=10.,
             helicity: float=1.) -> dict[str,object]:
    if not np.isfinite(helicity) or abs(helicity)>1:
        raise ValueError('The lepton polarization/helicity must lie in [-1,1].')
    ps=kinematics(sqrt_s,Q2,mass,theta,phi,lepton_energy)
    l,lp,q,k,p1,p2=(ps[key] for key in ('l','lp','q','k','p1','p2'))
    d1=dot(p1-q,p1-q)-mass**2
    d2=dot(p1-k,p1-k)-mass**2
    if min(abs(d1),abs(d2)) < 1e-12:
        raise ValueError('Propagator denominator is too close to zero for a stable test.')
    v=np.zeros((4,4,4,4),complex)
    for mu,alpha in itertools.product(range(4),repeat=2):
        v[mu,alpha]=(GAMMA[mu] @ (slash(p1-q)+mass*I4) @ GAMMA[alpha]/d1
                     +GAMMA[alpha] @ (slash(p1-k)+mass*I4) @ GAMMA[mu]/d2)
    u=slash(p1)+mass*I4
    w=slash(p2)-mass*I4
    lepton=leptonic(l,lp,helicity)
    b=np.zeros((2,2),complex)
    for i,j in itertools.product(range(2),repeat=2):
        b[i,j]=sum(lepton[mu,nu]*np.trace(u @ v[mu,i+1] @ w @ bar(v[nu,j+1]))
                   for mu,nu in itertools.product(range(4),repeat=2))
    # The following are amplitude Ward tests between on-shell spin projectors,
    # not the stronger (and false) demand that the bare vertex itself vanish.
    qward=max(np.linalg.norm(u @ sum((q[mu]*METRIC[mu,mu]*v[mu,a]
                  for mu in range(4)),np.zeros((4,4),complex)) @ w) for a in range(4))
    kward=max(np.linalg.norm(u @ sum((k[a]*METRIC[a,a]*v[mu,a]
                  for a in range(4)),np.zeros((4,4),complex)) @ w) for mu in range(4))
    amp_scale=max(1.,max(np.linalg.norm(u @ v[mu,a] @ w)
                           for mu,a in itertools.product(range(4),repeat=2)))
    bscale=max(1.,float(np.linalg.norm(b)))
    herm_res=float(np.linalg.norm(b-b.conj().T)/bscale)
    bherm=(b+b.conj().T)/2
    eigenvalues=np.linalg.eigvalsh(bherm)
    lorentz_errors={
        'l_squared':abs(dot(l,l)), 'lp_squared':abs(dot(lp,lp)),
        'k_squared':abs(dot(k,k)), 'q_squared_plus_Q2':abs(dot(q,q)+Q2),
        'p1_squared_minus_m2':abs(dot(p1,p1)-mass**2),
        'p2_squared_minus_m2':abs(dot(p2,p2)-mass**2),
        'momentum_conservation':float(np.linalg.norm(l+k-lp-p1-p2))}
    relative_wards={'photon':float(qward/(amp_scale*max(1.,np.linalg.norm(q)))),
                    'gluon':float(kward/(amp_scale*max(1.,np.linalg.norm(k))))}
    checks={
        'on_shell_and_conservation':max(lorentz_errors.values()) < 1e-9,
        'photon_Ward':relative_wards['photon'] < 1e-11,
        'gluon_Ward':relative_wards['gluon'] < 1e-11,
        'Hermitian_hard_matrix':herm_res < 1e-11,
        'positive_semidefinite_hard_matrix':bool(eigenvalues.min() >= -1e-11*bscale)}
    if not all(checks.values()):
        raise ArithmeticError(f'Physical hard-response check failed: {checks}')
    return {
        'description':'Coupling-stripped Born hard matrix; multiply by e^4 e_Q^2 g_s^2 T_F / Q^4.',
        'inputs':dict(sqrt_s=sqrt_s,Q2=Q2,mass=mass,theta=theta,phi=phi,
                      lepton_energy=lepton_energy,helicity=helicity),
        'momenta_GeV':{name:p.tolist() for name,p in ps.items()},
        'propagators_GeV2':[d1,d2],
        'B_real':b.real.tolist(),'B_imag':b.imag.tolist(),
        'Stokes':{'b_U':float((b[0,0].real+b[1,1].real)/2),
                  'b_G':float(b[0,1].imag),
                  'b_C':float((b[0,0].real-b[1,1].real)/2),
                  'b_S':float(b[0,1].real)},
        'B_eigenvalues':eigenvalues.tolist(),
        'relative_Ward_residuals':relative_wards,
        'relative_Hermiticity_residual':herm_res,
        'kinematic_residuals':lorentz_errors,
        'checks':checks}

def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--sqrt-s',type=float,default=5.)
    parser.add_argument('--Q2',type=float,default=4.)
    parser.add_argument('--mass',type=float,default=1.5)
    parser.add_argument('--theta',type=float,default=.8,help='Heavy-quark polar angle in radians.')
    parser.add_argument('--phi',type=float,default=.4,help='Heavy-quark azimuth in radians.')
    parser.add_argument('--lepton-energy',type=float,default=10.)
    parser.add_argument('--helicity',type=float,default=1.)
    parser.add_argument('--output',type=Path,default=Path('results/gluon_born_report.json'))
    args=parser.parse_args()
    try:
        report=evaluate(sqrt_s=args.sqrt_s,Q2=args.Q2,mass=args.mass,theta=args.theta,
                        phi=args.phi,lepton_energy=args.lepton_energy,helicity=args.helicity)
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    except (ValueError,ArithmeticError,OSError) as exc:
        parser.exit(1,f'Error: {exc}\n')
    print(json.dumps(report['Stokes'],indent=2))
    print('All five hard-response checks PASS.')
    print(f'Report: {args.output}')

if __name__=='__main__':
    main()
