"""Read-only diagnostics for ordered lepton currents and the Born trace.

No function in this module changes the historical Born implementation.
The candidate contraction is evaluated separately from its source-labelled
reference.  All spinors are anchored before a tensor label is compared.
"""
from __future__ import annotations

import numpy as np
import sympy as s
from itertools import product
from born_direct import G, SIGMA, SIGN, massive, massless, fourvectors, direct
from gluon_born_response import leptonic, physical_amplitude_first_leptonic, evaluate


G5 = 1j*G[0]@G[1]@G[2]@G[3]


def exact_current_identity():
    """Generic eight-component gamma trace, independent of spinor phases."""
    from process_dirac import gamma
    g=gamma();g5=s.I*g[0]*g[1]*g[2]*g[3]
    l=s.symbols('l0:4');p=s.symbols('p0:4');h=s.symbols('h')
    signs=(1,-1,-1,-1)
    slash=lambda v:sum((signs[i]*v[i]*g[i] for i in range(4)),s.zeros(4))
    lo=[signs[i]*l[i] for i in range(4)]
    po=[signs[i]*p[i] for i in range(4)]
    dot=sum(lo[i]*p[i] for i in range(4))
    for mu,nu in product(range(4),repeat=2):
        actual=s.trace(slash(p)*(signs[mu]*g[mu])*(s.eye(4)+h*g5)*
                       slash(l)*(signs[nu]*g[nu])/2)
        expected=2*(lo[mu]*po[nu]+lo[nu]*po[mu]-
                    (signs[mu] if mu==nu else 0)*dot-
                    s.I*h*sum(-s.LeviCivita(mu,nu,a,b)*l[a]*p[b]
                              for a,b in product(range(4),repeat=2)))
        if s.expand(actual-expected)!=0:
            raise AssertionError(f'exact ordered gamma trace {mu},{nu}')
    if s.trace(g[0]*g[1]*g[2]*g[3]*g5)!=-4*s.I:
        raise AssertionError('exact gamma5 trace orientation')
    return {'components':16,'generic_four_vectors':'eight independent symbols',
            'J_antisymmetric_sign':'-i h epsilon_lower',
            'reverse_antisymmetric_sign':'+i h epsilon_lower',
            'gamma5_trace':'-4*I'}


def _slash(p):
    return sum((SIGN[a]*p[a]*G[a] for a in range(4)),np.zeros((4,4),complex))


def spinor_anchor(p, h):
    u=massless(p,h)
    unit=p[1:]/np.linalg.norm(p[1:])
    hel=np.block([[sum(unit[a]*SIGMA[a] for a in range(3)),np.zeros((2,2))],
                  [np.zeros((2,2)),sum(unit[a]*SIGMA[a] for a in range(3))]])
    projector=np.outer(u,u.conj())@G[0]
    expected=(np.eye(4)+h*G5)@_slash(p)/2
    return {'h':h,'helicity_eigen_residual':float(np.linalg.norm(hel@u-h*u)),
            'gamma5_residual':float(np.linalg.norm(G5@u-h*u)),
            'normalization':float(np.vdot(u,u).real),
            'projector_residual':float(np.linalg.norm(projector-expected))}


def annihilation_spinor_anchor():
    """Report DY electron and antiparticle spinor branch conventions separately.

    The antiparticle *spinor* has H and gamma5 eigenvalue -h in this
    implementation.  The physical positron-state helicity convention
    involves the creation-operator label and is not inferred here.
    """
    from process_normalization import helicity_spinor
    theta,phi,energy=1.1,.8,5.
    unit=np.array([np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi),np.cos(theta)])
    h2=sum(unit[a]*SIGMA[a] for a in range(3))
    hel=np.block([[h2,np.zeros((2,2))],[np.zeros((2,2)),h2]])
    rows=[]
    for anti in (False,True):
        for label in (-1,1):
            psi=helicity_spinor(label,theta,phi,anti,energy)
            eigen=-label if anti else label
            p=np.r_[energy,energy*unit]
            density=np.outer(psi,psi.conj())@G[0]
            expected=(np.eye(4)+eigen*G5)@_slash(p)/2
            row={'branch':'antiparticle_v' if anti else 'electron_u',
                 'input_label':label,'spinor_helicity_eigenvalue':eigen,
                 'gamma5_eigenvalue':eigen,
                 'helicity_residual':float(np.linalg.norm(hel@psi-eigen*psi)),
                 'gamma5_residual':float(np.linalg.norm(G5@psi-eigen*psi)),
                 'barred_projector_residual':float(np.linalg.norm(density-expected))}
            if max(row['helicity_residual'],row['gamma5_residual'],
                   row['barred_projector_residual'])>1e-12:
                raise AssertionError('annihilation spinor branch anchor')
            rows.append(row)
    return rows


def ordered_currents(l, lp, h):
    incoming=massless(l,h)
    lowered=[SIGN[mu]*G[mu] for mu in range(4)]
    J=np.zeros((4,4),complex)
    for hp in (-1,1):
        outgoing=massless(lp,hp)
        j=np.array([np.vdot(outgoing,G[0]@lowered[mu]@incoming)
                    for mu in range(4)])
        J+=np.outer(j,j.conj())
    reverse=J.T
    rho=(np.eye(4)+h*G5)@_slash(l)/2
    trace=np.array([[np.trace(_slash(lp)@lowered[mu]@rho@lowered[nu])
                     for nu in range(4)] for mu in range(4)])
    source=leptonic(l,lp,h)
    return {'J':J,'reverse':reverse,'source':source,
            'spinor_trace_residual':float(np.max(np.abs(J-trace))),
            'current_to_source_transpose':float(np.max(np.abs(J-source.T))),
            'current_to_literal_source':float(np.max(np.abs(J-source)))}


def sidis_ordering():
    """Trace(Phi gamma^i Delta gamma^j) has reverse hard indices.

    The two arbitrary complex spinors supply independent projectors; this
    is a matrix identity, not a fitted SIDIS structure-function fixture.
    """
    p=np.array([2.4,.7,.3,1.1]);q=np.array([2.8,-.4,.8,1.4])
    mass1=float(np.sqrt(p[0]**2-np.dot(p[1:],p[1:])))
    mass2=float(np.sqrt(q[0]**2-np.dot(q[1:],q[1:])))
    u=massive(p,mass1)[0]+.31j*massive(p,mass1)[1]
    v=massive(q,mass2)[0]-.23j*massive(q,mass2)[1]
    phi=np.outer(u,u.conj())@G[0]
    delta=np.outer(v,v.conj())@G[0]
    a=np.array([np.vdot(v,G[0]@G[i]@u) for i in (1,2)])
    w=np.array([[np.trace(phi@G[i]@delta@G[j]) for j in (1,2)] for i in (1,2)])
    hard=np.outer(a,a.conj())
    return {'trace_to_reverse_hard':float(np.max(np.abs(w-hard.T))),
            'trace_to_amplitude_first':float(np.max(np.abs(w-hard))),
            'complex_offdiagonal':float(abs(w[0,1].imag))}


def _matrix(report):
    return np.array(report['B_real'])+1j*np.array(report['B_imag'])


def born_case(inputs):
    h=inputs['helicity']
    physical=direct(**inputs)
    literal=_matrix(evaluate(**inputs))
    # Algebraic identity L_source(h)^T = L_source(-h); the latter call is
    # only a separate candidate trace evaluation, never a literal match.
    candidate=_matrix(evaluate(**dict(inputs,helicity=-h)))
    e=np.array([2,1+2j],complex);e=e/np.linalg.norm(e)
    circular=np.array([1,1j],complex)/np.sqrt(2)
    rate=lambda B,z:float(np.real(z@B@z.conj()))
    return {'inputs':inputs,'literal_max_abs':float(np.max(np.abs(physical-literal))),
            'candidate_max_abs':float(np.max(np.abs(physical-candidate))),
            'literal_b_G':float(literal[0,1].imag),
            'candidate_b_G':float(candidate[0,1].imag),
            'physical_b_G':float(physical[0,1].imag),
            'unchanged_UCS_max_abs':float(max(abs(physical[0,0].real-literal[0,0].real),
                                         abs(physical[1,1].real-literal[1,1].real),
                                         abs(physical[0,1].real-literal[0,1].real))),
            'circular_rate_literal':rate(literal,circular),
            'circular_rate_physical':rate(physical,circular),
            'elliptic_rate_literal':rate(literal,e),
            'elliptic_rate_physical':rate(physical,e),
            'physical_eigenvalues':np.linalg.eigvalsh(physical).tolist(),
            'literal_eigenvalues':np.linalg.eigvalsh(literal).tolist()}


def run():
    exact_identity=exact_current_identity()
    event_l=np.array([5.,0.,0.,5.]);event_lp=np.array([5.,3.,0.,4.])
    exact=ordered_currents(event_l,event_lp,1)
    if not (abs(exact['J'][1,2]+10j)<1e-12 and
            abs(exact['source'][1,2]-10j)<1e-12):
        raise AssertionError('exact ordered-current event')
    if abs(np.trace(G[0]@G[1]@G[2]@G[3]@G5)+4j)>1e-12:
        raise AssertionError('gamma5 trace orientation')
    if not (leptonic(event_l,event_lp,1)[1,2]==10j and
            abs(physical_amplitude_first_leptonic(event_l,event_lp,1)[1,2]+10j)<1e-12):
        raise AssertionError('epsilon lower-index orientation')
    cases=[];anchors=[]
    for theta,phi in ((.8,.4),(1.2,1.1),(.55,2.2)):
        ps=fourvectors(theta=theta,phi=phi)
        for az in (.0,.73):
            rotation=np.array([[np.cos(az),-np.sin(az),0],
                               [np.sin(az),np.cos(az),0],[0,0,1]])
            l=ps['l'].copy();lp=ps['lp'].copy()
            l[1:]=rotation@l[1:];lp[1:]=rotation@lp[1:]
            for h in (-1,1):
                anchors.append(spinor_anchor(l,h))
                currents=ordered_currents(l,lp,h)
                if max(currents['spinor_trace_residual'],
                       currents['current_to_source_transpose'])>2e-12:
                    raise AssertionError('current/trace/source transpose')
        for h in (-1,1):
            case=born_case(dict(sqrt_s=5.,Q2=4.,mass=1.5,theta=theta,
                                phi=phi,lepton_energy=10.,helicity=h))
            if not (case['literal_max_abs']>1 and case['candidate_max_abs']<1e-8):
                raise AssertionError('literal mismatch and candidate agreement required')
            cases.append(case)
    sidis=sidis_ordering()
    if sidis['trace_to_reverse_hard']>1e-12 or sidis['complex_offdiagonal']<1e-5:
        raise AssertionError('SIDIS hard-index order')
    if any(max(x['helicity_eigen_residual'],x['gamma5_residual'],
               x['projector_residual'])>1e-12 for x in anchors):
        raise AssertionError('physical spinor anchor')
    return {'exact_gamma_identity':exact_identity,
            'exact_event':{'l':event_l.tolist(),'lp':event_lp.tolist(),
                           'J_xy_imag':float(exact['J'][1,2].imag),
                           'source_L_xy_imag':float(exact['source'][1,2].imag)},
            'gamma5_trace_0123_imag':-4.0,
            'epsilon_upper_0123':1,'epsilon_lower_0123':-1,
            'anchors':anchors,'annihilation_spinors':annihilation_spinor_anchor(),
            'born_cases':cases,'sidis':sidis,
            'status':{'diagnostic_execution':'PASS',
                      'literal_physical_source_agreement':'MISMATCH',
                      'candidate_index_interchange':'AGREES',
                      'publication_eligibility':'BLOCKED_AUTHOR_REVIEW'}}
