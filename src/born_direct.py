"""Independent two-diagram spinor amplitudes for the reduced Born matrix.

The kinematics and gamma matrices here are constructed without the production
trace implementation.  ``compare`` imports that implementation only after the
amplitude sum has been formed.
"""
from __future__ import annotations
from itertools import product
import numpy as np

SIGMA=(np.array([[0,1],[1,0]],complex),
       np.array([[0,-1j],[1j,0]],complex),
       np.array([[1,0],[0,-1]],complex))
G=(np.diag([1,1,-1,-1]).astype(complex),)+tuple(
    np.block([[np.zeros((2,2)),a],[-a,np.zeros((2,2))]]) for a in SIGMA)
SIGN=np.array([1.,-1.,-1.,-1.])


def product4(a,b):
    return float(a[0]*b[0]-np.dot(a[1:],b[1:]))


def fourvectors(*,sqrt_s=5.,Q2=4.,mass=1.5,theta=.8,phi=.4,lepton_energy=10.):
    vals=np.array([sqrt_s,Q2,mass,theta,phi,lepton_energy],float)
    if not np.all(np.isfinite(vals)) or Q2<=0 or mass<=0 or sqrt_s<=2*mass:
        raise ValueError('finite Q2>0, mass>0, and sqrt_s>2*mass are required')
    if not 0<theta<np.pi:
        raise ValueError('an interior heavy-quark polar angle is required')
    if lepton_energy<=sqrt_s/2:
        raise ValueError('interior lepton energy must exceed sqrt_s/2')
    omega=(sqrt_s**2+Q2)/(2*sqrt_s)
    q0=(sqrt_s**2-Q2)/(2*sqrt_s)
    k=np.array([omega,0.,0.,omega])
    q=np.array([q0,0.,0.,-omega])
    p=np.sqrt(sqrt_s**2/4-mass**2)
    direction=np.array([np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi),np.cos(theta)])
    p1=np.r_[sqrt_s/2,p*direction]
    p2=np.r_[sqrt_s/2,-p*direction]
    lz=(-Q2/2-lepton_energy*q0)/omega
    lt2=lepton_energy**2-lz**2
    if lt2<=0 or lepton_energy-q0<=0:
        raise ValueError('no interior positive-energy massless lepton solution')
    l=np.array([lepton_energy,np.sqrt(lt2),0.,lz])
    return dict(l=l,lp=l-q,k=k,q=q,p1=p1,p2=p2)


def slash(v):
    return sum((SIGN[a]*v[a]*G[a] for a in range(4)),np.zeros((4,4),complex))


def pauli_dot(v):
    return sum((v[a]*SIGMA[a] for a in range(3)),np.zeros((2,2),complex))


def massive(p,m,anti=False):
    a=np.sqrt(p[0]+m)
    out=[]
    for chi in np.eye(2,dtype=complex).T:
        small=pauli_dot(p[1:])@chi/a
        out.append(np.r_[small,a*chi] if anti else np.r_[a*chi,small])
    return out


def massless(p,helicity):
    if helicity not in (-1,1):raise ValueError('definite helicity is +/-1')
    unit=p[1:]/np.linalg.norm(p[1:])
    _,vectors=np.linalg.eigh(pauli_dot(unit))
    chi=vectors[:,1 if helicity==1 else 0]
    return np.sqrt(p[0])*np.r_[chi,helicity*chi]


def _diagrams(ps,m,photon,gluon):
    q,k,p1,p2=(ps[x] for x in ('q','k','p1','p2'))
    us=massive(p1,m);vs=massive(p2,m,True)
    first=slash(p1-q)+m*np.eye(4)
    second=slash(p1-k)+m*np.eye(4)
    d1=product4(p1-q,p1-q)-m*m
    d2=product4(p1-k,p1-k)-m*m
    out=np.zeros((2,2,2),complex)
    for d,a,b in product(range(2),range(2),range(2)):
        for mu,alpha in product(range(4),range(4)):
            chain=(G[mu]@first@G[alpha]/d1 if d==0
                   else G[alpha]@second@G[mu]/d2)
            out[d,a,b]+=SIGN[mu]*SIGN[alpha]*photon[mu]*gluon[alpha]*np.vdot(us[a],G[0]@chain@vs[b])
    return out


def direct(*,helicity=1.,**inputs):
    if helicity not in (-1.,0.,1.):
        raise ValueError('direct route supports -1, 0, or +1 beam polarization')
    ps=fourvectors(**inputs)
    if helicity==0:
        plus=direct(helicity=1.,**inputs)
        minus=direct(helicity=-1.,**inputs)
        return (plus+minus)/2
    m=inputs.get('mass',1.5)
    return direct_vectors(ps,m,int(helicity))


def direct_vectors(ps,m,helicity):
    l,lp=(ps[x] for x in ('l','lp'))
    e_in=massless(l,int(helicity));e_out=massless(lp,int(helicity))
    current=np.array([np.vdot(e_out,G[0]@G[mu]@e_in) for mu in range(4)])
    pol=np.eye(4)[1:3]
    amplitudes=np.array([_diagrams(ps,m,current,e).sum(axis=0) for e in pol]).reshape(2,-1)
    return amplitudes@amplitudes.conj().T


def compare(**inputs):
    """Independent physical spinor amplitudes versus corrected trace at h."""
    from gluon_born_response import evaluate
    left=direct(**inputs)
    result=evaluate(**inputs)
    right=np.array(result['B_real'])+1j*np.array(result['B_imag'])
    residual=float(np.max(np.abs(left-right)))
    scale=float(max(np.max(np.abs(right)),1e-30))
    if residual>1e-9+1e-10*scale:
        raise AssertionError(f'direct/trace mismatch: {inputs}, {left}, {right}')
    return {'max_abs':residual,'max_relative':residual/scale,
            'stokes':result['Stokes']}


def compare_source_label(**inputs):
    """Historical v1 reproduction: old lambda label = -physical spinor h."""
    from gluon_born_response import evaluate_legacy_source_label
    source_lambda=inputs.get('helicity',1.)
    if source_lambda not in (-1.,0.,1.):raise ValueError('supported polarization is -1, 0, +1')
    left=direct(**dict(inputs,helicity=-source_lambda))
    report=evaluate_legacy_source_label(**inputs)
    right=np.array(report['B_real'])+1j*np.array(report['B_imag'])
    residual=float(np.max(np.abs(left-right)))
    scale=float(np.max(np.abs(right)))
    if residual>1e-9+1e-10*scale:
        raise AssertionError(f'direct/trace mismatch after explicit helicity mapping: {inputs}')
    return {'max_abs':residual,'max_relative':residual/scale,
            'source_lambda':source_lambda,'spinor_helicity':-source_lambda,
            'B_real':left.real.tolist(),'B_imag':left.imag.tolist()}


def spinor_and_ward_checks():
    ps=fourvectors();m=1.5
    for name,anti in (('p1',False),('p2',True)):
        p=ps[name];spinors=massive(p,m,anti)
        projector=sum((np.outer(a,a.conj())@G[0] for a in spinors),np.zeros((4,4),complex))
        expected=slash(p)+(-m if anti else m)*np.eye(4)
        assert np.max(np.abs(projector-expected))<1e-13
    for h in (-1,1):
        a=massless(ps['l'],h)
        assert np.linalg.norm(slash(ps['l'])@a)<1e-13
        assert abs(np.vdot(a,a)-2*ps['l'][0])<1e-13
    lepton=massless(ps['l'],-1);out=massless(ps['lp'],-1)
    current=np.array([np.vdot(out,G[0]@G[i]@lepton) for i in range(4)])
    ey=np.eye(4)[2]
    photon=_diagrams(ps,m,ps['q'],ey)
    gluon=_diagrams(ps,m,current,ps['k'])
    pq=np.max(np.abs(photon.sum(axis=0)))
    pk=np.max(np.abs(gluon.sum(axis=0)))
    reference=max(np.max(np.abs(_diagrams(ps,m,current,np.eye(4)[i]))) for i in (1,2))
    qscale=np.linalg.norm(ps['q']);kscale=np.linalg.norm(ps['k'])
    relative_q=float(pq/(reference*qscale))
    relative_k=float(pk/(reference*kscale))
    assert pq<1e-12 and pk<1e-11
    assert relative_q<1e-12 and relative_k<1e-12
    assert np.max(np.abs(photon[0]))>1 and np.max(np.abs(gluon[0]))>1
    return {'massive_spin_sums':2,'massless_helicities':2,
            'photon_ward_abs':float(pq),'gluon_ward_abs':float(pk),
            'photon_ward_relative':relative_q,'gluon_ward_relative':relative_k,
            'photon_per_diagram_abs':float(np.max(np.abs(photon[0]))),
            'gluon_per_diagram_abs':float(np.max(np.abs(gluon[0]))),
            'ward_reference_amplitude':float(reference),
            'per_diagram_nonzero':True}
