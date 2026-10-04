"""Leptonic kinematics, massless on-shell QED current, and Born measures."""
from __future__ import annotations
import math
import numpy as np
import sympy as s
from process_dirac import gamma,basis


def sidis_normalization():
    x,z,Q,y,alpha,lam=s.symbols('x z Q y alpha lam',positive=True)
    eps=(1-y)/(1-y+y*y/2)
    dep=y*(2-y)/(2*(1-y+y*y/2))
    # Breit frame: q=(0,0,0,-Q), P points +z, y=Q/(E+Q/2).
    energy=Q/y-Q/2
    lT2=s.simplify(energy**2-Q**2/4)
    assert s.simplify(lT2-Q*Q*(1-y)/y**2)==0
    assert s.simplify(2*Q*energy-Q*Q*(2-y)/y)==0
    L=s.Matrix([[4*lT2+Q*Q,-s.I*lam*Q*Q*(2-y)/y],
                [s.I*lam*Q*Q*(2-y)/y,Q*Q]])
    normalized=(s.eye(2)+eps*s.diag(1,-1)-s.I*lam*dep*s.Matrix([[0,1],[-1,0]]))/2
    assert s.simplify(L/s.trace(L)-normalized)==s.zeros(2)
    assert s.simplify(1-eps**2-dep**2)==0
    pref=s.simplify(alpha**2*y/(8*z*Q**4)*(2*z/x)*s.trace(L))
    assert s.simplify(pref-alpha**2*(1-y+y*y/2)/(x*y*Q*Q))==0
    # j=-zp and d²p=d²j/z²; this is distinct from changing d²P_h.
    j1,j2,p1,p2=s.symbols('j1 j2 p1 p2',real=True)
    jac=s.det(s.Matrix([[-z,0],[0,-z]]))
    assert jac==z*z
    assert s.det(s.Matrix([[-1/z,0],[0,-1/z]]))==1/z**2
    return {'source_labels':['eq:SIDISdelta','eq:Lnorm','eq:depolarization','eq:prefactor'],
            'leptonic_normalization':'Breit-frame on-shell l,lprime; L / Tr_transverse L',
            'prefactor':'alpha_em**2*(1-y+y**2/2)/(x*y*Q**2)',
            'fragmentation_jacobian':'d2p=d2j/z**2',
            'domain':'0<y<1, x,z,Q>0; y=0 is excluded as an input'}


_G=tuple(np.array(g,dtype=complex) for g in gamma())
_G0=_G[0]


def helicity_spinor(h,theta,phi,antiparticle=False,energy=1.):
    """Massless solutions normalized so sums equal slash(p), h=±1."""
    if h not in (-1,1) or energy<=0:raise ValueError('spinor domain')
    c=math.cos(theta/2);ss=math.sin(theta/2)
    def chi(sign):
        return (np.array([c,np.exp(1j*phi)*ss]) if sign==1
                else np.array([-np.exp(-1j*phi)*ss,c]))
    q=chi(-h if antiparticle else h)
    return math.sqrt(energy)*np.r_[-h*q,q] if antiparticle else math.sqrt(energy)*np.r_[q,h*q]


def slash(theta,phi,energy=1.):
    n=np.array([math.sin(theta)*math.cos(phi),math.sin(theta)*math.sin(phi),math.cos(theta)])
    return energy*(_G0-sum(n[i]*_G[i+1] for i in range(3)))


def qed_amplitude(ha,hb,he,hp,theta,phi):
    """Coupling-stripped q qbar -> l- l+; Q²=4 in the test units."""
    qa=helicity_spinor(ha,0,0)
    qb=helicity_spinor(hb,math.pi,0,True)
    lm=helicity_spinor(he,theta,phi)
    lp=helicity_spinor(hp,math.pi-theta,phi+math.pi,True)
    jq=np.array([qb.conj()@_G0@g@qa for g in _G])
    jl=np.array([lm.conj()@_G0@g@lp for g in _G])
    return (jq[0]*jl[0]-jq[1:]@jl[1:])/4


def qed_current_checks():
    cases=[]
    _,g5,_,gm,_=basis()
    axial=np.array(gm*g5,dtype=complex)
    number=np.array(gm,dtype=complex)
    for h in (-1,1):
        for anti,expected_sign in ((False,1),(True,-1)):
            psi=helicity_spinor(h,math.pi,0,anti)
            density=np.outer(psi,psi.conj())@_G0
            f=np.trace(density@number)/2
            g=np.trace(density@axial)/2
            assert abs(g-expected_sign*h*f)<2e-14
    for theta,phi in ((.7,.4),(1.1,.8),(1.4,1.3)):
        lepton_tensor=np.zeros((2,2),complex)
        for he in (-1,1):
            for hp in (-1,1):
                lm=helicity_spinor(he,theta,phi)
                lp=helicity_spinor(hp,math.pi-theta,phi+math.pi,True)
                j=np.array([lm.conj()@_G0@_G[i]@lp for i in (1,2)])
                lepton_tensor+=np.outer(j,j.conj())
        direction=np.array([math.cos(phi),math.sin(phi)])
        expected_tensor=8*(np.eye(2)-math.sin(theta)**2*np.outer(direction,direction))
        assert np.max(np.abs(lepton_tensor-expected_tensor))<3e-14
        for antiparticle in (False,True):
            mat=sum(np.outer(helicity_spinor(h,theta,phi,antiparticle),
                             helicity_spinor(h,theta,phi,antiparticle).conj())@_G0
                    for h in (-1,1))
            assert np.max(np.abs(mat-slash(theta,phi)))<2e-14
        powers={(ha,hb):sum(abs(qed_amplitude(ha,hb,he,hp,theta,phi))**2
                             for he in (-1,1) for hp in (-1,1))
                for ha in (-1,1) for hb in (-1,1)}
        target=1+math.cos(theta)**2
        assert abs(sum(powers.values())/4-target)<3e-14
        assert max(abs(powers[h,h]) for h in (-1,1))<3e-14
        assert max(abs(powers[h,-h]-2*target) for h in (-1,1))<3e-14
        # Physical helicities are measured along each incoming momentum.
        # Correlation is -1; the formal trace has +G_A barG_B before mapping.
        corr=sum(ha*hb*value for (ha,hb),value in powers.items())/sum(powers.values())
        assert abs(corr+1)<3e-14
        # Coherent transverse preparations probe the two-unit lepton phase.
        transverse=[]
        for sign in (1,-1):
            qa=(helicity_spinor(1,0,0)+helicity_spinor(-1,0,0))/math.sqrt(2)
            qb=(helicity_spinor(1,math.pi,0,True)+
                sign*helicity_spinor(-1,math.pi,0,True))/math.sqrt(2)
            rate=0.
            for he in (-1,1):
                for hp in (-1,1):
                    lm=helicity_spinor(he,theta,phi)
                    lp=helicity_spinor(hp,math.pi-theta,phi+math.pi,True)
                    jq=np.array([qb.conj()@_G0@g@qa for g in _G])
                    jl=np.array([lm.conj()@_G0@g@lp for g in _G])
                    rate+=abs((jq[0]*jl[0]-jq[1:]@jl[1:])/4)**2
            assert abs(rate-(target+sign*math.sin(theta)**2*math.cos(2*phi)))<3e-14
            transverse.append(float(rate))
        cases.append({'theta':theta,'phi':phi,'unpolarized':float(sum(powers.values())/4),
                      'expected':target,'helicity_correlation':float(corr),
                      'transverse_rates':transverse,
                      'lepton_tensor_residual':float(np.max(np.abs(lepton_tensor-expected_tensor)))})
    return {'cases':cases,'spin_sums':'sum uu_bar = sum vv_bar = slash(p)',
            'normalization':'spin averaged |M/e²|²=1+cos²(theta)',
            'annihilation_antiquark_axial_ratio':'bar G / bar F = - physical helicity',
            'physical_longitudinal_analyzing_sign':-1}


def dy_normalization():
    alpha,shat,Q,Nc=s.symbols('alpha s Q N_c',positive=True)
    theta=s.symbols('theta',real=True)
    costheta=s.symbols('costheta',real=True)
    angular=s.integrate(1+costheta**2,(costheta,-1,1))*2*s.pi
    assert angular==16*s.pi/3
    # Two-body massless phase space: d sigma_hat/dOmega=|M|²/(64 pi² Q²).
    e4=(4*s.pi*alpha)**2
    partonic=s.simplify(e4/(64*s.pi**2*Q**2)*(1+costheta**2))
    assert partonic==alpha**2*(1+costheta**2)/(4*Q**2)
    hadronic=s.simplify(partonic/(shat*Nc))
    assert hadronic==alpha**2*(1+costheta**2)/(4*shat*Q**2*Nc)
    assert s.simplify(hadronic.subs(Nc,3).integrate((costheta,-1,1))*2*s.pi-
                      4*s.pi*alpha**2/(9*shat*Q**2))==0
    return {'source_label':'eq:DYnormalization','angular_integral':str(angular),
            'color_factor':'1/N_c','hadronic_prefactor':'alpha_em**2/(4*s*Q**2)',
            'partonic_route':'on-shell spinor amplitude, two-body flux/phase space',
            'sidis_x_factor':False,'fragmentation_jacobian':False}


def recoil_angle(qx,qy):
    """Collins-Soper transverse X is recoil aligned only when q_T>0."""
    if qx==0 and qy==0:raise ValueError('recoil azimuth undefined at q_T=0')
    return math.atan2(qy,qx)
