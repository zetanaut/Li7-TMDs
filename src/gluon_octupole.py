"""Physical-cone gluon rates and the fourteen source projections."""
from __future__ import annotations
import sympy as s
from spin_foundations import cartesian,I4
from spin_state_foundations import x3
from correlator_foundations import catalogue,cartesian_value,target_component_operators
from gluon_angular import MODES


def cone_states(psi,O=s.Rational(1,3),c=s.Rational(3,5)):
    if O<=0 or O>s.Rational(1,2) or abs(c)>1:raise ValueError('physical cone requires 0<O<=1/2, |cos(theta)|<=1')
    J,*_=cartesian();sin=s.sqrt(1-c*c)
    n=(sin*s.cos(psi),sin*s.sin(psi),c)
    X=x3(n,J)
    return {eta:I4/4+eta*O*s.sqrt(5)/3*X for eta in (-1,1)}


def state_moments():
    psi=s.symbols('psi',real=True)
    states=cone_states(psi)
    ops=target_component_operators()
    return {(eta,m,part):s.trigsimp(s.trace(rho*ops[3,m,part]))
            for eta,rho in states.items() for m in range(4)
            for part in range(1 if m==0 else 2)}


def term_rate(label,eta,beam,phi,psi,beta,t,moments):
    M=s.symbols('M',positive=True)
    kx=M*t*s.cos(phi);ky=M*t*s.sin(phi)
    F=G=C=S=0
    for part in range(1 if label.m==0 else 2):
        w=moments[eta,label.m,part]
        f,g,c,z=cartesian_value(label,part,kx,ky,M)
        F+=w*f;G+=w*g;C+=w*c;S+=w*z
    return s.expand(F+beam*G+s.cos(2*beta)*C+s.sin(2*beta)*S)


def verify_physical_rates():
    phi,psi,beta,t=s.symbols('phi psi beta t',real=True)
    moments=state_moments()
    O=s.Rational(1,3)
    alpha=(-s.Rational(9,25),s.Rational(8,25),
           s.Rational(48,25),s.Rational(32,25))
    for eta in (-1,1):
        for m in range(4):
            assert s.trigsimp(moments[eta,m,0]-eta*O*alpha[m]*s.cos(m*psi))==0
            if m:assert s.trigsimp(moments[eta,m,1]-eta*O*alpha[m]*s.sin(m*psi))==0
    labels={z.id():z for z in catalogue('gluon') if z.K==3}
    out={}
    for row in MODES:
        name,eps,p,q,kind,shift=row
        channel=name[0];m=int(name.split('.')[1][1]);n=(m if channel!='h' else int(name.split('.')[2]))
        target=next(z for z in labels.values() if (z.channel,z.m,z.n)==(channel,m,n))
        rate={(beam,eta):term_rate(target,eta,beam,phi,psi,beta,t,moments)
              for beam in (-1,1) for eta in (-1,1)}
        reduced=s.expand(sum(eta*(beam if eps else 1)*rate[beam,eta]
                             for beam in (-1,1) for eta in (-1,1))/(4*O))
        phase=p*phi+q*psi+shift*beta
        trig=1 if kind=='one' else s.sin(phase) if kind=='sin' else s.cos(phase)
        factor=(1 if channel!='h' else {0:s.Rational(1,2),1:s.Rational(1,2),2:1,3:1}[m]
                if n==abs(m-2) else s.Rational(1,4))
        expected=alpha[m]*factor*t**n*trig
        residual=s.simplify(s.expand((reduced-expected).rewrite(s.exp),power_exp=True))
        if residual!=0:raise AssertionError(f'physical cone {name}: {residual}')
        out[name]={'source_label':'tab:oct_fourier','target':target.id(),
                   'signed_alpha':str(alpha[m]),'orbital_rank':n,
                   'source_factor':str(factor),'four_spin_rates':4}
    return out


def projection_checks():
    phi,beta,psi=s.symbols('phi beta psi',real=True)
    # Purely longitudinal preparation has alpha_0=1 and no transverse
    # moments.  Its two independent response kernels remain beam odd/even.
    c=s.Integer(1);z=s.Integer(0)
    longitudinal=(s.legendre(3,c),(5*c*c-1)*z/2,5*c*z*z,s.Rational(5,2)*z**3)
    assert longitudinal==(1,0,0,0)
    transverse=(0,-s.Rational(1,2),0,s.Rational(5,2))
    psi_ref=s.Rational(2,5)
    angles=[2*s.pi*j/8 for j in range(8)]
    projector=lambda f:s.simplify(sum(s.cos(3*(psi_ref-x))*f.subs(psi,x)
                                     for x in angles)/4)
    llt=transverse[1]*s.sin(phi-psi)
    if s.trigsimp(projector(llt))!=0:raise AssertionError('LLT leaked into TTT projector')
    modes={
      'f33':s.sin(3*phi-3*psi),
      'g33':s.cos(3*phi-3*psi),
      'h33_low':s.sin(3*psi-phi-2*beta),
      'h33_high':s.sin(5*phi-3*psi-2*beta),
    }
    for name,mode in modes.items():
        observed=projector(transverse[3]*mode)
        expected=transverse[3]*mode.subs(psi,psi_ref)
        if s.simplify(s.expand((observed-expected).rewrite(s.exp),power_exp=True))!=0:
            raise AssertionError(f'TTT projection: {name}')
    x=s.symbols('x',real=True)
    norms={n:s.simplify(2/(5*s.pi)*s.integrate(s.Rational(5,2)*s.sin(n*x)**2,
                                                    (x,0,2*s.pi))) for n in (1,3,5)}
    if any(value!=1 for value in norms.values()):raise AssertionError('TTT recoil moment normalization')
    return {'longitudinal_alpha':list(map(str,longitudinal)),
            'transverse_alpha':list(map(str,transverse)),
            'eight_angle_TTT_modes':list(modes),'LLT_leakage':'0',
            'recoil_moment_unit_norms':{str(k):str(v) for k,v in norms.items()},
            'source_labels':['eq:TTTprep_projector','eq:TTTgluon_moments','eq:longitudinal_gluon_observables']}
