"""Synthetic octupole rates from Cartesian covariants, fitted with Fourier modes."""
from __future__ import annotations
import numpy as np
import sympy as s
from correlator_foundations import catalogue,cartesian_value
from gluon_angular import MODES
from gluon_octupole import state_moments
from born_direct import compare_source_label


def analyzer():
    result=compare_source_label()
    B=np.array(result['B_real'])+1j*np.array(result['B_imag'])
    bu=(B[0,0].real+B[1,1].real)/2
    bg=B[0,1].imag
    bc=(B[0,0].real-B[1,1].real)/2
    bs=B[0,1].real
    return bu,bg,np.hypot(bc,bs),np.arctan2(bs,bc)/2


def alpha(c):
    z=np.sqrt(1-c*c)
    return np.array([(5*c**3-3*c)/2,(5*c*c-1)*z/2,5*c*z*z,2.5*z**3])


def _row_parts(row):
    name,eps,p,q,kind,shift=row
    ch=name[0];m=int(name.split('.')[1][1]);n=m if ch!='h' else int(name.split('.')[2])
    low=n==abs(m-2)
    factor=1 if ch!='h' else ({0:.5,1:.5,2:1.,3:1.}[m] if low else .25)
    return ch,m,n,factor


def matrix(samples,bu,bg,bl,beta,t,c,O,acceptance=None):
    a=alpha(c)
    X=np.empty((len(samples),14))
    for j,(name,eps,p,q,kind,shift) in enumerate(MODES):
        ch,m,n,factor=_row_parts(MODES[j])
        hard={'f':bu,'g':bg,'h':bl}[ch]
        for i,(beam,phi,psi) in enumerate(samples):
            angle=p*phi+q*psi+shift*beta
            wave=1 if kind=='one' else np.sin(angle) if kind=='sin' else np.cos(angle)
            X[i,j]=O*a[m]*factor*t**n*hard*beam**eps*wave
    if acceptance is not None:
        X*=acceptance[:,None]
    return X


def _cartesian_functions():
    phi=s.symbols('phi',real=True)
    t=s.symbols('t',positive=True)
    labels={(x.channel,x.m,x.n):x for x in catalogue('gluon') if x.K==3}
    functions=[]
    for row in MODES:
        label=labels[_row_parts(row)[:3]]
        functions.append([s.lambdify((phi,t),cartesian_value(label,part,t*s.cos(phi),
                                    t*s.sin(phi),s.Integer(1)),'numpy')
                          for part in range(1 if label.m==0 else 2)])
    return functions


def forward(samples,coefficients,bu,bg,bl,beta,t,c,O):
    """Calculate rates via actual Cartesian traces, independent of fit matrix."""
    if c!=.6 or O!=1/3:
        raise ValueError('forward density fixture uses cos(theta)=3/5 and O=1/3')
    funcs=_cartesian_functions()
    symbolic_moments=state_moments()
    psi_symbol=s.symbols('psi',real=True)
    state={key:s.lambdify(psi_symbol,value,'numpy') for key,value in symbolic_moments.items()}
    rates=np.empty((len(samples),2))
    for i,(beam,phi,psi) in enumerate(samples):
        for k,eta in enumerate((-1,1)):
            value=bu
            for j,row in enumerate(MODES):
                ch,m,n,_=_row_parts(row)
                for part,fun in enumerate(funcs[j]):
                    F,G,C,S=(float(x) for x in fun(phi,t))
                    moment=float(state[eta,m,part](psi))
                    value+=coefficients[j]*moment*(bu*F+beam*bg*G+
                                                    bl*np.cos(2*beta)*C+bl*np.sin(2*beta)*S)
            rates[i,k]=value
    return rates


def _grid():
    return [(beam,2*np.pi*i/11,2*np.pi*j/8)
            for beam in (-1,1) for i in range(11) for j in range(8)]


def run_reconstruction():
    bu,bg,bl,beta=analyzer()
    t=.6;c=.6;O=1/3
    x_bj=.05;M_pair2=25.;Q2=4.
    x_g=x_bj*(1+M_pair2/Q2)
    if not 0<x_g<1:raise AssertionError('synthetic hadronic fraction outside support')
    samples=_grid()
    rng=np.random.default_rng(314159)
    coefficients=np.array([.015*(-1)**i*(1+.1*i) for i in range(14)])
    rates=forward(samples,coefficients,bu,bg,bl,beta,t,c,O)
    if rates.min()<=0:raise AssertionError('synthetic rates must be nonnegative')
    X=matrix(samples,bu,bg,bl,beta,t,c,O)
    observed=(rates[:,1]-rates[:,0])/2
    fitted=np.linalg.lstsq(X,observed,rcond=None)[0]
    residual=float(np.max(np.abs(fitted-coefficients)))
    if residual>1e-10:raise AssertionError('Cartesian synthetic recovery failed')
    moment_rows=[]
    denominator=float(np.sum(rates))
    for j,row in enumerate(MODES):
        name,eps,p,q,kind,shift=row
        values=[]
        for beam,phi,psi in samples:
            angle=p*phi+q*psi+shift*beta
            wave=1 if kind=='one' else np.sin(angle) if kind=='sin' else np.cos(angle)
            values.append(beam**eps*wave)
        Z=1 if kind=='one' else 2
        measured=Z/O*float(np.dot(values,rates[:,1]-rates[:,0]))/denominator
        ch,m,n,factor=_row_parts(row)
        expected=coefficients[j]*alpha(c)[m]*factor*t**n*{'f':bu,'g':bg,'h':bl}[ch]/bu
        if abs(measured-expected)>1e-12:raise AssertionError('normalized four-rate moment failed')
        moment_rows.append({'id':name,'Z':Z,'measured':measured,'expected_ratio':float(expected)})
    acceptance=np.array([1+.2*np.cos(phi)+.1*np.sin(psi)+.05*np.cos(phi+psi)
                         for beam,phi,psi in samples])
    assert acceptance.min()>.6
    folded=acceptance*observed
    Xweighted=matrix(samples,bu,bg,bl,beta,t,c,O,acceptance)
    if np.linalg.matrix_rank(Xweighted)!=14:raise AssertionError('accepted response lost rank')
    corrected=np.linalg.lstsq(Xweighted,folded,rcond=None)[0]
    naive=np.linalg.lstsq(X,folded,rcond=None)[0]
    if np.max(np.abs(corrected-coefficients))>1e-10:raise AssertionError('acceptance recovery failed')
    if np.max(np.abs(naive-coefficients))<1e-4:raise AssertionError('acceptance leakage undetected')
    singular=np.linalg.svd(X,compute_uv=False)
    normalized=X/np.linalg.norm(X,axis=0)
    rank_losses={
      'bG_zero':np.linalg.matrix_rank(matrix(samples,bu,0,bl,beta,t,c,O)),
      'bLin_zero':np.linalg.matrix_rank(matrix(samples,bu,bg,0,beta,t,c,O)),
      'alpha1_zero':np.linalg.matrix_rank(matrix(samples,bu,bg,bl,beta,t,1/np.sqrt(5),O)),
      'beam_averaged':np.linalg.matrix_rank((X[:88]+X[88:])/2),
      'fixed_psi':np.linalg.matrix_rank(matrix([(b,p,0.) for b,p,_ in samples],bu,bg,bl,beta,t,c,O)),
      'aliased_4x4':np.linalg.matrix_rank(matrix([(b,2*np.pi*i/4,2*np.pi*j/4)
                  for b in (-1,1) for i in range(4) for j in range(4)],bu,bg,bl,beta,t,c,O)),
      'O_zero':np.linalg.matrix_rank(matrix(samples,bu,bg,bl,beta,t,c,0)),
      't_zero':np.linalg.matrix_rank(matrix(samples,bu,bg,bl,beta,0,c,O)),
    }
    return {'cases':len(samples),'coefficients':coefficients.tolist(),
            'max_recovery_error':residual,'minimum_rate':float(rates.min()),
            'acceptance_min':float(acceptance.min()),
            'naive_acceptance_error':float(np.max(np.abs(naive-coefficients))),
            'corrected_acceptance_error':float(np.max(np.abs(corrected-coefficients))),
            'rank':int(np.linalg.matrix_rank(X)),'rank_losses':{k:int(v) for k,v in rank_losses.items()},
            'condition_raw':float(singular[0]/singular[-1]),
            'condition_column_normalized':float(np.linalg.cond(normalized)),
            'condition_known_acceptance_raw':float(np.linalg.cond(Xweighted)),
            'normalized_moments':moment_rows,
            'analyzer':{'b_U':bu,'b_G':bg,'b_lin':bl,'phi_B':beta},
            'fraction_relation':{'x_Bj':x_bj,'M_pair2_GeV2':M_pair2,
                                 'Q2_GeV2':Q2,'x_g':x_g,
                                 'hard_event':'collinear Born; recoil t is separate leading-power TMD input'},
            'source_labels':['eq:gluon_spin_cross','tab:oct_fourier','eq:gluon_measured_moment']}
