"""Independent Born reference, grid, scan, normalization, and domains."""
from __future__ import annotations
import numpy as np
import mpmath as mp
from born_direct import compare,direct,fourvectors,direct_vectors,product4
from born_precision import evaluate_strings
from gluon_born_response import evaluate

BASE={'sqrt_s':5.,'Q2':4.,'mass':1.5,'theta':.8,'phi':.4,'lepton_energy':10.}


def _matrix(report):return np.array(report['B_real'])+1j*np.array(report['B_imag'])


def raw_matrix_diagnostic(B,scale):
    if B.shape!=(2,2) or not np.all(np.isfinite(B)):
        raise ValueError('nonfinite or malformed hard matrix')
    residual=float(np.max(np.abs(B-B.conj().T)))
    if residual>1e-10*scale:raise ValueError('raw hard matrix is not Hermitian')
    eigen=np.linalg.eigvalsh(B)
    if eigen.min()<-1e-10*scale:raise ValueError('hard matrix is not positive')
    return residual,eigen


def case(inputs):
    mapped=compare(**inputs)
    production=evaluate(**inputs)
    B=_matrix(production)
    ps=fourvectors(**{k:v for k,v in inputs.items() if k!='helicity'})
    momentum=float(np.linalg.norm(ps['l']+ps['k']-ps['lp']-ps['p1']-ps['p2']))
    shells={
      'l':abs(product4(ps['l'],ps['l'])),'lp':abs(product4(ps['lp'],ps['lp'])),
      'k':abs(product4(ps['k'],ps['k'])),
      'q_plus_Q2':abs(product4(ps['q'],ps['q'])+inputs['Q2']),
      'p1_minus_m2':abs(product4(ps['p1'],ps['p1'])-inputs['mass']**2),
      'p2_minus_m2':abs(product4(ps['p2'],ps['p2'])-inputs['mass']**2),
    }
    energy_scale=max(np.linalg.norm(v) for v in ps.values())
    if momentum>1e-10*energy_scale or max(shells.values())>1e-10*energy_scale**2:
        raise AssertionError('independent four-vector kinematics failed')
    direct_B=direct(**inputs)
    polarizations=(np.array([1,1j])/np.sqrt(2),np.array([1,-1j])/np.sqrt(2),
                   np.array([2,1+2j])/3)
    complex_rates=[float(np.real(e@B@e.conj())) for e in polarizations]
    direct_rates=[float(np.real(e@direct_B@e.conj())) for e in polarizations]
    if max(abs(x-y) for x,y in zip(complex_rates,direct_rates))>1e-9:
        raise AssertionError('complex-polarization Born rates disagree')
    raw_herm,_=raw_matrix_diagnostic(B,float(np.max(np.abs(B))))
    return {'convention':production['convention'],
            'convention_digest':production['convention_digest'],
            'reference_source_sha256':production['reference_source_sha256'],
            'inputs':inputs,'B_real':production['B_real'],'B_imag':production['B_imag'],
            'stokes':production['Stokes'],'eigenvalues':production['B_eigenvalues'],
            'direct_residual_abs':mapped['max_abs'],
            'direct_residual_relative':mapped['max_relative'],
            'complex_polarization_rates':complex_rates,
            'raw_hermiticity_abs':raw_herm,
            'relative_wards':production['relative_Ward_residuals'],
            'momentum_residual_GeV':momentum,'mass_shell_residuals_GeV2':shells,
            'kinematic_scale_GeV':float(energy_scale),
            'checks':production['checks']}


def grid():
    cases=[case(dict(BASE,theta=theta,phi=phi,helicity=helicity))
           for theta in (.4,.8,1.7,2.7) for phi in (0.,.4,1.2)
           for helicity in (-1.,0.,1.)]
    expected={(x,y,z) for x in (.4,.8,1.7,2.7) for y in (0.,.4,1.2)
              for z in (-1.,0.,1.)}
    actual={(row['inputs']['theta'],row['inputs']['phi'],row['inputs']['helicity']) for row in cases}
    if len(cases)!=36 or actual!=expected:raise AssertionError('incomplete Born grid')
    return {'convention':'born-current-v2-physical-h',
            'cases':cases,'case_count':len(cases),'independent_direct_count':len(cases),
            'max_direct_residual_abs':max(x['direct_residual_abs'] for x in cases)}


def dense_scan():
    angles=np.linspace(.05,np.pi-.05,161)
    cases=[case(dict(BASE,theta=float(theta),helicity=helicity))
           for theta in angles for helicity in (-1.,0.,1.)]
    expected={(i,h) for i in range(161) for h in (-1.,0.,1.)}
    actual={(i//3,row['inputs']['helicity']) for i,row in enumerate(cases)
            if row['inputs']['theta']==float(angles[i//3])}
    if len(cases)!=483 or actual!=expected:raise AssertionError('incomplete dense scan')
    return {'convention':'born-current-v2-physical-h',
            'angle_array':angles.tolist(),'phi':.4,'helicities':[-1.,0.,1.],
            'cases':cases,'case_count':483,'independent_direct_count':483,
            'max_direct_residual_abs':max(x['direct_residual_abs'] for x in cases),
            'max_direct_residual_relative':max(x['direct_residual_relative'] for x in cases)}


def verify_scan_payload(payload):
    from validation_evidence import REFERENCE_SOURCE_SHA256,digest
    angles=np.linspace(.05,np.pi-.05,161).tolist()
    if payload.get('convention')!='born-current-v2-physical-h' or \
            payload.get('angle_array')!=angles or payload.get('phi')!=.4 or \
            payload.get('helicities')!=[-1.,0.,1.] or payload.get('case_count')!=483 or \
            payload.get('independent_direct_count')!=483:
        raise ValueError('dense scan specification mismatch')
    cases=payload.get('cases')
    if not isinstance(cases,list) or len(cases)!=483:raise ValueError('dense scan cases incomplete')
    for i,row in enumerate(cases):
        expected=dict(BASE,theta=angles[i//3],helicity=(-1.,0.,1.)[i%3])
        if (row.get('inputs')!=expected or row.get('convention')!=payload['convention'] or
                row.get('convention_digest')!=digest(payload['convention']) or
                row.get('reference_source_sha256')!=REFERENCE_SOURCE_SHA256):
            raise ValueError(f'dense scan input tuple or convention mismatch at {i}')
        if row.get('direct_residual_abs',float('inf'))>1e-8 or \
                row.get('direct_residual_relative',float('inf'))>1e-10 or \
                not all(row.get('checks',{}).values()) or \
                len(row.get('B_real',[]))!=2 or len(row.get('B_imag',[]))!=2:
            raise ValueError(f'dense scan diagnostic failed at {i}')
    if abs(payload['max_direct_residual_abs']-max(x['direct_residual_abs'] for x in cases))>1e-20:
        raise ValueError('dense scan residual summary mismatch')
    return True


def precision_set():
    # Exact decimal inputs are passed as strings at each independent precision.
    inputs=(
      {'sqrt_s':'5','Q2':'4','mass':'1.5','theta':'0.8','phi':'0.4','lepton_energy':'10'},
      {'sqrt_s':'5','Q2':'4','mass':'1.5','theta':'1.2','phi':'1.1','lepton_energy':'10'},
      {'sqrt_s':'3.01','Q2':'2','mass':'1.5','theta':'0.7','phi':'0.6','lepton_energy':'6'},
      {'sqrt_s':'5','Q2':'4','mass':'1.5','theta':'0.05','phi':'0.4','lepton_energy':'2.501'},
    )
    out=[]
    for case_inputs in inputs:
        a=evaluate_strings(**case_inputs,spinor_helicity=1,dps=50)
        b=evaluate_strings(**case_inputs,spinor_helicity=1,dps=80)
        p=evaluate(**{k:float(v) for k,v in case_inputs.items()},helicity=1.)
        st=[float(v) for v in b['stokes']]
        error=max(abs(st[i]-p['Stokes'][k]) for i,k in enumerate(('b_U','b_G','b_C','b_S')))
        with mp.workdps(85):
            stable=max(abs(mp.mpf(x)-mp.mpf(y)) for x,y in zip(a['stokes'],b['stokes']))
            stable_text=mp.nstr(stable,12)
            if stable>mp.mpf('1e-40'):raise AssertionError(f'precision instability: {case_inputs}')
        if error>1e-8:raise AssertionError(f'precision comparison: {case_inputs}')
        out.append({'inputs':case_inputs,'dps50':a,'dps80':b,
                    'double_max_abs':error,'precision_difference':stable_text})
    return {'convention':'born-current-v2-physical-h',
            'cases':out,'case_count':len(out),'precision_digits':[50,80],
            'convention':'physical electron helicity h=+1 in spinor and trace'}


def broader_cases():
    specs=(
      dict(sqrt_s=6.,Q2=1.5,mass=1.2,theta=.5,phi=2.,lepton_energy=9.,helicity=-1.),
      dict(sqrt_s=4.2,Q2=6.,mass=2.,theta=2.4,phi=-.7,lepton_energy=7.,helicity=0.),
      dict(sqrt_s=8.,Q2=12.,mass=1.1,theta=.12,phi=1.3,lepton_energy=9.,helicity=1.),
      dict(sqrt_s=5.,Q2=4.,mass=1.5,theta=3.0,phi=2.4,lepton_energy=2.501,helicity=0.),
    )
    return {'cases':[case(x) for x in specs],'case_count':len(specs)}


def normalization_checks():
    # Explicit normalized SU(3) fundamental generators.
    mats=[]
    for i,j in ((0,1),(0,2),(1,2)):
        x=np.zeros((3,3),complex);y=x.copy();x[i,j]=x[j,i]=.5
        y[i,j]=-.5j;y[j,i]=.5j;mats.extend((x,y))
    mats.extend((np.diag([.5,-.5,0]),np.diag([1,1,-2])/(2*np.sqrt(3))))
    TF=sum(np.trace(a@a).real for a in mats)/8
    assert len(mats)==8 and abs(TF-.5)<1e-14
    ps=fourvectors(**BASE)
    flux=2*(ps['l'][0]*ps['k'][0]-np.dot(ps['l'][1:],ps['k'][1:]))
    p=evaluate(**dict(BASE,helicity=1.))
    B=_matrix(p);r=1.7
    scaled=evaluate(**dict(BASE,sqrt_s=r*BASE['sqrt_s'],Q2=r*r*BASE['Q2'],
                           mass=r*BASE['mass'],lepton_energy=r*BASE['lepton_energy'],helicity=1.))
    Bscaled=_matrix(scaled)
    err=float(np.max(np.abs(Bscaled-r*r*B)))
    if err>1e-8:raise AssertionError('reduced Born dimension scaling failed')
    full=B/BASE['Q2']**2
    full_scaled=Bscaled/(r*r*BASE['Q2'])**2
    full_err=float(np.max(np.abs(full_scaled-full/r**2)))
    if full_err>1e-10:raise AssertionError('full Born dimension scaling failed')
    plus=B
    minus=_matrix(evaluate(**dict(BASE,helicity=-1.)))
    partial=_matrix(evaluate(**dict(BASE,helicity=.3)))
    beam_conjugation=float(np.max(np.abs(minus-plus.conj())))
    mixture=float(np.max(np.abs(partial-(.65*plus+.35*minus))))
    if beam_conjugation>1e-10 or mixture>1e-10:
        raise AssertionError('beam polarization mixture failed')
    period=_matrix(evaluate(**dict(BASE,phi=BASE['phi']+2*np.pi,helicity=1.)))
    periodic=float(np.max(np.abs(period-B)))
    if periodic>1e-10:raise AssertionError('azimuthal periodicity failed')
    angle=.73;c=np.cos(angle);ss=np.sin(angle)
    R=np.array([[c,-ss],[ss,c]])
    R3=np.array([[c,-ss,0],[ss,c,0],[0,0,1]])
    rotated={key:np.r_[v[0],R3@v[1:]] for key,v in ps.items()}
    original=direct_vectors(ps,BASE['mass'],-1)
    transformed=direct_vectors(rotated,BASE['mass'],-1)
    rotation=float(np.max(np.abs(transformed-R@original@R.T)))
    if rotation>1e-10:raise AssertionError('active hard-event rotation failed')
    return {'T_F':float(TF),'color_generators':8,'initial_spin_density':'I/2',
            'unpolarized_reduced_rate':float(np.trace(B).real/2),
            'partonic_flux_GeV2':float(flux),'inverse_flux_GeV_minus2':float(1/(2*flux)),
            'reduced_B_scaling_r2_max_abs':err,'full_B_scaling_r_minus2_max_abs':full_err,
            'beam_conjugation_max_abs':beam_conjugation,'partial_beam_mixture_max_abs':mixture,
            'azimuth_periodicity_max_abs':periodic,'active_rotation_max_abs':rotation,
            'source_labels':['eq:BornB','eq:Bornphase']}


def domain_checks():
    invalid=(dict(BASE,Q2=0.),dict(BASE,Q2=-1.),dict(BASE,mass=-1.),
             dict(BASE,sqrt_s=2.9),dict(BASE,lepton_energy=2.5),
             dict(BASE,lepton_energy=2.),dict(BASE,theta=float('nan')))
    outcomes=[]
    for x in invalid:
        try:direct(**dict(x,helicity=1.))
        except ValueError as exc:outcomes.append(str(exc))
        else:raise AssertionError('invalid Born input accepted')
    return {'rejected':outcomes,'count':len(outcomes)}
