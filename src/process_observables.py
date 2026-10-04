"""Physical target preparations and unequal-flavor synthetic convolutions."""
from __future__ import annotations
import sympy as s
from spin_foundations import cartesian,I4
from spin_state_foundations import x3
from correlator_foundations import target_component_operators
from process_convolutions import exact_integral,FIXTURE
from response_fixtures import SIDIS_ROWS,DY_ROWS


def cone_moments(cosine,octupole=s.Rational(1,3),psi=0):
    if not (0<=cosine<=1):raise ValueError('cone domain in this fixture')
    sine=s.sqrt(1-cosine*cosine)
    J,*_=cartesian()
    n=(sine*s.cos(psi),sine*s.sin(psi),cosine)
    X=x3(n,J)
    rho_plus=I4/4+octupole*s.sqrt(5)/3*X
    rho_minus=I4/4-octupole*s.sqrt(5)/3*X
    ops=target_component_operators()
    alpha=(s.legendre(3,cosine),
           (5*cosine*cosine-1)*sine/2,
           5*cosine*sine*sine,
           s.Rational(5,2)*sine**3)
    for eta,rho in ((1,rho_plus),(-1,rho_minus)):
        assert s.trace(rho)==1
        assert min(s.N(v) for v in rho.eigenvals())>=0
        for m in range(4):
            observed=s.simplify(s.trace(rho*ops[3,m,0]))
            assert s.simplify(observed-eta*octupole*alpha[m]*s.cos(m*psi))==0
            if m:
                observed_im=s.simplify(s.trace(rho*ops[3,m,1]))
                assert s.simplify(observed_im-eta*octupole*alpha[m]*s.sin(m*psi))==0
    return tuple(s.simplify(value) for value in alpha)


def normalized_sine_moment(numerator,reference,octupole):
    if reference==0 or octupole==0:raise ValueError('undefined reference normalization')
    return 2*numerator/(octupole*reference)


def target_checks():
    generic=cone_moments(s.Rational(3,5))
    transverse=cone_moments(s.Integer(0))
    assert generic==(-s.Rational(9,25),s.Rational(8,25),
                     s.Rational(48,25),s.Rational(32,25))
    assert transverse==(0,-s.Rational(1,2),0,s.Rational(5,2))
    assert normalized_sine_moment(s.Rational(3,10),s.Rational(3,2),s.Rational(2,5))==1
    # A beam helicity difference for one state sees residual vector g_10.
    vector_contamination=s.Rational(1,7)
    octupole_signal=s.Rational(1,3)
    assert vector_contamination+octupole_signal!=octupole_signal
    return {'generic_alpha':list(map(str,generic)),
            'transverse_alpha':list(map(str,transverse)),
            'sine_projection_factor':2,'vector_contamination_rejected':True,
            'source_labels':['eq:cone_alpha','eq:DYnumbermoment','eq:octdouble']}


def flavor_checks():
    # Two inequivalent flavors; both annihilation orderings have independent
    # widths and amplitudes. xA != xB and all masses/fractions stay nonunit.
    charges={'u':s.Rational(4,9),'d':s.Rational(1,9)}
    xA=s.Rational(7,20);xB=s.Rational(11,50)
    inputs={
      'u':{'Aq':s.Rational(7,5),'Abar':s.Rational(2,5),'Bq':s.Rational(6,5),
           'Bbar':s.Rational(3,5),'waq':s.Rational(7,5),'wabar':s.Rational(9,10),
           'wbq':s.Rational(11,10),'wbbar':s.Rational(8,5)},
      'd':{'Aq':s.Rational(4,5),'Abar':s.Rational(3,10),'Bq':s.Rational(9,10),
           'Bbar':s.Rational(1,2),'waq':s.Rational(6,5),'wabar':s.Rational(13,10),
           'wbq':s.Rational(17,10),'wbbar':s.Rational(4,5)},
    }
    row=DY_ROWS[0];total=0;first=0;second=0
    for flavor,inputs_f in inputs.items():
        for order,amp_a,amp_b,wa,wb in (
          ('first','Aq','Bbar','waq','wbbar'),
          ('second','Abar','Bq','wabar','wbq')):
            fixture=dict(FIXTURE,width_a=inputs_f[wa],width_b=inputs_f[wb])
            # Independent synthetic fraction dependence on the two parents.
            value=(charges[flavor]*(1+xA)*inputs_f[amp_a]*
                   (1+2*xB)*inputs_f[amp_b]*exact_integral('DY',row,fixture)/3)
            if order=='first':first+=value
            else:second+=value
            total+=value
    assert s.simplify(first)!=0 and s.simplify(second)!=0 and s.simplify(first-second)!=0
    assert s.simplify(total-first)==second and s.simplify(total-(first+2*second))==-second
    return {'flavors':list(inputs),'x_A':str(xA),'x_B':str(xB),
            'first_ordering':float(s.N(first)),'second_ordering':float(s.N(second)),
            'combined':float(s.N(total)),'color_factor':'1/3',
            'source_label':'eq:DYconv'}


def sidis_flavor_checks():
    x=s.Rational(7,20)
    row=SIDIS_ROWS[20]  # octupole number response
    data={
        'u':(s.Rational(4,9),s.Rational(6,5),s.Rational(7,10),
             s.Rational(7,5),s.Rational(11,10)),
        'd':(s.Rational(1,9),s.Rational(4,5),s.Rational(9,10),
             s.Rational(13,10),s.Rational(8,5)),
    }
    terms=[]
    for charge,distribution,fragmentation,wa,wb in data.values():
        fixture=dict(FIXTURE,width_a=wa,width_b=wb)
        terms.append(charge*distribution*fragmentation*exact_integral('SIDIS',row,fixture))
    total=x*sum(terms)
    assert total!=0 and total!=sum(terms)
    return {'flavors':list(data),'x_Bj':str(x),'flavor_factor':'x_Bj e_a**2',
            'distribution_and_fragmentation_independent':True,
            'with_x':float(s.N(total)),'without_x':float(s.N(sum(terms))),
            'source_label':'eq:convolution'}


def spin_difference_checks():
    O,N,dep,eps,phi,F,G,H,L=s.symbols('O N dep eps phi F G H L',real=True)
    # Physical +/- preparations; L is a residual vector moment shared by
    # both states. An ordinary beam difference retains L.
    sidis=lambda lam,eta:N*(F+lam*dep*(L+eta*O*G)+eta*O*eps*s.sin(2*phi)*H)
    double=s.expand(sum(lam*eta*sidis(lam,eta)
                        for lam in (-1,1) for eta in (-1,1))/4)
    assert s.simplify(double-N*O*dep*G)==0
    ordinary=s.expand((sidis(1,1)-sidis(-1,1))/2)
    assert s.simplify(ordinary-N*dep*(L+O*G))==0
    theta=s.symbols('theta',real=True)
    dy=lambda lam,eta:N*(1+s.cos(theta)**2)*(F-lam*eta*O*G)
    dy_double=s.expand(sum(lam*eta*dy(lam,eta)
                           for lam in (-1,1) for eta in (-1,1))/4)
    assert s.simplify(dy_double+N*O*(1+s.cos(theta)**2)*G)==0
    return {'sidis_double':'N*O*D_ell*G30',
            'sidis_ordinary_includes_vector':'N*D_ell*(L+O*G30)',
            'dy_double':'-N_DY*O*(1+cos(theta_l)**2)*G30',
            'source_labels':['eq:octdouble','eq:DYlongdouble']}
