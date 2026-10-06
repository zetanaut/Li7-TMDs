"""Test-local wrong conventions must fail their named scientific diagnostic."""
import sys
from pathlib import Path
import unittest
from unittest.mock import patch
import sympy as s

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from response_fixtures import SIDIS_ROWS,DY_ROWS
from process_responses import derived_expression,expected_expression,H,T
from process_convolutions import momentum_sign_control,exact_integral,quadrature,FIXTURE
from process_normalization import qed_current_checks,recoil_angle
from process_observables import flavor_checks,target_checks,normalized_sine_moment,sidis_flavor_checks
from process_reversal import Link,Coefficient,reverse_density,time_reversal_unitary,derived_coefficient_sign
from spin_foundations import I4
from correlator_foundations import catalogue


class ProcessFailureTests(unittest.TestCase):
    def test_collins_sign_and_high_branch(self):
        for row in (SIDIS_ROWS[1],SIDIS_ROWS[-1]):
            good=expected_expression('SIDIS',row)
            self.assertEqual(s.cancel(derived_expression('SIDIS',row)-good),0)
            self.assertNotEqual(s.cancel(-derived_expression('SIDIS',row)-good),0)
        row=SIDIS_ROWS[-1]
        self.assertNotEqual(s.cancel(2*derived_expression('SIDIS',row)-expected_expression('SIDIS',row)),0)

    def test_target_factor_and_azimuth_conversion(self):
        row=SIDIS_ROWS[24]
        expression=derived_expression('SIDIS',row)
        self.assertNotEqual(s.cancel(s.Rational(3,5)*expression-expression),0)
        self.assertNotEqual(s.cancel(s.Rational(3,5)**2*expression-s.Rational(3,5)*expression),0)
        # Rotating only the hadron axis changes the target phase in a mixed state.
        self.assertNotEqual(s.cancel(expression.subs(T,T*H)-expression),0)

    def test_momentum_sign_and_jacobian(self):
        control=momentum_sign_control()
        self.assertTrue(control['unweighted_same'])
        self.assertGreater(control['odd_residual'],0)
        row=SIDIS_ROWS[-1]
        actual=float(s.N(exact_integral('SIDIS',row)))
        with patch('process_convolutions.sigma_p',return_value=1):
            wrong=quadrature('SIDIS',row,12,'cartesian')
        self.assertGreater(abs(actual-wrong),1e-5)
        z=FIXTURE['z']
        self.assertNotEqual(1/z**2,1)
        sidis=sidis_flavor_checks()
        self.assertNotEqual(sidis['with_x'],sidis['without_x'])

    def test_antiquark_beam_and_exchange(self):
        current=qed_current_checks()
        self.assertEqual(current['physical_longitudinal_analyzing_sign'],-1)
        self.assertNotEqual(current['physical_longitudinal_analyzing_sign'],1)
        row=DY_ROWS[-1]
        self.assertNotEqual(s.cancel(-derived_expression('DY',row)-expected_expression('DY',row)),0)
        flavor=flavor_checks()
        self.assertNotAlmostEqual(flavor['combined'],float(s.Rational(7,20))*flavor['combined'])
        self.assertNotEqual(flavor['second_ordering'],0)
        self.assertNotAlmostEqual(flavor['combined'],flavor['first_ordering'])
        self.assertNotAlmostEqual(flavor['combined'],flavor['first_ordering']+2*flavor['second_ordering'])

    def test_antiunitary_and_links(self):
        U=time_reversal_unitary()
        A=s.Matrix([[1,s.I,0,1],[s.I,2,1,0],[0,1,1,s.I],[1,0,s.I,1]])
        rho=A*A.H/s.trace(A*A.H)
        self.assertNotEqual(reverse_density(rho),U*rho*U.H)
        odd=next(label for label in catalogue('quark') if label.K==3 and label.channel=='f' and label.n%2)
        self.assertNotEqual(derived_coefficient_sign(odd),derived_coefficient_sign(odd)*(-1)**odd.n)
        link=Link('gluon',('+','-'))
        self.assertNotEqual(link.reverse(),Link('gluon',('-','-')))
        self.assertNotEqual(link.class_name,Link('gluon',('+','+')).class_name)
        coefficient=Coefficient('quark','f',3,1,1,'A','u',False,
                                Link('quark',('+',)),'2 GeV','4 GeV2','scheme-A')
        self.assertTrue(coefficient.comparable(coefficient.reverse()))
        from dataclasses import replace
        with self.assertRaisesRegex(ValueError,'mismatched'):
            coefficient.comparable(replace(coefficient.reverse(),mu='3 GeV'))
        with self.assertRaises(NotImplementedError):
            Link('gluon',('+','+'),'loop_trace')

    def test_transverse_prep_and_denominator(self):
        moments=target_checks()['transverse_alpha']
        self.assertNotEqual(moments[1],'0')
        with self.assertRaisesRegex(ValueError,'undefined reference normalization'):
            normalized_sine_moment(1,0,s.Rational(1,3))
        with self.assertRaisesRegex(ValueError,'recoil azimuth undefined'):
            recoil_angle(0,0)

    def test_whole_asymmetry_sign(self):
        # Imposed target f_DY=-f_SIDIS with distinct partners and kernels.
        target=s.Rational(2,7);sidis_companion=s.Rational(3,5)
        dy_companion=s.Rational(11,10)
        sidis=target*sidis_companion
        dy=-target*dy_companion
        self.assertNotEqual(dy,-sidis)


if __name__=='__main__':unittest.main()
