"""Scientific mutations for the finite checkpoint-5 checks."""
import sys
from pathlib import Path
import unittest
import sympy as s
import mpmath as mp

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from collinear_limits import (selection,EXPECTED_CANDIDATES,EXPECTED_COLLINEAR,
                              circle_average,X,Y,r,M,operations_and_tail)
from fourier_limits import normalized_component,bessel_quadrature,exact_gaussian
from local_moments import tensor_weights,spin_multiplicities,moment_sign,nuclear_bookkeeping
from conditional_positivity import (fixed_target_checks,source_fixed_matrix,
                                    matrix_algorithm_counterexamples,
                                    source_collinear_matrix,joint_gram_checks,collinear_blocks)
from source_joint_positivity import source_mapping_check

class LimitsFailureTests(unittest.TestCase):
    def test_angular_link_and_zero_momentum(self):
        for species in ('quark','gluon'):
            result=selection(species)
            self.assertEqual(set(result['candidates']),set(EXPECTED_CANDIDATES[species]))
            self.assertEqual(set(result['survivors']),set(EXPECTED_COLLINEAR[species]))
            self.assertEqual(len(set(result['candidates'])-set(result['survivors'])),1)
        self.assertEqual(circle_average(X*X+Y*Y),r*r)
        self.assertNotEqual(circle_average(X*X+Y*Y).subs(r,1),
                            (X*X+Y*Y).subs({X:0,Y:0}))
        self.assertIn('gluon.h.32[0]',selection('gluon')['candidates'])
        self.assertNotIn('gluon.h.32[0]',selection('gluon')['survivors'])

    def test_fourier_phase_mass_and_measure(self):
        exact_gaussian()
        value=bessel_quadrature(1,1,(.45,.63),1.7,2.3,'gaussian')
        wrong_phase=-value
        self.assertGreater(abs(value-wrong_phase),1e-2)
        wrong_mass=value*2.3
        self.assertGreater(abs(value-wrong_mass),1e-2)
        rank5=bessel_quadrature(5,1,(.45,.63),1.7,2.3,'gaussian')
        wrong_stf=rank5*16
        self.assertGreater(abs(rank5-wrong_stf),1e-3)
        radius=(.45**2+.63**2)**.5
        wrong_j0=complex(2*mp.pi*1j*((.45+1j*.63)/radius)*mp.quad(
            lambda k:k*(k/2.3)*mp.besselj(0,k*radius)*
                     mp.exp(-(k/1.7)**2)/(mp.pi*1.7**2),[0,13.6]))
        self.assertGreater(abs(value-wrong_j0),1e-2)
        self.assertEqual(normalized_component(0),1)

    def test_fourier_inverse_and_cutoff(self):
        data=exact_gaussian()
        self.assertEqual(data['inverse_normalization'],'1/(pi*Lambda**2)')
        self.assertNotEqual(s.simplify(1/(s.pi*s.Symbol('Lambda')**2)-
                                       4*s.pi/s.Symbol('Lambda')**2),0)
        tail=operations_and_tail()
        self.assertEqual(tail['unbounded_limit'],'logarithmic')
        scale,R=s.symbols('L R',positive=True)
        self.assertEqual(s.limit(s.log(1+R*R/scale**2),R,s.oo),s.oo)

    def test_local_rank_antiquark_and_nuclear(self):
        self.assertEqual(max(spin_multiplicities(tensor_weights(1,'tensor'))),1)
        self.assertNotIn(2,spin_multiplicities(tensor_weights(1,'tensor')))
        self.assertNotIn(3,spin_multiplicities(tensor_weights(2,'axial')))
        self.assertIn(3,spin_multiplicities(tensor_weights(3,'axial')))
        self.assertEqual(moment_sign('axial',1),1)
        self.assertNotEqual(moment_sign('axial',1),moment_sign('vector',1))
        book=nuclear_bookkeeping()
        self.assertEqual(book['count_only_momentum_counterexample'],'1/2')
        self.assertEqual(book['conditional_total_tensor_witness']['total'],'0')
        self.assertNotEqual(book['conditional_total_tensor_witness']['quark_flavor_1'],'0')
        self.assertIn('all partonic',book['total_tensor_constraint'])

    def test_fixed_trace_and_principal_minors(self):
        fixed_target_checks(); examples=matrix_algorithm_counterexamples()
        self.assertEqual(examples['pairwise_minors'],['7/16']*3)
        self.assertEqual(examples['pairwise_determinant'],'-49/32')
        self.assertEqual(examples['leading_minors'],['0']*3)
        self.assertEqual(examples['singular_positive'],'boundary')
        negative=source_fixed_matrix('quark',-3,1,1,1)
        self.assertEqual(s.trace(negative),-3)
        self.assertGreater((-3)**2,1**2+1**2+1**2)

    def test_joint_transpose_and_block_factor(self):
        self.assertTrue(joint_gram_checks()['partial_transpose_detected'])
        f,h=s.symbols('f h',real=True)
        M=source_collinear_matrix('gluon',f,0,0,0,h,0)
        self.assertNotEqual(M.subs({f:1,h:1}).det(),M.subs({f:1,h:s.sqrt(3)/6}).det())
        self.assertEqual(M.subs({f:1,h:s.sqrt(3)/6}).det(),0)

    def test_source_index_conversion_and_partial_transpose(self):
        mapping=source_mapping_check()
        self.assertEqual(len(mapping['minimal_entries']),2)
        self.assertEqual(mapping['Bell_partial_transpose_negative_eigenvalue'],'-1/2')
        self.assertFalse(mapping['partial_transpose_psd_preserving'])
        self.assertNotEqual(mapping['minimal_entries'][1]['source_twice'],
                            mapping['minimal_entries'][1]['old_auxiliary'])
        identity=source_collinear_matrix('quark',1,0,0,0,0,0)
        self.assertEqual(identity,s.eye(8)/2)
        self.assertNotEqual(identity,s.eye(8)/8)  # inverse-density weights

    def test_collinear_bound_factor_mutation(self):
        h=s.symbols('h',real=True)
        q=collinear_blocks('quark')
        self.assertEqual(len(q['determinants']),3)
        M=source_collinear_matrix('quark',1,0,0,0,h,0)
        determinant=s.factor(M.extract((0,5),(0,5)).det())
        self.assertEqual(s.simplify(determinant-(1-3*h*h)/4),0)
        self.assertNotEqual(s.simplify(determinant-(1-4*h*h)/4),0)
        self.assertEqual(len(collinear_blocks('gluon')['determinants']),2)

if __name__=='__main__':unittest.main()
