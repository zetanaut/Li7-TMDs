"""Wrong-order and wrong-interpretation controls for checkpoint 6."""
import sys
from pathlib import Path
import unittest
import numpy as np
import sympy as s

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from born_direct import G,SIGMA,massless,direct
from current_index_conventions import (G5,exact_current_identity,spinor_anchor,annihilation_spinor_anchor,ordered_currents,
                                       sidis_ordering,born_case)
from gluon_born_response import leptonic,epsilon_lower
from source_joint_positivity import (source_mapping_check,partial_parton_transpose,
                                     spectral_index_check,spectral_recovery)
from validation_evidence import CONVENTION_REVIEW,assert_convention_review,EvidenceError


class ConventionClosureTests(unittest.TestCase):
    def test_spinor_eigenvalue_not_adapter(self):
        p=np.array([5.,1.5,2.,np.sqrt(25-1.5**2-2.**2)])
        for h in (-1,1):
            anchor=spinor_anchor(p,h)
            self.assertLess(anchor['helicity_eigen_residual'],1e-13)
            self.assertLess(anchor['gamma5_residual'],1e-13)
            self.assertLess(anchor['projector_residual'],1e-13)
            self.assertAlmostEqual(anchor['normalization'],10.)
            self.assertGreater(spinor_anchor(p,-h)['helicity_eigen_residual']+
                               np.linalg.norm(massless(p,h)-massless(p,-h)),1.)
        rows=annihilation_spinor_anchor()
        self.assertEqual(len(rows),4)
        self.assertTrue(all(row['spinor_helicity_eigenvalue']==
                            (-row['input_label'] if row['branch']=='antiparticle_v'
                             else row['input_label']) for row in rows))

    def test_current_order_and_epsilon_lowering(self):
        self.assertEqual(exact_current_identity()['components'],16)
        l=np.array([5.,0.,0.,5.]);lp=np.array([5.,3.,0.,4.])
        x=ordered_currents(l,lp,1)
        self.assertAlmostEqual(x['J'][1,2].imag,-10.)
        self.assertAlmostEqual(x['reverse'][1,2].imag,10.)
        self.assertAlmostEqual(x['source'][1,2].imag,10.)
        self.assertLess(x['current_to_source_transpose'],1e-12)
        self.assertGreater(x['current_to_literal_source'],19.)
        self.assertEqual(epsilon_lower(0,1,2,3),-1)
        self.assertEqual(epsilon_lower(1,0,2,3),1)
        self.assertAlmostEqual(np.trace(G[0]@G[1]@G[2]@G[3]@G5).imag,-4.)

    def test_literal_mismatch_candidate_is_separate(self):
        case=born_case(dict(sqrt_s=5.,Q2=4.,mass=1.5,theta=.8,
                            phi=.4,lepton_energy=10.,helicity=1))
        self.assertGreater(case['literal_max_abs'],200.)
        self.assertLess(case['candidate_max_abs'],1e-8)
        self.assertGreater(abs(case['circular_rate_physical']-
                               case['circular_rate_literal']),100.)
        self.assertGreater(abs(case['elliptic_rate_physical']-
                               case['elliptic_rate_literal']),1.)
        self.assertLess(case['unchanged_UCS_max_abs'],1e-8)

    def test_sidis_order_is_distinct(self):
        x=sidis_ordering()
        self.assertLess(x['trace_to_reverse_hard'],1e-12)
        self.assertGreater(x['trace_to_amplitude_first'],1e-2)
        self.assertGreater(x['complex_offdiagonal'],1e-2)

    def test_source_auxiliary_mapping_and_rank(self):
        x=source_mapping_check()
        self.assertEqual(x['source_indexed_columns'],{'quark':32,'gluon':32})
        self.assertTrue(x['coefficient_rank_preserved'])
        self.assertFalse(x['matrix_rank_or_psd_preserved'])
        self.assertEqual(len(x['minimal_entries']),2)
        self.assertEqual(spectral_index_check()['columns_checked'],{'quark':32,'gluon':32})

    def test_partial_transpose_not_psd(self):
        bell=s.zeros(8,1);bell[0]=bell[5]=1/s.sqrt(2)
        gram=bell*bell.H
        self.assertEqual(gram.rank(),1)
        pt=partial_parton_transpose(gram)
        self.assertIn(-s.Rational(1,2),pt.eigenvals())
        self.assertEqual(pt.rank(),4)

    def test_finite_k_partial_transpose_negative(self):
        for species in ('quark','gluon'):
            x=spectral_recovery(species)
            self.assertEqual(x['source_map_rank'],32)
            self.assertFalse(x['partial_transpose_applied_to_Gram'])
            self.assertLess(x['finite_k_partial_transpose_min_eigenvalue'],-.1)

    def test_review_status_cannot_be_suppressed(self):
        self.assertEqual(CONVENTION_REVIEW['physical_helicity_anchor_status'],
                         'VERIFIED_EIGENVALUE')
        self.assertEqual(CONVENTION_REVIEW['source_formula_agreement_status'],
                         'LITERAL_PHYSICAL_HELICITY_MISMATCH')
        self.assertEqual(CONVENTION_REVIEW['publication_eligibility'],
                         'BLOCKED_AUTHOR_REVIEW')
        self.assertTrue(CONVENTION_REVIEW['unresolved_scientific_issues'])
        manifest={'convention_review':dict(CONVENTION_REVIEW)}
        payload={'status':{'diagnostic_execution':'PASS',
                           'literal_physical_source_agreement':'MISMATCH',
                           'candidate_index_interchange':'AGREES',
                           'publication_eligibility':'BLOCKED_AUTHOR_REVIEW'}}
        assert_convention_review(manifest,payload)
        bad={'convention_review':dict(CONVENTION_REVIEW,
                                      publication_eligibility='ELIGIBLE')}
        with self.assertRaises(EvidenceError):assert_convention_review(bad,payload)
        bad_payload={'status':dict(payload['status'],
                                   literal_physical_source_agreement='AGREES')}
        with self.assertRaises(EvidenceError):assert_convention_review(manifest,bad_payload)


if __name__=='__main__':unittest.main()
