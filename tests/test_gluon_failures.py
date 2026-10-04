"""Test-local scientific mutations with named diagnostics."""
import copy
import sys
from pathlib import Path
import unittest
from unittest.mock import patch
import numpy as np
import sympy as s

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from gluon_stokes import density,reconstruct
from gluon_responses import verify_rows
import gluon_responses
from gluon_angular import certificate,MODES
import gluon_angular
from born_direct import _diagrams,fourvectors,massless,G,direct,compare,compare_source_label
from born_validation import dense_scan,verify_scan_payload,normalization_checks,raw_matrix_diagnostic
from gluon_reconstruction import run_reconstruction
from gluon_octupole import projection_checks
from gluon_born_response import evaluate


class GluonFailureTests(unittest.TestCase):
    def test_transposed_gram_and_helicity_sign(self):
        B=reconstruct(5,2,1,s.Rational(1,2))
        D=density(1,s.Rational(1,2),0,0)
        self.assertNotEqual(s.trace(B*D),s.trace(B.T*D))
        self.assertNotEqual(s.trace(B*D),s.trace(B*D.T))
        self.assertNotEqual(s.trace(B*D),5-2*s.Rational(1,2))
        e=np.array([1,1j])/np.sqrt(2)
        A=np.array([[1+2j,.3-.4j],[.7+.2j,-.6+1.1j]])
        gram=A@A.conj().T
        self.assertGreater(abs(np.trace(gram@np.outer(e.conj(),e))-
                               np.trace(gram@np.outer(e,e.conj()))),.1)
        with self.assertRaisesRegex(AssertionError,'direct/trace mismatch'):
            compare()
        self.assertLess(compare_source_label()['max_abs'],1e-9)

    def test_low_high_coefficients_and_odd_dual(self):
        rows=list(gluon_responses.ROWS)
        index=next(i for i,x in enumerate(rows) if x[:4]==(3,2,'h',0))
        bad=list(rows[index]);bad[4]=s.Rational(1,2);rows[index]=tuple(bad)
        with patch.object(gluon_responses,'ROWS',tuple(rows)),self.assertRaisesRegex(AssertionError,'gluon.h.32'):
            verify_rows()
        rows=list(gluon_responses.ROWS)
        index=next(i for i,x in enumerate(rows) if x[:4]==(3,3,'h',5))
        bad=list(rows[index]);bad[4]=s.Rational(1,2);rows[index]=tuple(bad)
        with patch.object(gluon_responses,'ROWS',tuple(rows)),self.assertRaisesRegex(AssertionError,'gluon.h.33'):
            verify_rows()
        rows=list(gluon_responses.ROWS)
        index=next(i for i,x in enumerate(rows) if x[:4]==(3,0,'h',2))
        bad=list(rows[index]);bad[7]='cos';rows[index]=tuple(bad)
        with patch.object(gluon_responses,'ROWS',tuple(rows)),self.assertRaisesRegex(AssertionError,'gluon.h.30'):
            verify_rows()

    def test_angular_row_and_beam_mutations(self):
        bad=list(MODES);bad[0]=('g.30',0,0,0,'one',0)
        with patch.object(gluon_angular,'MODES',bad),self.assertRaisesRegex(AssertionError,'beam-helicity'):
            certificate()
        bad=list(MODES);bad[-1]=bad[-2]
        with patch.object(gluon_angular,'MODES',bad),self.assertRaises(AssertionError):
            certificate()

    def test_born_diagram_and_normalization_mutations(self):
        ps=fourvectors();l,lp=(ps[k] for k in ('l','lp'))
        a=massless(l,-1);b=massless(lp,-1)
        current=np.array([np.vdot(b,G[0]@G[i]@a) for i in range(4)])
        amps=np.array([_diagrams(ps,1.5,current,np.eye(4)[i]).reshape(2,4) for i in (1,2)])
        correct=(amps.sum(axis=1))@amps.sum(axis=1).conj().T
        dropped=amps[:,0]@amps[:,0].conj().T
        flipped=(amps[:,0]-amps[:,1])@(amps[:,0]-amps[:,1]).conj().T
        incoherent=dropped+amps[:,1]@amps[:,1].conj().T
        for wrong in (dropped,flipped,incoherent,correct/2):
            self.assertGreater(np.max(np.abs(wrong-correct)),1.)
        photon=_diagrams(ps,1.5,ps['q'],np.eye(4)[2])
        gluon=_diagrams(ps,1.5,current,ps['k'])
        self.assertLess(np.max(np.abs(photon.sum(axis=0))),1e-12)
        self.assertLess(np.max(np.abs(gluon.sum(axis=0))),1e-11)
        self.assertGreater(np.max(np.abs(photon[0])),1.)
        self.assertGreater(np.max(np.abs(gluon[0])),1.)
        self.assertGreaterEqual(np.linalg.eigvalsh(incoherent).min(),-1e-10)
        self.assertEqual(normalization_checks()['T_F'],.5)
        self.assertNotEqual(normalization_checks()['T_F']/2,.5)
        corrupted=correct.copy();corrupted[1,0]+=3j
        with self.assertRaisesRegex(ValueError,'raw hard matrix'):
            raw_matrix_diagnostic(corrupted,float(np.max(np.abs(correct))))
        corrupted=correct.copy();corrupted[0,0]+=2j
        with self.assertRaisesRegex(ValueError,'raw hard matrix'):
            raw_matrix_diagnostic(corrupted,float(np.max(np.abs(correct))))

    def test_invalid_domain_and_scan_identity(self):
        with self.assertRaisesRegex(ValueError,'interior lepton energy'):
            direct(lepton_energy=2.5)
        payload=dense_scan()
        payload['cases'][17]=copy.deepcopy(payload['cases'][16])
        with self.assertRaisesRegex(ValueError,'input tuple mismatch'):
            verify_scan_payload(payload)

    def test_moment_prep_and_rank_loss(self):
        projected=projection_checks()
        self.assertEqual(projected['transverse_alpha'],['0','-1/2','0','5/2'])
        self.assertNotEqual(projected['transverse_alpha'][1],'0')
        result=run_reconstruction()
        self.assertEqual(result['rank'],14)
        self.assertLess(result['rank_losses']['aliased_4x4'],14)
        self.assertEqual(result['rank_losses']['O_zero'],0)
        row=next(x for x in result['normalized_moments'] if x['id']=='h.30.2')
        self.assertEqual(row['Z'],2)
        self.assertGreater(abs(row['measured']/2-row['expected_ratio']),1e-6)
        self.assertEqual(result['rank_losses']['bG_zero'],10)
        self.assertEqual(result['rank_losses']['bLin_zero'],7)


if __name__=='__main__':unittest.main()
