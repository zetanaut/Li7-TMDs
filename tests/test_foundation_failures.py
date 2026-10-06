"""Injected scientific mistakes must fail their named mathematical diagnostic."""
from __future__ import annotations
import copy
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import sympy as S
from spin_foundations import cartesian,reconstruct,zero,I4,R
from spin_state_foundations import x3
from transverse_foundations import cartesian_tensor,dyadic_regression
from correlator_foundations import catalogue,map_matrix,cartesian_value,helicity_value
from foundation_certificates import build_certificate,verify_certificate
from gluon_dictionary_foundations import B,FIXTURE,check_difference,reference_matrix,minimal_matrix,SLL
from projector_foundations import analytic_recover

class FoundationFailureTests(unittest.TestCase):
    def test_subtraction_coefficients(self):
        J,Q,O,cq,co=cartesian()
        self.assertFalse(zero(sum((Q[a,a]+I4 for a in range(3)),S.zeros(4))))
        bad=sum((O[a,a,2]+R(1,10)*(J[2]+2*int(a==2)*J[a]) for a in range(3)),S.zeros(4))
        self.assertFalse(zero(bad))

    def test_imaginary_density_and_multiplicity(self):
        J,Q,O,_,_=cartesian()
        Xmat=S.zeros(4);Xmat[0,1]=S.I;Xmat[1,0]=-S.I
        real_only=S.re(reconstruct(Xmat,J,Q,O))
        self.assertFalse(zero(real_only-Xmat))
        Ymat=Q[0,1]
        wrong=S.trace(Ymat)*I4/4+sum((S.trace(Ymat*z)*z for z in J),S.zeros(4))/5
        wrong+=sum((S.trace(Ymat*Q[a,b])*Q[a,b] for a in range(3) for b in range(a,3)),S.zeros(4))/6
        wrong+=R(2,9)*sum((S.trace(Ymat*z)*z for z in O.values()),S.zeros(4))
        self.assertFalse(zero(wrong-Ymat))

    def test_reversed_rotation_phase(self):
        J,_,_,_,_=cartesian()
        phi=S.pi/2;theta=2*S.atan(R(1,2))
        from sympy.physics.wigner import wigner_d_small
        from spin_foundations import HELICITIES
        uy=wigner_d_small(R(3,2),-theta)
        good=S.diag(*(S.exp(-S.I*phi*m) for m in HELICITIES))*uy
        bad=S.diag(*(S.exp(S.I*phi*m) for m in HELICITIES))*uy
        target=R(4,5)*J[1]+R(3,5)*J[2]
        self.assertTrue(zero(good*J[2]*good.H-target))
        self.assertFalse(zero(bad*J[2]*bad.H-target))

    def test_stf_dyadic_and_mass(self):
        target,pure,dyadic,cstf,dstf=dyadic_regression(S.Integer(2),S.Integer(1))
        self.assertEqual(cstf,S.zeros(2))
        self.assertNotEqual(dstf,S.zeros(2))
        normalized=cartesian_tensor(3,S.Integer(2),S.Integer(1),S.Integer(3))
        missing=cartesian_tensor(3,S.Integer(2),S.Integer(1),S.Integer(1))
        self.assertNotEqual(normalized,missing)

    def test_parity_dual_and_epsilon(self):
        label=next(z for z in catalogue('quark') if z.channel=='h' and z.K==0)
        good=cartesian_value(label,0,S.Integer(2),S.Integer(1),S.Integer(3))
        raw=good[3]-S.I*good[2]  # inverse of clockwise E action
        self.assertNotEqual(good[2]+S.I*good[3],raw)
        self.assertNotEqual(good[3],-good[3])

    def test_duplicate_omit_and_rank_preserving_sign(self):
        _,matrix=map_matrix('quark')
        self.assertLess(matrix[:,1:].rank(),32)
        duplicate=matrix.copy();duplicate[:,1]=duplicate[:,0]
        self.assertLess(duplicate.rank(),32)
        flip=matrix.copy();flip[:,1]=-flip[:,1]
        self.assertEqual(flip.rank(),32)
        self.assertNotEqual(flip,matrix)

    def test_lt_shift_and_polarization_mapping(self):
        wrong=dict(FIXTURE);wrong[('f',2,1,1)]=B['f1LT']
        self.assertNotEqual(check_difference(wrong),S.zeros(2))
        wrong_lt=dict(FIXTURE);wrong_lt[('h',2,1,1)]=B['h1LT']
        self.assertNotEqual(check_difference(wrong_lt),S.zeros(2))
        wrong_circle=dict(FIXTURE);wrong_circle[('g',2,2,2)]=B['g1TT']
        self.assertNotEqual(check_difference(wrong_circle),S.zeros(2))
        ref=reference_matrix();named=minimal_matrix()
        self.assertEqual(ref,named)
        cartesian=ref.subs(SLL,R(3,2)*SLL)
        self.assertNotEqual(cartesian,ref)

    def test_certificate_entry_label_digest_and_sign(self):
        cert=build_certificate('quark')
        verify_certificate(cert)
        mutations=[]
        entry=copy.deepcopy(cert);entry['minor'][0][0]='2';mutations.append(entry)
        label=copy.deepcopy(cert);label['columns'][0],label['columns'][1]=label['columns'][1],label['columns'][0];mutations.append(label)
        digest=copy.deepcopy(cert);digest['scientific_source_digest']='0'*64;mutations.append(digest)
        sign=copy.deepcopy(cert)
        sign['minor']=[[str(-S.Rational(row[0]))]+row[1:] for row in sign['minor']]
        sign['determinant']=str(-S.Rational(sign['determinant']))
        self.assertNotEqual(S.Rational(sign['determinant']),0)
        mutations.append(sign)
        duplicate=copy.deepcopy(cert);duplicate['columns'][1]=duplicate['columns'][0];mutations.append(duplicate)
        point=copy.deepcopy(cert);point['point']['M_A']='1';mutations.append(point)
        for mutated in mutations:
            with self.assertRaisesRegex(ValueError,'minor entry|column order|source digest|evaluation point'):
                verify_certificate(mutated)

    def test_projector_singular_domain(self):
        with self.assertRaisesRegex(ValueError,'undefined at k_T=0'):
            analytic_recover('quark',[S.Integer(1)]*32,radial=0)
        with self.assertRaisesRegex(ValueError,'proven nonzero'):
            analytic_recover('gluon',[S.Integer(1)]*32,radial=S.symbols('t'))

    def test_transverse_octupole_has_llt(self):
        J,_,O,_,_=cartesian()
        rho=I4/4+S.sqrt(5)/12*x3((1,0,0),J)
        self.assertNotEqual(S.trace(rho*O[2,2,0]),0)

if __name__=='__main__':unittest.main()
