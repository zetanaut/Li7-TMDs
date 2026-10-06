"""Source-indexed exact collinear block witness and separate checker."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import sympy as s
from conditional_positivity import collinear_blocks
from validation_evidence import scientific_source_digest

def certificate():
    return {'schema_version':1,'source_labels':['eq:jointM','eq:Oraise','eq:qSoffer','eq:gSoffer'],
            'convention':'parton index first; upper off-diagonal +i; factor 1/2 retained',
            'source_digest':scientific_source_digest(),
            'quark':collinear_blocks('quark'),'gluon':collinear_blocks('gluon')}

def save(path):
    Path(path).write_text(json.dumps(certificate(),indent=2,sort_keys=True)+'\n')

def check(path):
    """Rebuild witness with explicit matrix entries, separate from its graph finder."""
    data=json.loads(Path(path).read_text())
    if data.get('schema_version')!=1 or data.get('source_digest')!=scientific_source_digest():
        raise ValueError('stale collinear block certificate')
    f,fq,g,go,h,ho=s.symbols('f f_Q g g_O h h_O',real=True)
    Uo,Ui=f+fq,f-fq
    Go,Gi=s.Rational(3,2)*g+s.Rational(3,10)*go,g/2-s.Rational(9,10)*go
    a=s.sqrt(3)*(h+s.Rational(2,5)*ho)
    b=2*h-s.Rational(6,5)*ho
    expected={
      'quark':[(Uo+Go,Ui-Gi,a),(Ui+Gi,Ui+Gi,b),(Uo+Go,Ui-Gi,a)],
      'gluon':[(Uo+Go,Ui+Gi,2*s.sqrt(3)*h)]*2,
    }
    expected_blocks={
      'quark':[s.Matrix([[Uo+Go,a],[a,Ui-Gi]])/2,
               s.Matrix([[Ui+Gi,b],[b,Ui+Gi]])/2,
               s.Matrix([[Ui-Gi,a],[a,Uo+Go]])/2,
               s.Matrix([[(Uo-Go)/2]]),s.Matrix([[(Uo-Go)/2]])],
      'gluon':[s.Matrix([[Uo+Go,2*s.sqrt(3)*h],[2*s.sqrt(3)*h,Ui+Gi]])/2,
               s.Matrix([[Ui+Gi,2*s.sqrt(3)*h],[2*s.sqrt(3)*h,Uo+Go]])/2,
               s.Matrix([[(Ui-Gi)/2]]),s.Matrix([[(Uo-Go)/2]]),
               s.Matrix([[(Uo-Go)/2]]),s.Matrix([[(Ui-Gi)/2]])],
    }
    wanted_groups={'quark':[[0,5],[1,6],[2,7],[3],[4]],
                   'gluon':[[0,6],[1,7],[2],[3],[4],[5]]}
    for species in ('quark','gluon'):
        row=data[species]
        groups=row['groups'];order=row['permutation']
        if sorted(order)!=list(range(8)) or sum(map(len,groups))!=8:
            raise ValueError('invalid simultaneous permutation')
        blocks=[s.Matrix([[s.sympify(v,locals={'f':f,'f_Q':fq,'g':g,'g_O':go,'h':h,'h_O':ho})
                          for v in line] for line in block]) for block in row['blocks']]
        circular=s.Matrix([[1,1],[-s.I,s.I]])/s.sqrt(2)
        expected_basis=s.eye(8) if species=='quark' else s.kronecker_product(circular,s.eye(4))
        encoded_basis=s.Matrix([[s.sympify(v) for v in line] for line in row['basis']])
        if encoded_basis!=expected_basis or groups!=wanted_groups[species] or order!=sum(groups,[]):
            raise ValueError('source basis or block permutation mismatch')
        if len(blocks)!=len(expected_blocks[species]) or any(
             s.simplify(got-want)!=s.zeros(got.rows) for got,want in zip(blocks,expected_blocks[species])):
            raise ValueError('source amplitude or polarization block entry mismatch')
        pair=[B for B in blocks if B.rows==2]
        if len(pair)!=len(expected[species]) or any(B!=B.H for B in blocks):
            raise ValueError('block dimensions or Hermiticity')
        wanted=sorted([s.factor((x*y-z*z)/4) for x,y,z in expected[species]],key=str)
        got=sorted([s.factor(B.det()) for B in pair],key=str)
        if wanted!=got or list(map(str,got))!=sorted(row['determinants']):
            raise ValueError('independent source bound determinant mismatch')
        if any(B[0,0]==0 and B[1,1]==0 and B[0,1]!=0 for B in pair):
            raise ValueError('impossible block')
    return {'certificate_sha256':hashlib.sha256(Path(path).read_bytes()).hexdigest(),
            'quark_blocks':len(expected['quark']),'gluon_blocks':len(expected['gluon']),
            'method':'explicit source amplitudes and determinants versus stored block graph'}
