"""Exact, independently checkable minor witnesses for the 64x32 maps."""
from __future__ import annotations
import json
from dataclasses import asdict
import re
from pathlib import Path
import sympy as S
from sympy.polys.matrices import DomainMatrix
from correlator_foundations import map_matrix, row_definitions
from validation_evidence import (REFERENCE_SOURCE_SHA256, scientific_source_digest,
                                 file_digest, digest)

ROOT=Path(__file__).resolve().parents[1]
CERTIFICATE_VERSION=2
POINT={'kx':'2','ky':'1','M_A':'3'}
RATIONAL=re.compile(r'-?(?:0|[1-9][0-9]*)(?:/[1-9][0-9]*)?\Z')


def rational(value):
    if not isinstance(value,str) or not RATIONAL.fullmatch(value):
        raise ValueError('certificate entry is not a canonical rational')
    result=S.Rational(value)
    if str(result)!=value:raise ValueError('certificate entry is not reduced')
    return result


def matrix_digest(matrix):
    return digest([[str(matrix[i,j]) for j in range(matrix.cols)] for i in range(matrix.rows)])


def build_certificate(species):
    if species not in ('quark','gluon'):raise ValueError('unknown species')
    labels,M=map_matrix(species,route='cartesian')
    selected=list(M.T.rref()[1])
    if len(selected)!=32:raise ValueError('map lacks 32 independent rows')
    minor=M.extract(selected,list(range(32)))
    determinant=DomainMatrix.from_Matrix(minor).det()
    if determinant==0:raise ValueError('selected minor vanishes')
    return {
        'format_version':CERTIFICATE_VERSION,'species':species,
        'reference_sha256':REFERENCE_SOURCE_SHA256,
        'scientific_source_digest':scientific_source_digest(),
        'exact_field':'Q','point':POINT,
        'columns':[label.id() for label in labels],
        'column_semantics':[asdict(label) for label in labels],
        'rows':[[K,m,part,channel] for K,m,part,channel in row_definitions()],
        'selected_rows':selected,
        'minor':[[str(minor[i,j]) for j in range(32)] for i in range(32)],
        'determinant':str(determinant),
        'method':'transpose-RREF row selection; exact DomainMatrix determinant',
    }


def verify_certificate(value,*,require_current=True):
    if not isinstance(value,dict):raise ValueError('certificate is not an object')
    if value.get('format_version')!=CERTIFICATE_VERSION or value.get('species') not in ('quark','gluon'):
        raise ValueError('certificate version or species mismatch')
    if value.get('reference_sha256')!=REFERENCE_SOURCE_SHA256:raise ValueError('reference digest mismatch')
    if require_current and value.get('scientific_source_digest')!=scientific_source_digest():
        raise ValueError('source digest mismatch')
    if value.get('exact_field')!='Q' or value.get('point')!=POINT:
        raise ValueError('exact field or evaluation point mismatch')
    species=value['species']
    labels,M=map_matrix(species,route='helicity')
    if (value.get('columns')!=[label.id() for label in labels] or
        value.get('column_semantics')!=[asdict(label) for label in labels]):
        raise ValueError('semantic column order mismatch')
    if value.get('rows')!=[[K,m,part,channel] for K,m,part,channel in row_definitions()]:
        raise ValueError('row convention mismatch')
    selected=value.get('selected_rows')
    if (not isinstance(selected,list) or len(selected)!=32 or
        any(type(x)!=int or x<0 or x>=64 for x in selected) or len(set(selected))!=32):
        raise ValueError('invalid independent-row selection')
    raw=value.get('minor')
    if not isinstance(raw,list) or len(raw)!=32 or any(not isinstance(row,list) or len(row)!=32 for row in raw):
        raise ValueError('malformed minor')
    minor=S.Matrix([[rational(entry) for entry in row] for row in raw])
    expected=M.extract(selected,list(range(32)))
    if minor!=expected:raise ValueError('minor entry disagrees with independent helicity construction')
    determinant=DomainMatrix.from_Matrix(minor).det()
    if determinant==0 or str(determinant)!=value.get('determinant'):
        raise ValueError('determinant mismatch or zero')
    return {'species':species,'rank_lower_bound':len(selected),'matrix_shape':[64,32],
            'determinant':str(determinant),'selected_rows':selected,
            'certificate_sha256':digest(value),'method':'recomputed helicity entries and exact minor'}


def read_certificate(path):
    return json.loads(Path(path).read_text(encoding='utf-8'),
                      object_pairs_hook=_unique_pairs)


def _unique_pairs(pairs):
    out={}
    for key,v in pairs:
        if key in out:raise ValueError('duplicate certificate key')
        out[key]=v
    return out


def save_certificates(directory):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    for species in ('quark','gluon'):
        value=build_certificate(species)
        (directory/f'{species}_rank.json').write_text(json.dumps(value,indent=2,sort_keys=True)+'\n',encoding='utf-8')
