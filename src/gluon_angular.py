"""Independent finite-mode and product-integral octupole Gram checks."""
from __future__ import annotations
import json
import hashlib
from pathlib import Path
import sympy as s

# id, beam parity, Fourier (varphi,psi), trigonometric type, analyzer phase.
MODES=[('g.30',1,0,0,'one',0),('h.30.2',0,2,0,'sin',-2)]
for m in (1,2,3):
    MODES.extend(((f'f.3{m}',0,m,-m,'sin',0),
                  (f'g.3{m}',1,m,-m,'cos',0),
                  (f'h.3{m}.{abs(m-2)}',0,2-m,m,'sin',-2),
                  (f'h.3{m}.{m+2}',0,m+2,-m,'sin',-2)))


def coefficients(mode):
    _,eps,p,q,kind,shift=mode
    if kind=='one':return {(eps,0,0):s.Integer(1)}
    c=s.exp(s.I*shift*s.Symbol('beta',real=True))
    if kind=='cos':return {(eps,p,q):c/2,(eps,-p,-q):s.conjugate(c)/2}
    return {(eps,p,q):c/(2*s.I),(eps,-p,-q):-s.conjugate(c)/(2*s.I)}


def fourier_gram():
    out=s.zeros(len(MODES))
    for i,a in enumerate(MODES):
        ac=coefficients(a)
        for j,b in enumerate(MODES):
            bc=coefficients(b)
            out[i,j]=s.simplify(sum(x*y for (ea,pa,qa),x in ac.items()
                for (eb,pb,qb),y in bc.items()
                if ea==eb and pa==-pb and qa==-qb))
    return out


def product_integral_gram():
    # Independent trigonometric orthogonality: beam parity is averaged;
    # each sine/cosine product is integrated on the two-torus.
    out=s.zeros(len(MODES))
    for i,(_,ea,pa,qa,ta,sa) in enumerate(MODES):
        for j,(_,eb,pb,qb,tb,sb) in enumerate(MODES):
            if ea!=eb:continue
            if ta=='one' or tb=='one':out[i,j]=int(ta==tb=='one');continue
            if (pa,qa)!=(pb,qb) and (pa,qa)!=(-pb,-qb):continue
            sign=1 if (pa,qa)==(pb,qb) else -1
            phase=(sa-sign*sb)*s.Symbol('beta',real=True)
            if ta==tb:out[i,j]=s.cos(phase)/2
            elif ta=='sin':out[i,j]=s.sin(phase)/2
            else:out[i,j]=-s.sin(phase)/2
    return out


def certificate():
    from validation_evidence import scientific_source_digest
    from gluon_response_fixture import OCT_MODES
    assert len(MODES)==14 and len({x[0] for x in MODES})==14
    if any(eps!=int(name.startswith('g.')) for name,eps,*_ in MODES):
        raise AssertionError('beam-helicity mode assignment changed')
    if tuple(MODES)!=OCT_MODES:raise AssertionError('octupole source-mode fixture mismatch')
    A=fourier_gram();B=product_integral_gram()
    if A!=B:raise AssertionError('independent angular integrals disagree')
    expected=s.diag(1,*([s.Rational(1,2)]*13))
    if A!=expected:raise AssertionError('octupole Gram witness changed')
    freqs=[(x[2],x[3]) for x in MODES]
    assert max(abs(p) for p,q in freqs)==5 and max(abs(q) for p,q in freqs)==3
    return {'source_label':'tab:oct_fourier','modes':[list(x) for x in MODES],
            'gram':[[str(A[i,j]) for j in range(14)] for i in range(14)],
            'determinant':str(A.det()),'independent_routes':['Fourier convolution','trigonometric product integral'],
            'no_alias_grid':[11,8],
            'convention':'conjugate-amplitude-first; physical electron h; 11x8 finite Fourier grid',
            'scientific_source_digest':scientific_source_digest(),
            'mode_fixture_sha256':hashlib.sha256((Path(__file__).parent/'gluon_response_fixture.py').read_bytes()).hexdigest()}


def save(path):
    Path(path).write_text(json.dumps(certificate(),sort_keys=True,indent=2)+'\n')


def check(path):
    actual=json.loads(Path(path).read_text())
    if actual!=certificate():raise AssertionError('angular certificate mismatch')
    return True
