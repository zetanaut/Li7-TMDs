"""Direct lower-spin B tensor construction versus minimal K<=2 correlator.

Reference: Boer et al., arXiv:1607.01654v2, eqs. 11,24-25,36-38,95.
Approved manuscript: eq:referenceSTF,literal_gluon_traces,gproject,
  gluon_scalar_dictionary,gluon_linear_dictionary.
All transverse Lorentz contractions lower indices with g_T=-delta.
"""
from __future__ import annotations
from itertools import product
import sympy as S
from transverse_foundations import cartesian_tensor,E,stf2
from correlator_foundations import catalogue,cartesian_value

R=S.Rational
kx,ky,x=S.symbols('kx ky x',real=True)
M=S.symbols('M_A',positive=True)
SL,STx,STy,SLL,SLTx,SLTy,STTxx,STTxy=S.symbols('S_L S_Tx S_Ty S_LL S_LTx S_LTy S_TTxx S_TTxy',real=True)
TARGET={(0,0,0):S.Integer(1),(1,0,0):SL,(1,1,0):STx,(1,1,1):STy,
        (2,0,0):SLL,(2,1,0):SLTx,(2,1,1):SLTy,(2,2,0):STTxx,(2,2,1):STTxy}
B_NAMES=('f1','f1T_perp','f1LL','f1LT','f1TT','g1','g1T','g1LT','g1TT',
         'h1_perp','h1L_perp','h1','h1T_perp','h1LL_perp','h1LT','h1LT_perp',
         'h1TT','h1TT_perp','h1TT_perpperp')
B={name:S.Symbol('B_'+name,real=True) for name in B_NAMES}

# Independently reviewed dictionary fixture from the approved manuscript.
FIXTURE={
    ('f',0,0,0):B['f1'],('f',1,1,1):B['f1T_perp'],
    ('f',2,0,0):B['f1LL'],('f',2,1,1):B['f1LT']+B['h1LT'],
    ('f',2,2,2):B['f1TT']-B['h1TT_perp'],
    ('g',1,0,0):B['g1'],('g',1,1,1):B['g1T'],
    ('g',2,1,1):B['g1LT'],('g',2,2,2):-B['g1TT'],
    ('h',0,0,2):B['h1_perp'],('h',1,0,2):-B['h1L_perp'],
    ('h',1,1,1):-B['h1'],('h',1,1,3):-B['h1T_perp'],
    ('h',2,0,2):B['h1LL_perp'],('h',2,1,1):2*B['h1LT'],
    ('h',2,1,3):-B['h1LT_perp'],('h',2,2,0):B['h1TT'],
    ('h',2,2,4):B['h1TT_perpperp'],
}

def _ktensors():
    return {n:cartesian_tensor(n,kx,ky,M) for n in range(5)}

def _mat(t):return S.Matrix(2,2,lambda a,b:t[a,b])

def reference_matrix():
    """Literal reference free-index tensors, before any dictionary use."""
    k=_ktensors();v=S.Matrix([kx/M,ky/M]);K2=_mat(k[2]);K4=k[4]
    ST=S.Matrix([STx,STy]);LT=S.Matrix([SLTx,SLTy]);TT=S.Matrix([[STTxx,STTxy],[STTxy,-STTxx]])
    I=S.eye(2)
    Cst=(ST.T*E*v)[0]
    Clt=(LT.dot(v))
    Ctt=sum(K2[a,b]*TT[a,b] for a,b in product(range(2),repeat=2))
    Stt=-sum(E[b,g]*K2[g,a]*TT[a,b] for a,b,g in product(range(2),repeat=3))
    # Stt is the literal circular contraction: epsilon^b{}_g=-E[b,g].
    V3ST=S.Matrix(2,2,lambda a,b:sum(k[3][a,b,c]*ST[c] for c in range(2)))
    V3LT=S.Matrix(2,2,lambda a,b:sum(k[3][a,b,c]*LT[c] for c in range(2)))
    V4TT=S.Matrix(2,2,lambda a,b:sum(K4[a,b,c,d]*TT[c,d] for c,d in product(range(2),repeat=2)))
    # epsilon^i{}_alpha=-E[i,alpha], S_alpha=-S[alpha].
    Ldual=-E*K2
    Tlow=-(E*v*ST.T+ST*(E*v).T+E*ST*v.T+v*(E*ST).T)/4
    Thigh=-E*V3ST
    # B source LT symmetrization has no factor 1/2.
    LTlinear=LT*v.T+v*LT.T
    # TT mixed Lorentz index is lowered: S^i{}_alpha=-S^{i alpha}.
    TTmixed=-(TT*K2+(TT*K2).T)
    matrix=(I*B['f1']+K2*B['h1_perp']+
        S.I*E*SL*B['g1']+Ldual*SL*B['h1L_perp']+
        I*Cst*B['f1T_perp']+S.I*E*ST.dot(v)*B['g1T']+
        Tlow*B['h1']+Thigh*B['h1T_perp']+
        I*SLL*B['f1LL']+K2*SLL*B['h1LL_perp']+
        I*Clt*B['f1LT']+S.I*E*(LT.T*E*v)[0]*B['g1LT']+
        LTlinear*B['h1LT']-V3LT*B['h1LT_perp']+
        I*Ctt*B['f1TT']+S.I*E*Stt*B['g1TT']+
        TT*B['h1TT']+TTmixed*B['h1TT_perp']+V4TT*B['h1TT_perpperp'])
    return S.expand(x*matrix/2)


def minimal_matrix(mapping=FIXTURE,tensor_mass=M):
    """Independent minimal Cartesian correlator, with separate fixture."""
    result=S.zeros(2)
    for label in catalogue('gluon'):
        if label.K>2:continue
        key=(label.channel,label.K,label.m,label.n)
        coefficient=mapping[key]
        for part in range(1 if label.m==0 else 2):
            F,G,Hx,Hy=cartesian_value(label,part,kx,ky,tensor_mass)
            result += TARGET[label.K,label.m,part]*coefficient*(S.eye(2)*F+S.I*E*G+
                       S.Matrix([[Hx,Hy],[Hy,-Hx]]))
    return S.expand(x*result/2)


def projections(matrix):
    return (S.expand(S.trace(matrix)/x),
            S.expand(-S.I*sum(E[a,b]*matrix[a,b] for a,b in product(range(2),repeat=2))/x),
            S.expand(stf2(matrix)*2/x))


def check_difference(mapping=FIXTURE):
    ref=reference_matrix();minimal=minimal_matrix(mapping)
    return S.Matrix(2,2,lambda a,b:S.cancel(ref[a,b]-minimal[a,b]))


def derive_conversion_at_exact_point():
    """Project the direct reference tensors, then solve for minimal columns."""
    from correlator_foundations import map_matrix, coordinates
    reference=reference_matrix()
    ref_columns=[]
    substitutions={kx:S.Integer(2),ky:S.Integer(1),M:S.Integer(3)}
    target_variables=tuple(value for value in TARGET.values() if isinstance(value,S.Symbol))
    for name in B_NAMES:
        F,G,H=projections(S.diff(reference,B[name]))
        column=[]
        for K,m,part in coordinates():
            if K>2:continue
            target=TARGET[K,m,part]
            components=(F,G,H[0,0],H[0,1])
            for component in components:
                projected=component if target==1 else S.diff(component,target)
                column.append(S.cancel(projected.subs({v:0 for v in target_variables}).subs(substitutions)))
        ref_columns.append(S.Matrix(column))
    reference_map=S.Matrix.hstack(*ref_columns)
    labels,minimal=map_matrix('gluon',S.Integer(2),S.Integer(1),S.Integer(3))
    minimal_lower=minimal[:36,:18]
    rows=minimal_lower.T.rref()[1]
    if len(rows)!=18:raise AssertionError('gluon.dictionary_minimal_rank')
    conversion=minimal_lower.extract(rows,list(range(18))).inv()*reference_map.extract(rows,list(range(19)))
    if minimal_lower*conversion!=reference_map:
        raise AssertionError('gluon.dictionary_spanning')
    return conversion


def run_checks():
    diff=check_difference()
    if diff!=S.zeros(2):raise AssertionError(f'gluon.dictionary_complete: {diff}')
    ref=reference_matrix()
    # Rank of the reference-to-minimal coefficient map (19 -> 18).
    min_labels=[l for l in catalogue('gluon') if l.K<=2]
    conversion=S.Matrix([[S.diff(FIXTURE[l.channel,l.K,l.m,l.n],B[name]) for name in B_NAMES]
                         for l in min_labels])
    derived=derive_conversion_at_exact_point()
    if derived!=conversion:raise AssertionError('gluon.dictionary_derived_fixture_disagreement')
    if conversion.rank()!=18:raise AssertionError('gluon.dictionary_rank')
    kernel=conversion.nullspace()
    expected=S.Matrix([int(name in ('f1TT','h1TT_perp')) for name in B_NAMES])
    if len(kernel)!=1 or conversion*expected!=S.zeros(18,1):raise AssertionError('gluon.dictionary_kernel')
    if S.simplify(ref.subs({B['f1TT']:B['f1TT']+1,B['h1TT_perp']:B['h1TT_perp']+1})-ref)!=S.zeros(2):
        raise AssertionError('gluon.dictionary_reference_null')
    m0=S.symbols('M_0',positive=True)
    scaled={key:(m0/M)**key[3]*coefficient for key,coefficient in FIXTURE.items()}
    old=minimal_matrix()
    scaled_at_new=minimal_matrix(scaled,tensor_mass=m0)
    if any(S.cancel(scaled_at_new[a,b]-old[a,b])!=0 for a,b in product(range(2),repeat=2)):
        raise AssertionError('gluon.dictionary_mass_conversion')
    ll_scaled=dict(FIXTURE)
    for key in (('f',2,0,0),('h',2,0,2)):
        ll_scaled[key]=R(3,2)*ll_scaled[key]
    equal_cartesian=ref.subs(SLL,R(3,2)*SLL)
    if any(S.cancel(equal_cartesian[a,b]-minimal_matrix(ll_scaled)[a,b])!=0
           for a,b in product(range(2),repeat=2)):
        raise AssertionError('gluon.dictionary_cartesian_polarization')
    cgamma=S.symbols('C_Gamma',real=True)
    normalized=minimal_matrix({key:cgamma*value for key,value in FIXTURE.items()})
    if any(S.cancel((cgamma*ref-normalized)[a,b])!=0 for a,b in product(range(2),repeat=2)):
        raise AssertionError('gluon.dictionary_overall_normalization')
    return {
       'gluon.dictionary_complete':{'reference_coefficients':19,'minimal_coefficients':18,'matrix_entries':4,
            'domain':'Q[x,kx,ky,M_A^-1,target components]; x generic, M_A>0'},
       'gluon.dictionary_mass_polarization':{'mass':'symbolic M_A,M_0 positive','polarization':'named versus equal Cartesian S_LL','overall':'generic C_Gamma and x/2'},
       'gluon.dictionary_derived_conversion':{'reference_columns':19,'minimal_rows':36,'derived_entries':18*19,'point':['2','1','3'],'method':'independent direct-tensor projection'},
       'gluon.dictionary_rank_kernel':{'rank':18,'kernel':['B_f1TT + B_h1TT_perp'],
            'normalization':'common x/2 removed by gproject'},
    }
