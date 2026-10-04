"""Reviewed Cartesian covariants and independent helicity/parity construction.

Rows are (target K,m,real/imag, parton F,G,H_re,H_im), in that order.
Target coordinates are formal response derivatives, not individual states.
Momentum components are physical; every rank-n tensor carries M_A^-n.
"""
from __future__ import annotations
from dataclasses import dataclass, asdict
from itertools import product
import sympy as S
from transverse_foundations import cartesian_tensor, E, stf2

R=S.Rational

@dataclass(frozen=True)
class Label:
    species:str
    channel:str
    K:int
    m:int
    n:int
    parent:str='A'
    flavor:str='unspecified'
    antiquark:bool|None=None
    link:str='unspecified'
    color:str='unspecified'
    def id(self):return f'{self.species}.{self.channel}.{self.K}{self.m}[{self.n}]'

# Route A: reviewed, explicit slot lists transcribed from Ffull/Gfull,
# T0--T3 and H0--B3. These are not produced by the helicity rule.
SCALAR_F=((0,0),(1,1),(2,0),(2,1),(2,2),(3,1),(3,2),(3,3))
SCALAR_G=((1,0),(1,1),(2,1),(2,2),(3,0),(3,1),(3,2),(3,3))
QUARK_H=(((0,0,1),),((1,0,1),(1,1,0),(1,1,2)),
         ((2,0,1),(2,1,0),(2,1,2),(2,2,1),(2,2,3)),
         ((3,0,1),(3,1,0),(3,1,2),(3,2,1),(3,2,3),(3,3,2),(3,3,4)))
GLUON_H=(((0,0,2),),((1,0,2),(1,1,1),(1,1,3)),
          ((2,0,2),(2,1,1),(2,1,3),(2,2,0),(2,2,4)),
          ((3,0,2),(3,1,1),(3,1,3),(3,2,0),(3,2,4),(3,3,1),(3,3,5)))

def catalogue(species):
    h=QUARK_H if species=='quark' else GLUON_H
    labels=[Label(species,'f',K,m,m) for K,m in SCALAR_F]
    labels += [Label(species,'g',K,m,m) for K,m in SCALAR_G]
    labels += [Label(species,'h',K,m,n) for group in h for K,m,n in group]
    return tuple(sorted(labels,key=lambda z:(z.K,z.m,{'f':0,'g':1,'h':2}[z.channel],z.n)))

def helicity_catalogue(species):
    """Route B: enumerate orbital helicities and reflection signs."""
    spin_flip=1 if species=='quark' else 2
    out=[]
    for K in range(4):
        for m in range(K+1):
            for channel in ('f','g'):
                if m or ((K%2==0)==(channel=='f')):
                    out.append(Label(species,channel,K,m,m))
            ranks=[spin_flip] if m==0 else [abs(m-spin_flip),m+spin_flip]
            for n in ranks:out.append(Label(species,'h',K,m,n))
    return tuple(out)

def target_tensor(m,part):
    if m==0:return {():S.Integer(1)}
    return {idx:S.Integer((-1)**(sum(idx)//2))*(int(part==0) if sum(idx)%2==0 else int(part==1))
            for idx in product(range(2),repeat=m)}

def _contract(A,B):return sum(A[idx]*B[idx] for idx in A)

def _dual(A):return {(a,)+idx:sum(E[a,b]*A[(b,)+idx] for b in range(2))
                       for a in range(2) for idx in product(range(2),repeat=len(next(iter(A)))-1)}

def cartesian_value(label,part,kx,ky,mass):
    """Route A explicit Cartesian tensor contractions from displayed terms."""
    K,m,n=label.K,label.m,label.n
    A=target_tensor(m,part)
    k={j:cartesian_tensor(j,kx,ky,mass) for j in range(6)}
    if label.channel in ('f','g'):
        value=_contract(A,k[m])
        ordinary=(K%2==0)==(label.channel=='f')
        if not ordinary:value=-_contract(_dual(A),k[m])
        return (value,0,0,0) if label.channel=='f' else (0,value,0,0)
    if label.species=='quark':
        if m==0:v=S.Matrix([k[1][(a,)] for a in range(2)])
        elif n==m-1:
            v=S.Matrix([sum(A[(a,)+tail]*k[m-1][tail] for tail in product(range(2),repeat=m-1)) for a in range(2)])
        else:
            v=S.Matrix([sum(k[m+1][(a,)+tail]*A[tail] for tail in product(range(2),repeat=m)) for a in range(2)])
        if K%2==0:v=E*v
        return (0,0,S.expand(v[0]),S.expand(v[1]))
    if m==0:B=S.Matrix(2,2,lambda a,b:k[2][a,b])
    elif n==abs(m-2):
        if m==1:
            a=S.Matrix([A[(j,)] for j in range(2)])
            v=S.Matrix([k[1][(j,)] for j in range(2)])
            B=stf2(a*v.T)
        elif m==2:B=S.Matrix(2,2,lambda a,b:A[a,b])
        else:B=S.Matrix(2,2,lambda a,b:sum(A[(a,b,c)]*k[1][(c,)] for c in range(2)))
    else:B=S.Matrix(2,2,lambda a,b:sum(k[m+2][(a,b)+tail]*A[tail]
                      for tail in product(range(2),repeat=m)))
    if K%2==1:B=E*B
    return (0,0,S.expand(B[0,0]),S.expand(B[0,1]))

def helicity_value(label,part,kx,ky,mass):
    """Route B uses complex helicity branch formulae and reflection parity."""
    K,m,n=label.K,label.m,label.n
    z=(kx+S.I*ky)/mass
    A=S.Integer(1) if part==0 else S.I
    if label.channel in ('f','g'):
        wave=S.expand(S.conjugate(A)*z**m)
        ordinary=(K%2==0)==(label.channel=='f')
        value=S.re(wave) if ordinary else S.im(wave)
        return (value,0,0,0) if label.channel=='f' else (0,value,0,0)
    if label.species=='quark':
        if m==0:v=z
        elif n==m-1:v=A*S.conjugate(z)**(m-1)
        else:v=S.conjugate(A)*z**(m+1)/2
        if K%2==0:v=-S.I*v
    else:
        if m==0:v=z*z/2
        elif n==abs(m-2):v={1:A*z/2,2:A,3:A*S.conjugate(z)}[m]
        else:v=S.conjugate(A)*z**(m+2)/4
        if K%2==1:v=-S.I*v
    return (0,0,S.expand(S.re(S.expand_complex(v))),S.expand(S.im(S.expand_complex(v))))

def coordinates():
    return tuple((K,m,part) for K in range(4) for m in range(K+1)
                 for part in range(1 if m==0 else 2))

def row_definitions():
    return tuple((K,m,part,channel) for K,m,part in coordinates()
                 for channel in ('F','G','H_re','H_im'))

def map_matrix(species,kx=S.Integer(2),ky=S.Integer(1),mass=S.Integer(3),route='cartesian'):
    labels=catalogue(species) if route=='cartesian' else helicity_catalogue(species)
    value=cartesian_value if route=='cartesian' else helicity_value
    columns=[]
    for label in labels:
        col=[]
        for K,m,part in coordinates():
            col.extend(value(label,part,kx,ky,mass) if (K,m)==(label.K,label.m) else (0,0,0,0))
        columns.append(S.Matrix(col))
    return labels,S.Matrix.hstack(*columns)

def run_checks():
    out={}
    for species in ('quark','gluon'):
        ca,A=map_matrix(species,route='cartesian')
        cb,B=map_matrix(species,route='helicity')
        if ca!=cb:raise AssertionError(f'{species}.semantic_catalogue')
        out[f'{species}.semantic_catalogue']={'labels':[asdict(x) for x in ca]}
        if A!=B:
            for row,col in product(range(64),range(32)):
                if A[row,col]!=B[row,col]:raise AssertionError(f'{species}.helicity_comparison row={row} col={col}')
        out[f'{species}.helicity_comparison']={'entries':A.rows*A.cols,'point':['2','1','3']}
        counts=[sum(x.K==K for x in ca) for K in range(4)]
        chans=[sum(x.channel==ch for x in ca) for ch in ('f','g','h')]
        k3=[sum(x.K==3 and x.channel==ch for x in ca) for ch in ('f','g','h')]
        if counts!=[2,6,10,14] or chans!=[8,8,16] or k3!=[3,4,7]:
            raise AssertionError('catalogue_counts')
        out[f'{species}.channel_counts']={'rank_counts':counts,'channels':chans,'K3_channels':k3}
    # Exact symbolic mass/momentum comparison, not just one specialization.
    from transverse_foundations import X,Y,M
    for species in ('quark','gluon'):
        labels=catalogue(species)
        cases=0
        for label in labels:
            for part in range(1 if label.m==0 else 2):
                left=cartesian_value(label,part,X,Y,M)
                right=helicity_value(label,part,X,Y,M)
                if any(S.cancel(a-b)!=0 for a,b in zip(left,right)):
                    raise AssertionError(f'{species}.symbolic_covariant {label.id()} part={part}')
                cases+=1
        out[f'{species}.symbolic_covariants']={'component_cases':cases,'domain':'Q(kx,ky,M_A), M_A>0'}
        _,zero_matrix=map_matrix(species,S.Integer(0),S.Integer(0),S.Integer(3))
        if any(zero_matrix[:,i]!=S.zeros(64,1) for i,label in enumerate(labels) if label.n>0):
            raise AssertionError(f'{species}.zero_momentum')
        out[f'{species}.zero_momentum']={'surviving_orbital_rank_zero':sum(label.n==0 for label in labels)}
        _,other=map_matrix(species,S.Integer(1),S.Integer(2),S.Integer(5))
        if other.rank()!=32:raise AssertionError(f'{species}.second_exact_point')
        out[f'{species}.second_exact_point']={'point':['1','2','5'],'rank':32}
    from spin_foundations import hermitian_basis
    ops=target_component_operators()
    transform=S.Matrix([[S.trace(H*ops[coord]) for coord in coordinates()]
                        for H in hermitian_basis()])
    if transform.rank()!=16:raise AssertionError('target.moment_coordinates')
    out['target.moment_coordinates']={'dimension':16,'method':'actual Cartesian target components on complete Hermitian matrix basis'}
    return out


def target_component_operators():
    """Operators whose expectation values are the source's A_Km components."""
    from spin_foundations import cartesian, I4
    J,Q,O,_,_=cartesian()
    return {
        (0,0,0):I4,
        (1,0,0):J[2],(1,1,0):J[0],(1,1,1):J[1],
        (2,0,0):Q[2,2],(2,1,0):2*Q[2,0],(2,1,1):2*Q[2,1],
        (2,2,0):Q[0,0]-Q[1,1],(2,2,1):2*Q[0,1],
        (3,0,0):O[2,2,2],(3,1,0):O[2,2,0],(3,1,1):O[2,2,1],
        (3,2,0):2*(O[2,0,0]-O[2,1,1]),(3,2,1):4*O[2,0,1],
        (3,3,0):O[0,0,0]-3*O[0,1,1],
        (3,3,1):3*O[0,0,1]-O[1,1,1],
    }


def joint_covariant(species,column,kx=S.Integer(2),ky=S.Integer(0),mass=S.Integer(3)):
    from sympy.physics.matrices import msigma
    labels=helicity_catalogue(species)
    label=labels[column]
    ops=target_component_operators()
    pauli=[S.eye(2),msigma(3),msigma(1),msigma(2)] if species=='quark' else [S.eye(2),msigma(2),msigma(3),msigma(1)]
    W=S.zeros(8)
    for K,m,part in coordinates():
        vals=helicity_value(label,part,kx,ky,mass) if (K,m)==(label.K,label.m) else (0,0,0,0)
        W+=sum((S.kronecker_product(vals[a]*pauli[a],ops[K,m,part]) for a in range(4)),S.zeros(8))
    return W


def parity_bound(species):
    """Independent unitary reflection on the actual parton-target space."""
    from sympy.physics.wigner import wigner_d_small
    from sympy.physics.matrices import msigma
    uy=wigner_d_small(S.Rational(3,2),-S.pi)
    up=-S.I*msigma(2) if species=='quark' else msigma(3)
    U=S.kronecker_product(up,uy)
    if species=='gluon':U=S.I*U
    if U**2 != S.eye(8):raise AssertionError('parity involution')
    pplus=(S.eye(8)+U)/2;pminus=(S.eye(8)-U)/2
    dimensions=(int(S.trace(pplus)),int(S.trace(pminus)))
    if dimensions!=(4,4):raise AssertionError('parity eigenspaces')
    for col in range(32):
        W=joint_covariant(species,col)
        if W.H!=W or U*W*U.H!=W:
            raise AssertionError(f'{species}.parity_invariance column={col}')
    return {'eigenspaces':list(dimensions),'real_invariant_dimension':sum(d*d for d in dimensions),
            'covariants_checked':32,'method':'independent unitary reflection on 8x8 Hermitian space'}


def rotation_covariance(species):
    from spin_foundations import HELICITIES
    target=S.diag(*(S.exp(-S.I*S.pi*m/2) for m in HELICITIES))
    parton=S.diag(S.exp(-S.I*S.pi/4),S.exp(S.I*S.pi/4)) if species=='quark' else S.Matrix([[0,-1],[1,0]])
    U=S.kronecker_product(parton,target)
    for col in range(32):
        before=joint_covariant(species,col,S.Integer(2),S.Integer(1),S.Integer(3))
        after=joint_covariant(species,col,S.Integer(-1),S.Integer(2),S.Integer(3))
        if S.simplify(U*before*U.H-after)!=S.zeros(8):
            raise AssertionError(f'{species}.rotation_covariance col={col}')
    return {'columns':32,'angle':'pi/2','point':['2','1','3'],'representation':'parton and target active unitary'}


def joint_expectation_check(species):
    """Contract independent Cartesian moments of complex states with W."""
    from sympy.physics.matrices import msigma
    labels=catalogue(species);ops=target_component_operators()
    A=S.Matrix([[1,S.I,0,1],[S.I,2,1,0],[0,1,1,S.I],[1,0,S.I,1]])
    rho=A*A.H/S.trace(A*A.H)
    if all(S.im(rho[i,j])==0 for i in range(4) for j in range(4)):
        raise AssertionError('complex coherence fixture is real')
    pauli=[S.eye(2),msigma(3),msigma(1),msigma(2)] if species=='quark' else [S.eye(2),msigma(2),msigma(3),msigma(1)]
    for col,label in enumerate(labels):
        FGH=[S.Integer(0)]*4
        for K,m,part in coordinates():
            if (K,m)!=(label.K,label.m):continue
            moment=S.trace(rho*ops[K,m,part])
            value=cartesian_value(label,part,S.Integer(2),S.Integer(1),S.Integer(3))
            FGH=[FGH[a]+moment*value[a] for a in range(4)]
        expected=sum((FGH[a]*pauli[a] for a in range(4)),S.zeros(2))
        joint=joint_covariant(species,col,S.Integer(2),S.Integer(1),S.Integer(3))
        actual=S.Matrix(2,2,lambda i,j:S.trace(rho*joint[4*i:4*(i+1),4*j:4*(j+1)]))
        if S.simplify(actual-expected)!=S.zeros(2):
            raise AssertionError(f'{species}.joint_expectation col={col}')
    return {'columns':32,'state':'positive exact complex coherence','index_order':'parton then target'}
