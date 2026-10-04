"""Source-indexed, conditional parton-target positivity calculations."""
from __future__ import annotations
from itertools import combinations
import sympy as s
from spin_foundations import cartesian
from correlator_foundations import catalogue,joint_covariant

R=s.Rational
J,Q,O,_,_=cartesian()
I4=s.eye(4)

def source_collinear_matrix(species, f, fq, g, go, h=0, ho=0):
    """Literal eq:jointM from the source's defining target operators."""
    F=f*I4+fq*Q[2,2]
    G=g*J[2]+go*O[2,2,2]
    if species=='quark':
        Tx=h*J[0]+ho*O[2,2,0]
        Ty=h*J[1]+ho*O[2,2,1]
        return s.BlockMatrix([[(F+G)/2,(Tx+s.I*Ty)/2],
                              [(Tx-s.I*Ty)/2,(F-G)/2]]).as_explicit()
    if species=='gluon':
        Hxx=h*(Q[0,0]-Q[1,1])
        Hxy=2*h*Q[0,1]
        return s.BlockMatrix([[(F+Hxx)/2,(Hxy+s.I*G)/2],
                              [(Hxy-s.I*G)/2,(F-Hxx)/2]]).as_explicit()
    raise ValueError(species)

def source_fixed_matrix(species,F,G,X,Y):
    if species=='quark':
        M=s.Matrix([[F+G,X+s.I*Y],[X-s.I*Y,F-G]])/2
    elif species=='gluon':
        M=s.Matrix([[F+X,Y+s.I*G],[Y-s.I*G,F-X]])/2
    else:raise ValueError(species)
    return M

def exact_psd(M):
    if M!=M.H:raise ValueError('raw matrix is not Hermitian')
    values=M.eigenvals()
    if any(s.ask(s.Q.negative(v)) is True for v in values):return 'indefinite'
    if all(s.ask(s.Q.nonnegative(v)) is True for v in values):
        return 'boundary' if any(v==0 for v in values) else 'positive_semidefinite'
    return 'undetermined'

def fixed_target_checks():
    F,G,X,Y=s.symbols('F G X Y',real=True)
    rows=[]
    for species in ('quark','gluon'):
        M=source_fixed_matrix(species,F,G,X,Y)
        t=s.symbols('t')
        cp=s.factor(M.charpoly(t).as_expr())
        expected=s.expand(t*t-F*t+R(1,4)*(F*F-G*G-X*X-Y*Y))
        if s.expand(cp-expected)!=0:raise AssertionError('fixed-target characteristic polynomial')
        cases={
          'interior':(s.Integer(3),s.Integer(1),s.Integer(1),s.Integer(1)),
          'boundary':(s.Integer(1),s.Integer(1),s.Integer(0),s.Integer(0)),
          'violating':(s.Integer(1),s.Integer(2),s.Integer(0),s.Integer(0)),
          'negative_trace':(s.Integer(-3),s.Integer(1),s.Integer(1),s.Integer(1)),
        }
        for name,values in cases.items():
            verdict=exact_psd(M.subs(dict(zip((F,G,X,Y),values))))
            wanted={'interior':'positive_semidefinite','boundary':'boundary',
                    'violating':'indefinite','negative_trace':'indefinite'}[name]
            if verdict!=wanted:raise AssertionError(f'fixed-target {species} {name}')
            rows.append({'species':species,'case':name,'values':list(map(str,values)),
                         'predicate':verdict})
    return {'characteristic_polynomial':str(expected),'cases':rows,
            'inequality':'F >= sqrt(G**2+X**2+Y**2), including F>=0',
            'source_labels':['eq:qD','eq:gD','eq:fixedrhobounds']}

def _components(M):
    unseen=set(range(M.rows)); groups=[]
    while unseen:
        group={unseen.pop()};front=list(group)
        while front:
            i=front.pop()
            for j in list(unseen):
                if M[i,j]!=0 or M[j,i]!=0:
                    unseen.remove(j);group.add(j);front.append(j)
        groups.append(tuple(sorted(group)))
    return tuple(sorted(groups,key=lambda x:(-len(x),x)))

def collinear_blocks(species):
    f,fq,g,go,h,ho=s.symbols('f f_Q g g_O h h_O',real=True)
    M=source_collinear_matrix(species,f,fq,g,go,h,ho)
    if species=='gluon':
        B=s.Matrix([[1,1],[-s.I,s.I]])/s.sqrt(2)
        U=s.kronecker_product(B,I4)
        if U.H*U!=s.eye(8):raise AssertionError('circular basis unitary')
        M=s.simplify(U.H*M*U)
    else:U=s.eye(8)
    M=s.simplify(M)
    groups=_components(M)
    order=sum((list(gp) for gp in groups),[])
    perm=s.eye(8)[:,order]
    block=s.simplify(perm.T*M*perm)
    if any(block[i,j]!=0 for i in range(8) for j in range(8)
           if next(z for z,gp in enumerate(groups) if i<sum(len(x) for x in groups[:z+1]))!=
              next(z for z,gp in enumerate(groups) if j<sum(len(x) for x in groups[:z+1]))):
        raise AssertionError('collinear block permutation')
    Uo,Ui=f+fq,f-fq
    Go,Gi=R(3,2)*g+R(3,10)*go,g/2-R(9,10)*go
    if species=='quark':
        expected=(3*(h+R(2,5)*ho)**2,(Uo+Go)*(Ui-Gi),
                  (2*h-R(6,5)*ho)**2,(Ui+Gi)**2)
    else:
        expected=(12*h*h,(Uo+Go)*(Ui+Gi))
    determinants=[s.factor(M.extract(gp,gp).det()) for gp in groups if len(gp)==2]
    # Match the separately transcribed source inequalities through the
    # actual block determinants (their matrix entries include factor 1/2).
    if species=='quark':
        required=[s.factor((expected[1]-expected[0])/4),
                  s.factor((expected[3]-expected[2])/4),
                  s.factor((expected[1]-expected[0])/4)]
    else:required=[s.factor((expected[1]-expected[0])/4)]*2
    if sorted(map(str,determinants))!=sorted(map(str,required)):
        raise AssertionError(f'{species} source collinear bound factor mismatch: {determinants} vs {required}')
    return {'species':species,'basis':[[str(U[i,j]) for j in range(8)] for i in range(8)],
            'groups':[list(x) for x in groups],'permutation':order,
            'blocks':[[[str(M[i,j]) for j in gp] for i in gp] for gp in groups],
            'determinants':list(map(str,determinants)),
            'source_labels':['eq:jointM','eq:UG','eq:diagbound',
                             'eq:qSoffer' if species=='quark' else 'eq:gSoffer']}

def matrix_algorithm_counterexamples():
    C=s.Matrix([[1,-R(3,4),-R(3,4)],[-R(3,4),1,-R(3,4)],[-R(3,4),-R(3,4),1]])
    Z=s.Matrix([[0,0,0],[0,1,2],[0,2,1]])
    pair=[C.extract(index,index).det() for index in combinations(range(3),2)]
    leading=[Z[:i,:i].det() for i in range(1,4)]
    if set(C.eigenvals())!={-R(1,2),R(7,4)} or pair!=[R(7,16)]*3 or C.det()!=-R(49,32):
        raise AssertionError('pairwise-minor counterexample')
    if leading!=[0,0,0] or Z.extract((1,2),(1,2)).det()!=-3:
        raise AssertionError('leading-minor counterexample')
    boundary=s.diag(1,0,2)
    if exact_psd(boundary)!='boundary':raise AssertionError('singular PSD boundary')
    return {'pairwise_spectrum':{str(k):v for k,v in C.eigenvals().items()},
            'pairwise_minors':list(map(str,pair)),'pairwise_determinant':str(C.det()),
            'leading_minors':list(map(str,leading)),
            'nonleading_negative_minor':'-3','singular_positive':'boundary'}

def joint_gram_checks():
    # Conjugate-amplitude-first convention of eq:Gram.  Complex entries
    # deliberately expose a partial transpose or discarded coherence.
    A=s.Matrix(3,8,lambda row,col:
       R((row+2)*(col+1)%7-3,5)+s.I*R((2*row+col)%5-2,7))
    M=A.H*A
    psi=s.Matrix([1,s.I,2,1-s.I])
    rho=psi*psi.H/(psi.H*psi)[0]
    partial=s.Matrix(2,2,lambda a,b:s.trace(rho*M[4*a:4*a+4,4*b:4*b+4]))
    amplitude=s.Matrix(3,2,lambda x,a:sum(A[x,4*a+i]*psi[i] for i in range(4)))
    expected=amplitude.H*amplitude/(psi.H*psi)[0]
    if s.simplify(partial-expected)!=s.zeros(2):
        raise AssertionError('spectral target contraction index order')
    if s.simplify(partial-s.Matrix(2,2,lambda a,b:
           s.trace(rho.T*M[4*a:4*a+4,4*b:4*b+4])))==s.zeros(2):
        raise AssertionError('complex target coherence failed to detect transpose')
    # Exact factor A.H*A certifies PSD for every complex vector without
    # expanding an eighth-degree characteristic polynomial.
    if M!=M.H or s.simplify(M-A.H*A)!=s.zeros(8):
        raise AssertionError('spectral Gram factorization')
    return {'amplitudes':'3 exact rational-complex rows on 8 parton-target indices',
            'Gram_rank':M.rank(),'partial_target_state':'exact complex pure rho',
            'partial_spectrum':[str(v) for v in partial.eigenvals()],
            'index_formula':'M_(ai,bj)=sum_X conjugate(A_Xai)*A_Xbj; D_ab=Tr(rho M_ab)',
            'partial_transpose_detected':True,'source_labels':['eq:jointM','eq:Gram']}

def collinear_witnesses():
    rows=[]
    examples={
      'quark':{
       'interior':(3,0,0,0,R(1,2),0),
       'middle_boundary':(1,0,0,0,R(1,2),0),
       'zero_boundary':(0,0,0,0,0,0),
       'outer_violation':(1,0,0,0,1,R(5,3)),
       'middle_violation':(1,0,0,0,-2,5),
       'diagonal_violation':(1,0,2,0,0,0),
      },
      'gluon':{
       'interior':(3,0,0,0,R(1,2),0),
       'saturated':(1,0,0,0,s.sqrt(3)/6,0),
       'zero_boundary':(0,0,0,0,0,0),
       'flip_violation':(1,0,0,0,1,0),
       'diagonal_violation':(1,0,2,0,0,0),
      }}
    for species,cases in examples.items():
        for name,values in cases.items():
            M=s.simplify(source_collinear_matrix(species,*values))
            verdict=exact_psd(M)
            wanted='indefinite' if 'violation' in name else 'boundary' if name in ('zero_boundary','middle_boundary','saturated') else 'positive_semidefinite'
            if verdict!=wanted:raise AssertionError(f'{species} {name}: {verdict} vs {wanted}')
            rows.append({'species':species,'case':name,'coefficients':list(map(str,values)),
                         'predicate':verdict,'minimum_eigenvalue':str(min(M.eigenvals(),key=lambda x:float(x)))})
    return {'cases':rows,'case_count':len(rows),'scheme_domain':'conditional positive spectral/input prescription'}

def convention_reproducer():
    """Expose the parton-index transposition in the old basis explicitly."""
    for species,label_id in (('quark','quark.h.11[0]'),('gluon','gluon.g.10[0]')):
        col=next(i for i,x in enumerate(catalogue(species)) if x.id()==label_id)
        old=joint_covariant(species,col,0,0,1)
        if species=='quark':source=source_collinear_matrix(species,0,0,0,0,1,0)
        else:source=source_collinear_matrix(species,0,0,1,0,0,0)
        if old==source:raise AssertionError('old source convention discrepancy disappeared')
        # Existing map is twice the correlator and uses the opposite sigma2.
        if old==2*source:raise AssertionError('parton sigma2 sign discrepancy disappeared')
    return {'affected':['quark.h.11[0]','gluon.g.10[0]'],
            'source':'eq:qD, eq:gD, eq:jointM upper off-diagonal has +i',
            'existing_map':'joint_covariant uses +msigma(2), whose upper off-diagonal is -i',
            'status':'resolved by explicit source-index construction and checked parton-index conversion; original auxiliary map retained'}
