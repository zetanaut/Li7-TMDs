"""String-input arbitrary-precision direct Born amplitudes (mpmath 1.3.0)."""
from __future__ import annotations
import mpmath as mp


def _gamma():
    z=mp.mpc(0);o=mp.mpc(1);i=mp.j
    sx=mp.matrix([[z,o],[o,z]])
    sy=mp.matrix([[z,-i],[i,z]])
    sz=mp.matrix([[o,z],[z,-o]])
    eye=mp.eye(2);zero=mp.zeros(2)
    def block(a,b,c,d):
        return mp.matrix([[a[r,s] if s<2 else b[r,s-2] for s in range(4)] for r in range(2)]
                        +[[c[r-2,s] if s<2 else d[r-2,s-2] for s in range(4)] for r in range(2,4)])
    return (block(eye,zero,zero,-eye),)+tuple(block(zero,a,-a,zero) for a in (sx,sy,sz)),(sx,sy,sz)


def _slash(v,g):
    return sum((sgn*v[a]*g[a] for a,sgn in enumerate((1,-1,-1,-1))),mp.zeros(4))


def _dot(a,b):return a[0]*b[0]-sum(a[j]*b[j] for j in (1,2,3))


def _vectors(values):
    root,Q2,m,theta,phi,E=values
    if Q2<=0 or m<=0 or root<=2*m or E<=root/2:
        raise ValueError('interior high-precision Born domain')
    omega=(root*root+Q2)/(2*root)
    q0=(root*root-Q2)/(2*root)
    k=[omega,0,0,omega];q=[q0,0,0,-omega]
    p=mp.sqrt(root*root/4-m*m)
    p1=[root/2,p*mp.sin(theta)*mp.cos(phi),p*mp.sin(theta)*mp.sin(phi),p*mp.cos(theta)]
    p2=[root/2,-p1[1],-p1[2],-p1[3]]
    lz=(-Q2/2-E*q0)/omega
    lt2=E*E-lz*lz
    if lt2<=0 or E-q0<=0:raise ValueError('no massless lepton solution')
    l=[E,mp.sqrt(lt2),0,lz]
    lp=[l[j]-q[j] for j in range(4)]
    return l,lp,k,q,p1,p2


def _column(values):return mp.matrix([[a] for a in values])


def _spinors(p,m,sig,anti=False):
    a=mp.sqrt(p[0]+m)
    pd=sum((p[j+1]*sig[j] for j in range(3)),mp.zeros(2))
    out=[]
    for chi in (_column((1,0)),_column((0,1))):
        low=pd*chi/a
        upper=low if anti else a*chi
        lower=a*chi if anti else low
        out.append(_column((upper[0],upper[1],lower[0],lower[1])))
    return out


def _lepton(p,h):
    theta=mp.atan2(p[1],p[3]) # the test CM lepton plane is x-z
    chi=(_column((mp.cos(theta/2),mp.sin(theta/2))) if h==1 else
         _column((-mp.sin(theta/2),mp.cos(theta/2))))
    a=mp.sqrt(p[0]);return _column((a*chi[0],a*chi[1],h*a*chi[0],h*a*chi[1]))


def evaluate_strings(*,sqrt_s='5',Q2='4',mass='1.5',theta='0.8',phi='0.4',
                     lepton_energy='10',spinor_helicity=-1,dps=50):
    if spinor_helicity not in (-1,1):raise ValueError('definite spinor helicity required')
    with mp.workdps(dps):
        values=[mp.mpf(v) for v in (sqrt_s,Q2,mass,theta,phi,lepton_energy)]
        l,lp,k,q,p1,p2=_vectors(values)
        g,sig=_gamma();m=values[2]
        ein=_lepton(l,spinor_helicity);eout=_lepton(lp,spinor_helicity)
        current=[(eout.conjugate().transpose()*g[0]*g[mu]*ein)[0] for mu in range(4)]
        us=_spinors(p1,m,sig);vs=_spinors(p2,m,sig,True)
        d1=_dot([p1[j]-q[j] for j in range(4)],[p1[j]-q[j] for j in range(4)])-m*m
        d2=_dot([p1[j]-k[j] for j in range(4)],[p1[j]-k[j] for j in range(4)])-m*m
        numerator1=_slash([p1[j]-q[j] for j in range(4)],g)+m*mp.eye(4)
        numerator2=_slash([p1[j]-k[j] for j in range(4)],g)+m*mp.eye(4)
        amplitudes=mp.zeros(2,4)
        for a in range(2):
            for si in range(2):
                for sj in range(2):
                    value=mp.mpc(0)
                    for mu in range(4):
                        chain=g[mu]*numerator1*g[a+1]/d1+g[a+1]*numerator2*g[mu]/d2
                        value+=(-1 if mu else 1)*current[mu]*(us[si].conjugate().transpose()*g[0]*chain*vs[sj])[0]
                    amplitudes[a,2*si+sj]=value
        B=amplitudes*amplitudes.conjugate().transpose()
        stokes=((B[0,0]+B[1,1])/2,(B[0,1]-B[1,0])/(2j),
                (B[0,0]-B[1,1])/2,(B[0,1]+B[1,0])/2)
        return {'dps':dps,'inputs':list(map(str,(sqrt_s,Q2,mass,theta,phi,lepton_energy))),
                'spinor_helicity':spinor_helicity,
                'stokes':[mp.nstr(mp.re(x),dps-5) for x in stokes],
                'B_real':[[mp.nstr(mp.re(B[i,j]),dps-5) for j in range(2)] for i in range(2)],
                'B_imag':[[mp.nstr(mp.im(B[i,j]),dps-5) for j in range(2)] for i in range(2)]}
