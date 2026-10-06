"""Reviewed source Fourier rows; no production-generator imports.

Tuple fields are (K,m,channel,n,coefficient,p,q,angular parity).  The
phase is p*phi_k+q*psi, with -2*phi_B also present for h rows.
"""
from fractions import Fraction

ROWS=(
 (0,0,'f',0,1,0,0,'one'),(0,0,'h',2,Fraction(1,2),2,0,'cos'),
 (1,0,'g',0,1,0,0,'one'),(1,0,'h',2,Fraction(1,2),2,0,'sin'),
 (1,1,'f',1,1,1,-1,'sin'),(1,1,'g',1,1,1,-1,'cos'),
 (1,1,'h',1,Fraction(1,2),1,1,'sin'),(1,1,'h',3,Fraction(1,4),3,-1,'sin'),
 (2,0,'f',0,1,0,0,'one'),(2,0,'h',2,Fraction(1,2),2,0,'cos'),
 (2,1,'f',1,1,1,-1,'cos'),(2,1,'g',1,1,1,-1,'sin'),
 (2,1,'h',1,Fraction(1,2),1,1,'cos'),(2,1,'h',3,Fraction(1,4),3,-1,'cos'),
 (2,2,'f',2,1,2,-2,'cos'),(2,2,'g',2,1,2,-2,'sin'),
 (2,2,'h',0,1,0,2,'cos'),(2,2,'h',4,Fraction(1,4),4,-2,'cos'),
 (3,0,'g',0,1,0,0,'one'),(3,0,'h',2,Fraction(1,2),2,0,'sin'),
 (3,1,'f',1,1,1,-1,'sin'),(3,1,'g',1,1,1,-1,'cos'),
 (3,1,'h',1,Fraction(1,2),1,1,'sin'),(3,1,'h',3,Fraction(1,4),3,-1,'sin'),
 (3,2,'f',2,1,2,-2,'sin'),(3,2,'g',2,1,2,-2,'cos'),
 (3,2,'h',0,1,0,2,'sin'),(3,2,'h',4,Fraction(1,4),4,-2,'sin'),
 (3,3,'f',3,1,3,-3,'sin'),(3,3,'g',3,1,3,-3,'cos'),
 (3,3,'h',1,1,-1,3,'sin'),(3,3,'h',5,Fraction(1,4),5,-3,'sin'),
)
assert len(ROWS)==32

# Independently transcribed tab:oct_fourier modes. The fifth field is the
# trigonometric function and the last is the coefficient of phi_B.
OCT_MODES=(
 ('g.30',1,0,0,'one',0),('h.30.2',0,2,0,'sin',-2),
 ('f.31',0,1,-1,'sin',0),('g.31',1,1,-1,'cos',0),
 ('h.31.1',0,1,1,'sin',-2),('h.31.3',0,3,-1,'sin',-2),
 ('f.32',0,2,-2,'sin',0),('g.32',1,2,-2,'cos',0),
 ('h.32.0',0,0,2,'sin',-2),('h.32.4',0,4,-2,'sin',-2),
 ('f.33',0,3,-3,'sin',0),('g.33',1,3,-3,'cos',0),
 ('h.33.1',0,-1,3,'sin',-2),('h.33.5',0,5,-3,'sin',-2),
)
