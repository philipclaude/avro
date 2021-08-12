import sympy

from sympy import diff, exp, sqrt, simplify, expand
from sympy.abc import x, y, z, t, pi
from sympy.printing.c import ccode

r = 0.4 + 0.4 * t
u = exp( -20 * ( r*r - x*x + y*y )**2 )

ux = diff( u , x )
uy = diff( u , y )
ut = diff( u , t )

uxx = diff( ux , x )
uxy = diff( ux , y )
uxt = diff( ux , t )

uyy = diff( uy , y )
uyt = diff( uy , t )

utt = diff( ut , t )

h = [ [uxx,uxy,uxy] , [uxy,uyy,uyt] , [uxt,uyt,utt] ]

print('real_t uxx = ' + ccode(uxx) + ';' )
print('real_t uxy = ' + ccode(uxy)+';')
print('real_t uxt = ' + ccode(uxt)+';')

print('real_t uyy = ' + ccode(uyy)+';')
print('real_t uyt = ' + ccode(uyt)+';')

print('real_t utt = ' + ccode(utt)+';')
