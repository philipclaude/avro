
import pylab as pyl
import numpy as npy

def smoothActivation_exp(x, a, eps):
  b = -npy.log(npy.expm1(a*eps))/a
  return x - b + (npy.log1p( npy.exp(-a*(x - b)) ) - npy.log1p( npy.exp(-a*(1-(x - b))) ))/a;

x = npy.linspace(-0.05,1.05,500)
a1 = npy.linspace(100,1000,7)
eps1 = npy.logspace(-10,-1,7)

eps = 1e-12
for a in a1:
  pyl.plot( x, smoothActivation_exp(x, a, eps), linewidth=2 );

pyl.legend( ["a = " + str(a) for a in a1] )
pyl.xlabel("x")
pyl.ylabel("smoothActivation_exp(x,a,"+str(eps)+")")

pyl.figure(2)
a = 500
for eps in eps1:
  pyl.plot( x, smoothActivation_exp(x, a, eps), linewidth=2 );
 
pyl.legend( ["eps = " + str(eps) for eps in eps1] )
pyl.xlabel("x")
pyl.ylabel("smoothActivation_exp(x,"+str(a)+",eps)")

pyl.figure(3)
x = npy.linspace(-0.5,1.5,500)

eps = 1e-10
a = 500
pyl.plot( x, smoothActivation_exp(x, a, eps), linewidth=2 );

pyl.xlabel("x")
pyl.ylabel("smoothActivation_exp(x,"+str(a)+","+str(eps)+")")

pyl.figure(4)
pyl.plot( x, -a*(x + npy.log(npy.exp(a*eps) - 1)/a), linewidth=2 );
pyl.plot( x, -a*(1-x), linewidth=2 );

b = -npy.log(npy.exp(a*eps) - 1)/a
x = -0.25
print -a*(x + b), npy.exp(-a*(x + b)), smoothActivation_exp(x, a, eps)
print -a*(1-x), npy.exp(-a*(1-x)), smoothActivation_exp(x, a, eps)
print
x = 1.25
print -a*(x + b), npy.exp(-a*(x + b)), smoothActivation_exp(x, a, eps)
print -a*(1-x), npy.exp(-a*(1-x)), smoothActivation_exp(x, a, eps)
print
x = 0.5
print -a*(x + b), npy.exp(-a*(x + b)), (npy.log1p( npy.exp(-a*(x - b)) ))/a;
print -a*(1-x), npy.exp(-a*(1-x)), (-npy.log1p( npy.exp(-a*(1-x)) ))/a;

pyl.show()
