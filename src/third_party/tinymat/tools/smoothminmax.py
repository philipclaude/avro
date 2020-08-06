
import pylab as pyl
import numpy as npy

def smoothmax(x, y, a):
  m = npy.maximum(x,y);
  return m + (1/a)*npy.log( npy.exp(a*(x-m)) + npy.exp(a*(y-m)) );

def smoothmin(x, y, a):
  m = npy.minimum(x,y);
  return m - (1/a)*npy.log( npy.exp(-a*(x-m)) + npy.exp(-a*(y-m)) );

s1 = npy.linspace(-1,1,500)
s2 = npy.linspace(-1,1,11)

a = 40.
r = []
for x in s2:
  pyl.plot( s1, smoothmin(s1, x, a), linewidth=2 );

pyl.legend( ["s2 = " + str(x) for x in s2] )
pyl.xlabel("s1")
pyl.ylabel("smoothmin(s1,s2)")

pyl.figure(2)
for x in s2:
  pyl.plot( s1, smoothmax(s1, x, a), linewidth=2 );

pyl.legend( ["s2 = " + str(x) for x in s2] )
pyl.xlabel("s1")
pyl.ylabel("smoothmax(s1,s2)")

pyl.figure(3)
a1 = npy.linspace(10,100,10)
for a in a1:
  pyl.plot( s1, smoothmax(0, s1, a), linewidth=2 );

pyl.legend( ["a = " + str(x) for x in a1] )
pyl.xlabel("s1")
pyl.ylabel("smoothmax(0,s1)")


pyl.show()
