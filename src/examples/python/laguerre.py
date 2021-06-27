# add the path to the pyavro module using sys
import random
import sys
sys.path.insert(0,'../../../build/release/lib/')

# import pyavro
import pyavro

number = 2
dim    = number
udim   = dim -1
ctx    = pyavro.Context(number,dim,udim)

nb_points = 10**4
x = [ random.uniform(0,1) for i in range(dim*nb_points) ]
w = [ 0.0 for i in range(nb_points) ]

x = ctx.compute_laguerre( x , w , 10 )
mass = [ 1./nb_points for i in range(nb_points) ]
w = ctx.compute_optimal_transport( x , mass , w , 10 )

vertices,polytopes,nv_per_elem = pyavro.retrieve_polytopes(ctx)

nb_vertices  = len(vertices)/dim
nb_polytopes = len(nv_per_elem)

'''
i = 0
for k in range(nb_polytopes):
    print('polytope[',k,']: (',end='')
    for j in range(nv_per_elem[k]):
        print(polytopes[i],end=' ')
        i += 1
    print(')')
'''

ctx.plot()
