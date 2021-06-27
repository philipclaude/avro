# add the path to the pyavro module using sys
import sys
sys.path.insert(0,'../../../build/debug/lib/')

# import pyavro
import pyavro

number = 3
dim    = number
udim   = dim -1
ctx    = pyavro.Context(number,dim,udim)

ctx.define_geometry("box")
ctx.define_mesh("CKF-3-3-3")
ctx.attach_geometry()

coordinates,connectivity = pyavro.retrieve_mesh(ctx)

nb_points = int(len(coordinates)/dim)
nb_rank   = int(number*(number+1)/2)

metric = [ 0.0 for i in range(nb_points*nb_rank) ]
h      = 0.25
for k in range(nb_points):
    idx = k*nb_rank
    for i in range(number):
        metric[idx] = 1./h**2
        idx += (number -i)

ctx.adapt(metric)

coordinates,connectivity = pyavro.retrieve_mesh(ctx)
nb_points = int(len(coordinates)/dim)
nb_elems  = int(len(connectivity)/(number+1))
print('adapted mesh has nb_points = ',nb_points,'nb_elems = ',nb_elems)

faces,geometry = pyavro.retrieve_boundary(ctx)
assert( len(faces) == len(geometry) )
k = 0
for bgroup in faces:
    nb_faces_in_group = int(len(bgroup)/number)
    print('boundary on entity',geometry[k],'has ',nb_faces_in_group,'faces')
    k += 1
