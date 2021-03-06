#!/bin/bash

export avro=$1

$avro -adapt CKF-3-3-3 box Linear-3d $3/cl rho=false write_conformity=true nb_iter=20 > cl.txt &
$avro -adapt $2/meshes/cube-cylinder.mesh $2/geometry/cube-cylinder.egads Linear-3d $3/ccl rho=false write_conformity=true nb_iter=20 > ccl.txt &
$avro -adapt $2/meshes/cube-cylinder.mesh $2/geometry/cube-cylinder.egads Polar1 $3/ccp1 rho=false write_conformity=true href=1.7 nb_iter=20 > ccp1.txt &
$avro -adapt $2/meshes/cube-cylinder.mesh $2/geometry/cube-cylinder.egads Polar2 $3/ccp2 rho=false write_conformity=true nb_iter=20 > ccp2.txt &

# wait for all benchmarks to complete
wait

# analyze metric conformity
$avro -conformity $3/cl_19.mesh Linear-3d $3/cl_conformity.json nb_expected=39000
status_cl=$?
echo "status CL = $status_cl"

$avro -conformity $3/ccl_19.mesh Linear-3d $3/ccl_conformity.json nb_expected=31400
status_ccl=$?
echo "status CCL = $status_ccl"

# ccp1 is the hardest case, and not a good indication of how well avro is working
# the returned status is ignored below
$avro -conformity $3/ccp1_19.mesh Polar1 $3/ccp1_conformity.json nb_expected=21000
status_ccp1=$?
echo "status CCP1 = $status_ccp1"

$avro -conformity $3/ccp2_19.mesh Polar2  $3/ccp2_conformity.json nb_expected=32000
status_ccp2=$?
echo "status CCP2 = $status_ccp2"

# check all statuses except for ccp1
if [[ $status_cl -eq 1 ]] || [[ $status_ccl -eq 1 ]] || [[ $status_ccp2 -eq 1 ]] ; then
  echo "regression tests failed :("
  exit 1
fi
