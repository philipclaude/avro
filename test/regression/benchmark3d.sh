#!/bin/bash

export avro=$1

$avro -adapt CKF-3-3-3 box Linear-3d $3/cl rho=false write_conformity=true > cl.txt &
$avro -adapt $2/meshes/cube-cylinder.mesh $2/geometry/cube-cylinder.egads Linear-3d $3/ccl rho=false write_conformity=true > ccl.txt &
$avro -adapt $2/meshes/cube-cylinder.mesh $2/geometry/cube-cylinder.egads Polar1 $3/ccp1 rho=false write_conformity=true > ccp1.txt &
$avro -adapt $2/meshes/cube-cylinder.mesh $2/geometry/cube-cylinder.egads Polar2 $3/ccp2 rho=false write_conformity=true > ccp2.txt &
