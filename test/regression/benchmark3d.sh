#!/bin/bash

export avro=/Users/pcaplan/Codes/gitlab/mach/avro/build/release/bin/avro

$avro -adapt CKF-3-3-3 box Linear-3d cl rho=false write_conformity=true > cl.txt &
$avro -adapt cube-cylinder.mesh cube-cylinder.egads Linear-3d ccl rho=false write_conformity=true > ccl.txt &
$avro -adapt cube-cylinder.mesh cube-cylinder.egads Polar1 ccp1 rho=false write_conformity=true > ccp1.txt &
$avro -adapt cube-cylinder.mesh cube-cylinder.egads Polar2 ccp2 rho=false write_conformity=true > ccp2.txt &
