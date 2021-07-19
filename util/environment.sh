#!/bin/bash

shopt -s nocasematch

dir=`pwd`
builddir=`basename $dir`

if [[ `hostname` == *"wazowski"* ]]; then

  if [[ $builddir == *"intel14"* ]]; then
    source /opt/intel/composer_xe_2013_sp1/bin/compilervars.sh intel64 > iccvars.out 2>&1
  elif [[ $builddir == *"intel15"* ]]; then
    source /opt/intel/composer_xe_2015/bin/compilervars.sh intel64 > iccvars.out 2>&1
  elif [[ $builddir == *"intel17"* ]]; then
    source /opt/intel/compilers_and_libraries_2017/linux/bin/compilervars.sh intel64 > iccvars.out 2>&1
  elif [[ $builddir == *"intel18"* ]]; then
    source /opt/intel/compilers_and_libraries_2018/linux/bin/compilervars.sh intel64 > iccvars.out 2>&1
  #else
  #  source /opt/intel/bin/iccvars.sh intel64 > iccvars.out 2>&1
  fi

  # this makes sure cmake finds open-mpi rather than intel mpi
  export PATH=/usr/bin:$PATH
  export LD_LIBRARY_PATH=/usr/lib:$LD_LIBRARY_PATH

  #export LAPACK_DIR=/home/gitlab-runner/lapack/lapack-3.6.1/install
  export ESP_DIR=/home/gitlab-runner/Codes/EngSketchPad
  export CAS_DIR=/home/gitlab-runner/Codes/OpenCASCADE-7.3.1
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CAS_DIR/lib

  #export PARMETIS_DIR=/home/gitlab-runner/parmetis/parmetis-4.0.3/install
  #export METIS_DIR=/home/gitlab-runner/parmetis/metis-5.1.0/install
  #export PATH=/home/gitlab-runner/lcov/lcov-1.13/bin:$PATH

elif [[ `hostname` == *"midd-19641"* ]]; then
  # hope that everything is fine
  echo "assuming libraries are installed"
else
  if [ -z "$LAPACK_DIR" ]; then
    echo "Please set LAPACK_DIR in your environment."
    exit 0
  fi
fi

env > env.log 2>&1
