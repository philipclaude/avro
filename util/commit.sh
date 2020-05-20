#!/bin/bash

workspace=`pwd`
config=$1


#Files might linger if a build was aborted
echo "Removing any lingering untracked files"
for f in `git ls-files --others --exclude-standard`; do
  rm -f $f
done

#Create the build directory if it does not exist
cmakedir=build/$config
mkdir -p $cmakedir
cd $cmakedir


if [[ $config == *"coverage"* ]]; then
  CMAKE_ARGS=""
else
  CMAKE_ARGS=""
fi
if [[ $config == *"memcheck"* ]]; then
  CMAKE_ARGS="$CMAKE_ARGS -DUSE_MPI=OFF"
fi

# there is a bug in gcc 4.8 which prevents from using some features needed by mpi wrapper
if [[ $config == *"gnu48"* ]]; then
  CMAKE_ARGS="$CMAKE_ARGS -DUSE_MPI=OFF"
fi

time source $workspace/util/jenkins/configure.sh $CMAKE_ARGS

# number of processors used in compile
nproc=1

if [[ $config == *"coverage"* ]]; then
  # just because we have to clean, use more processors for coverage
  nproc=6

  # coverage information always needs to be completely recompiled
  time make coverage_cleaner
fi

# build all libraries and executables
time make -j $nproc

if [[ $config == *"coverage"* ]]; then

  time make unit_coverage

  time make coverage_info
  if [[ `hostname` == "wazowski" ]]; then
    #time python $HOME/util/lcov/lcov_cobertura.py coverage.info -o $WORKSPACE/build/coverage.xml
  fi
  # generate html documents
  #time make genhtml
  # udpdate coverage information on acdl (@marshall do i need to rename this?)
  #time rsync -rtuc --delete avroCoverageHTML jenkins@acdl.mit.edu:/var/lib/jenkins/userContent/

fi

# warning if any files were generated.
# count the number of files
cd $workspace
tmpfiles=`git ls-files --others --exclude=build --exclude=app/WebViewer | wc -l`

if [ $tmpfiles -ne 0 ]; then
  echo "warning: Files should not be generated in code pushed to the repository."
  echo "warning: Files found :"
  for f in `git ls-files --others --exclude=build --exclude=external`; do
    echo "warning: $f"
  done
fi
