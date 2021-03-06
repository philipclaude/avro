#!/bin/bash

workspace=`pwd`
config=$1
nproc=$2


# files might linger if a build was aborted
echo "Removing any lingering untracked files"
for f in `git ls-files --others --exclude-standard`; do
  rm -f $f
done

# update submodules
./get_submodules.sh

# create the build directory if it does not exist
cmakedir=build/$config
mkdir -p $cmakedir
cd $cmakedir

CMAKE_ARGS="-DMACHII_LIBRARY_LOCATION='/home/gitlab-runner/Codes/mach-II/library' -DAVRO_HEADLESS_GRAPHICS=1"

if [[ $config == *"coverage"* ]]; then
  CMAKE_ARGS="$CMAKE_ARGS  -DAVRO_WITH_MPI=ON"
fi

# there is a bug in gcc 4.8 which prevents from using some features needed by mpi wrapper
if [[ $config == *"gnu48"* ]]; then
  CMAKE_ARGS="$CMAKE_ARGS -DUSE_MPI=OFF"
fi

if [[ $config == *"no_esp"* ]]; then
  CMAKE_ARGS="$CMAKE_ARGS -DAVRO_NO_ESP=ON"
fi

if [[ $config == *"mpi"* ]]; then
  CMAKE_ARGS="$CMAKE_ARGS -DAVRO_WITH_MPI=ON"
fi

time source $workspace/util/configure.sh $CMAKE_ARGS

if [[ $config == *"coverage"* ]]; then

  # coverage information always needs to be completely recompiled
  time make coverage_cleaner
fi

# build all libraries and executables
time make -j $nproc

if [[ $config == *"coverage"* ]]; then

  Xvfb :5 -screen 0 1280x1280x24 &
  export DISPLAY=:5

  time make unit_coverage
  time make coverage_info
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
