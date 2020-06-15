#!/bin/bash

WORKSPACE=$(git rev-parse --show-toplevel)

source $WORKSPACE/util/environment.sh

cmake $CMAKEARGS -DAVRO_WITH_MPI=ON \
      $@ \
      $WORKSPACE
