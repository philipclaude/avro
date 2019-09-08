// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include <vector>
#include <memory> //std::shared_ptr

#include "MPI/communicator_fwd.h"

#ifndef GLOBALCONTINUOUSMAP_H
#define GLOBALCONTINUOUSMAP_H

namespace SANS
{

struct GlobalContinuousMap
{
  GlobalContinuousMap() : nDOFpossessed(0), nDOFghost(0), nDOF_rank_offset(0) {}
  explicit GlobalContinuousMap(const std::shared_ptr<mpi::communicator>& comm)
    : comm(comm), nDOFpossessed(0), nDOFghost(0), nDOF_rank_offset(0) {}

  GlobalContinuousMap ( const GlobalContinuousMap & ) = default;
  GlobalContinuousMap ( GlobalContinuousMap && ) = default; //make sure there is a move constructor

  std::shared_ptr<mpi::communicator> comm;
  int nDOFpossessed;
  int nDOFghost;
  int nDOF_rank_offset;
  std::vector<int> remoteGhostIndex;
};

} // namespace SANS

#endif //GLOBALCONTINUOUSMAP_H
