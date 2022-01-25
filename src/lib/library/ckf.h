//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_LIBRARY_CKF_H_
#define avro_LIB_LIBRARY_CKF_H_

#include "avro_types.h"

#include "mesh/points.h"
#include "mesh/topology.h"

namespace avro
{

class CKF_Triangulation : public Topology<Simplex>
{
public:
  CKF_Triangulation( const std::vector<index_t>& dims );

  void generate();

private:
  void ndgrid( const std::vector<real_t>& dx , coord_t current );
  index_t find_vertex( const real_t* y ) const;

  Points points_;
  std::vector<real_t> p_;
  std::vector<index_t> u_;
  std::vector<index_t> dims_;

  Table<index_t> grid_;
};

template<typename type>
class CubeDomain : public Topology<type>
{
public:
  CubeDomain( coord_t number , coord_t dim , index_t n=10 );

private:
  Points points_;
};

} // avro

#endif
