//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef AVRO_LIB_ELEMENT_LAGRANGE_NODES_H_
#define AVRO_LIB_ELEMENT_LAGRANGE_NODES_H_

#include "common/error.h"
#include "avro_types.h"

namespace avro
{

template<typename type, int N, int P> class LagrangeNodes;

template<typename type,int P>
class LagrangeNodes<type,1,P> {
public:
  static const std::vector<real_t> coord_s_;
};

template<typename type,int P>
class LagrangeNodes<type,2,P> {
public:
  static const std::vector<real_t> coord_s_;
  static const std::vector<real_t> coord_t_;
};

template<typename type,int P>
class LagrangeNodes<type,3,P> {
public:
  static const std::vector<real_t> coord_s_;
  static const std::vector<real_t> coord_t_;
  static const std::vector<real_t> coord_u_;
};

template<typename type,int P>
class LagrangeNodes<type,4,P> {
public:
  static const std::vector<real_t> coord_s_;
  static const std::vector<real_t> coord_t_;
  static const std::vector<real_t> coord_u_;
  static const std::vector<real_t> coord_v_;
};

template<typename type,int N, int P>
struct GetLagrangeNodes {
  static void get( std::vector<real_t>& nodes );
};

template<typename type,int P>
struct GetLagrangeNodes<type,1,P> {
  static void get( std::vector<real_t>& nodes );
};

template<typename type,int P>
struct GetLagrangeNodes<type,2,P> {
  static void get( std::vector<real_t>& nodes );
};

template<typename type,int P>
struct GetLagrangeNodes<type,3,P> {
  static void get( std::vector<real_t>& nodes );
};

template<typename type,int P>
struct GetLagrangeNodes<type,4,P> {
  static void get( std::vector<real_t>& nodes );
};

template<typename type> void get_lagrange_nodes( coord_t number , coord_t order , std::vector<real_t>& nodes );

}
#endif
