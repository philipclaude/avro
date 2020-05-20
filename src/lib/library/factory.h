// avro: Adaptive Voronoi Remesher
// Copyright 2017-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef AVRO_LIBRARY_FACTORY_H_
#define AVRO_LIBRARY_FACTORY_H_

#include "common/types.h"

#include <memory>
#include <stdlib.h>
#include <string>
#include <vector>

namespace avro
{

class Context;
class Model;
class Points;
class Mesh;
template<typename type> class Topology;
class MetricAttachment;

namespace library
{

std::shared_ptr<MetricAttachment>
get_metric( const std::string& name , Points& points , bool& is_analytic ,
           const std::vector<real_t>& params = std::vector<real_t>() );

std::shared_ptr<Model>
get_geometry( const std::string& name , bool& curved );

std::shared_ptr<Mesh>
get_mesh( const std::string& name , std::shared_ptr<TopologyBase>& ptopology , coord_t number=-1 );

} // programs

} // avro

#endif
