// avro: Adaptive Voronoi Remesher
// Copyright 2017-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef AVRO_MESH_DELAUNAY_TRIANGULATION_H_
#define AVRO_MESH_DELAUNAY_TRIANGULATION_H_

#include "common/types.h"

#include "mesh/topology.h"
#include "mesh/types.h"
#include "mesh/vertices.h"

namespace avro
{

namespace delaunay
{

class RestrictedVoronoiDiagram;

class Triangulation : public Topology<Simplex>
{
public:
    // constructor from a set of vertices
    Triangulation( Vertices& v ) :
        Topology<Simplex>( v , v.dim() ) {}

    // constructor from a restricted Voronoi diagram
    Triangulation( const RestrictedVoronoiDiagram& rvd );

    void compute();

    void initialize();
    index_t nextClosest( const index_t p , const std::vector<index_t> skip ) const;

    void iterate();
    void step();

    bool insphere( const index_t n , const index_t p ) const;
    void add( const index_t n );

private:
    std::vector<index_t> facet_;
    Data<index_t> facets_;

};

} // delaunay

} // avro

#endif
