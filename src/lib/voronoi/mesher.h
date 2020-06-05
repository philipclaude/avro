//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef AVRO_MESH_DELAUNAY_MESHER_H_
#define AVRO_MESH_DELAUNAY_MESHER_H_

#if 0

#include "common/types.h"

#include "geometry/tessellation.h"

#include "mesh/topology.h"
#include "mesh/vertices.h"

namespace avro
{

class Mesh;
class Model;
class Entity;

template <typename type>
class Mesher
{
protected:
  typedef smart_ptr(Topology<type>) Topology_ptr;

  Mesher( Model* _model , Mesh& _mesh );
  Mesher( Topology<type>& boundary , Mesh<type>& _mesh );

  Mesher( Vertices& vertices , Mesh<type>& _mesh ) :
    model_(NULL),
    mesh_(_mesh),
    boundaryFacets_(vertices,mesh_.number()) {}

  Mesher( Mesh<type>& boundary , Mesh<type>& _mesh ) :
    model_(NULL),
    mesh_(_mesh),
    boundaryFacets_(boundaryVertices_,boundary.number()) {}

  Mesher( ModelTessellation& tess , Mesh<type>& _mesh ) :
    model_(NULL),
    mesh_(_mesh),
    boundaryFacets_(tess.vertices(),tess.number()) {}

  void setDimension( const coord_t dim )
  {
    dim_ = dim;
    boundaryVertices_.setDimension(dim_);
    internalVertices_.setDimension(dim_);
  }

  int label( Entity* entity ) const;

public:
  Mesh<type>& mesh() const { return mesh_; }

  void addInternalVertices( const Vertices& vi )
  {
    avro_assert( vi.dim()==dim_ );
    for (index_t k=0;k<vi.nb();k++)
      internalVertices_.create( vi[k] );
  }

protected:
  Model* model_;
  Mesh<type>& mesh_;
  coord_t dim_;
  real expectedVolume_;
  real actualVolume_;

  Topology<type> boundaryFacets_;
  Vertices boundaryVertices_;
  Vertices internalVertices_;

  std::vector<Entity*> entities_;
};

class Triangle : public Mesher<Simplex>
{
public:
  Triangle( Model* _model , Mesh<Simplex>& _mesh );
  Triangle( Topology<Simplex>& boundary , Mesh<Simplex>& _mesh );
  Triangle( Vertices& vertices , Mesh<Simplex>& _mesh );
  Triangle( Mesh<Simplex>& boundary , Mesh<Simplex>& _mesh );
  Triangle( ModelTessellation& tess , Mesh<Simplex>& _mesh );
  void call( const std::string& switches );
};

class TetGen : public Mesher<Simplex>
{
public:
  TetGen( Model* _model , Mesh<Simplex>& _mesh );
  TetGen( Topology<Simplex>& boundary , Mesh<Simplex>& _mesh );
  TetGen( Vertices& vertices , Mesh<Simplex>& _mesh );
  TetGen( Mesh<Simplex>& boundary , Mesh<Simplex>& _mesh );
  void call( const std::string& switches );
};

} // avro

#endif

#endif
