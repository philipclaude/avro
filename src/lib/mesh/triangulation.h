#ifndef avro_LIB_MESH_TRIANGULATION_H_
#define avro_LIB_MESH_TRIANGULATION_H_

#include "common/types.h"

#include "master/master.h"
#include "master/simplex.h"
#include "mesh/topology.h"

#include <map>

namespace avro
{

class TriangulationBase : public Topology<Simplex>
{
public:
  TriangulationBase( coord_t number , coord_t dim ) :
    Topology<Simplex>(points_,number,1),
    points_(dim)
  {}

  virtual void extract() = 0;
  virtual void get_triangles( std::vector<index_t>& triangles ) const = 0;

protected:
  Points points_;
};

template<typename type>
class Triangulation : public TriangulationBase
{
public:
  Triangulation( const Topology<type>& topology );

  void extract();

  index_t add_simplex( index_t number , const index_t* v );

  const Topology<type>& topology() const { return topology_; }

  void get_triangles( std::vector<index_t>& triangles ) const;

  index_t add_point( coord_t number , const index_t *v , index_t nv );

private:
  const Topology<type>& topology_;
  std::vector<index_t> parents_;

  std::vector< std::map<Element,index_t> > elements_;
  std::map<Element,index_t> centroids_;
  std::map<index_t,coord_t> centroid2dim_;
};

} // avro

#endif
