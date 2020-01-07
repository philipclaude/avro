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
  virtual void get_simplices( coord_t number , std::vector<index_t>& simplices , std::vector<index_t>& parents ) const = 0;

  const Table<real_t>& reference_coordinates() const { return reference_coordinates_; }
  const std::vector<index_t>& point_parents() const { return point_parents_; }

  bool native( index_t k ) const { return native_[k]; }

protected:
  Points points_;
  std::vector<index_t> native_;
  Table<real_t> reference_coordinates_;
  std::vector<index_t> point_parents_;
  std::vector< std::map<index_t,std::vector<index_t> > > parents_;
};

template<typename type>
class Triangulation : public TriangulationBase
{
public:
  Triangulation( const Topology<type>& topology );

  void extract();

  index_t add_simplex( index_t number , const index_t* v , index_t parent );

  const Topology<type>& topology() const { return topology_; }

  void get_simplices( coord_t number , std::vector<index_t>& simplices , std::vector<index_t>& parents ) const;

  index_t add_point( coord_t number , const index_t *v , index_t nv , index_t parent );

private:
  const Topology<type>& topology_;

  std::vector< std::map<Element,index_t> > elements_;
  std::map<Element,index_t> centroids_;
  std::map<index_t,coord_t> centroid2dim_;
};

} // avro

#endif
