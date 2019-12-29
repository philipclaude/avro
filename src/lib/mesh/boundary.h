#ifndef luma_MESH_BOUNDARY_H_
#define luma_MESH_BOUNDARY_H_

#include "common/types.h"

#include "mesh/topology.h"
#include "mesh/points.h"

#include "numerics/geometry.h"

#include <map>

namespace luma
{

namespace BoundaryUtils
{
  Entity*
  geometryFacet( const Points& points , index_t* v , index_t nv , bool elem=false );
}

class Entity;

class BoundaryPoints : public Points
{
public:
  BoundaryPoints( const Points& points , bool meta = true );

  index_t index( index_t k ) const { return local2global_.at(k); }

private:
  std::map<index_t,index_t> local2global_;
};

template<typename type>
class Boundary : public Topology<type>
{
public:
  Boundary( Topology<type>& _topology ) :
      Topology<type>( _topology.points() , _topology.number() ),
      topology_(_topology)
  {}

  void extractall();
  void extract(bool interior=false);
  void extractInterior( Entity* e , Topology<type>& T );

  real_t volume( const index_t k )
  {
    real_t vol = 0.;
    Topology<type>& t = this->child(k);
    for (index_t j=0;j<t.nb();j++)
    {
      std::vector<const real_t*> xk;
      t.get_elem(j,xk);
      vol += numerics::volume( xk , t.points().dim() );
    }
    return vol;
  }

  index_t indexof( Entity* e );

  // returns the entity attached to the k'th boundary topology
  Entity* entity( const index_t k );
  const std::vector<Entity*>& entities() const { return entity_; }

  bool check() const;

  void print() const;

private:
  Topology<type>& topology_;

  std::map< Entity* , index_t > entity2child_;
  std::map< index_t , index_t > bnd2child_;

  std::vector<Entity*> entity_;

};


} // luma

#endif
