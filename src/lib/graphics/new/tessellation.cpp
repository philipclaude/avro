#include "common/error.h"

#include "graphics/new/tessellation.h"

#include "mesh/points.h"
#include "mesh/topology.h"

#include "numerics/linear_algebra.h"
#include "numerics/mat.h"

#include "avro_types.h"

#include <vector>

namespace avro
{

namespace graphics
{

typedef struct
{
  coord_t dim;
  std::vector<index_t> indices;
} Facet;

struct CanonicalFacet : Facet
{
  index_t local;
};

struct MeshFacet : Facet
{
  std::vector<index_t> parent;
  std::vector<index_t> local;
  std::vector<int>     orientation;
};


// needed to create a set/map of Facet
bool
operator<( const Facet& f , const Facet& g )
{
  // first check the topological dimension
  if (f.dim < g.dim)
    return true;

  // lexicographically compare the indices
  return std::lexicographical_compare(f.indices.begin(), f.indices.end(),
                                      g.indices.begin(), g.indices.end());
}

// needed to create a set/map of elements
bool
operator==( const Facet& fx , const Facet& fy )
{
  // assumes fx and fy have the same topological dimension
  // and that the indices are sorted
  avro_assert( fx.dim == fy.dim );
  for (index_t j = 0; j < fx.dim; j++)
    if (fx.indices[j] != fy.indices[j])
      return false;
  return true;
}

void
get_canonical_simplex_facets( coord_t number , std::vector<CanonicalFacet>& facets )
{
	facets.clear();

  if (number == 0) {
    CanonicalFacet f0; f0.dim = 0; f0.indices = {0}; f0.local = 0;
    facets.push_back(f0);
  }
  else if (number == 1) {
    // vertices
    CanonicalFacet f0; f0.dim = 0; f0.indices = {0}; f0.local = 0;
    CanonicalFacet f1; f1.dim = 0; f1.indices = {1}; f1.local = 1;
    CanonicalFacet f1_0; f1_0.dim = 1; f1_0.indices = {0,1}; f1_0.local = 0;
    facets.push_back( f0 );
    facets.push_back( f1 );
    facets.push_back( f1_0 );
  }
  else if (number == 2) {
    // vertices
    CanonicalFacet f0; f0.dim = 0; f0.indices = {0}; f0.local = 0;
    CanonicalFacet f1; f1.dim = 0; f1.indices = {1}; f1.local = 1;
    CanonicalFacet f2; f2.dim = 0; f2.indices = {2}; f2.local = 2;

    // edges (normal facing outwards)
    CanonicalFacet f1_0; f1_0.dim = 1; f1_0.indices = {2,1}; f1_0.local = 0;
    CanonicalFacet f1_1; f1_1.dim = 1; f1_1.indices = {0,2}; f1_1.local = 1;
    CanonicalFacet f1_2; f1_2.dim = 1; f1_2.indices = {1,0}; f1_2.local = 2;

    CanonicalFacet f2_0; f2_0.dim = 2; f2_0.indices = {0,1,2}; f2_0.local = 0;

    facets.push_back( f0 );
    facets.push_back( f1 );
    facets.push_back( f2 );
    facets.push_back( f1_0 );
    facets.push_back( f1_1 );
    facets.push_back( f1_2 );
    facets.push_back( f2_0 );
  }
  else if (number == 3) {
    // vertices
    CanonicalFacet f0; f0.dim = 0; f0.indices = {0}; f0.local = 0;
    CanonicalFacet f1; f1.dim = 0; f1.indices = {1}; f1.local = 1;
    CanonicalFacet f2; f2.dim = 0; f2.indices = {2}; f2.local = 2;
    CanonicalFacet f3; f3.dim = 0; f3.indices = {3}; f3.local = 3;

    // edges (orientation doesn't matter, but this is consistent with SANS -- good for testing!)
    CanonicalFacet f1_0; f1_0.dim = 1; f1_0.indices = {0,1}; f1_0.local = 5;
    CanonicalFacet f1_1; f1_1.dim = 1; f1_1.indices = {1,2}; f1_1.local = 2;
    CanonicalFacet f1_2; f1_2.dim = 1; f1_2.indices = {0,2}; f1_2.local = 3;
    CanonicalFacet f1_3; f1_3.dim = 1; f1_3.indices = {0,3}; f1_3.local = 4;
    CanonicalFacet f1_4; f1_4.dim = 1; f1_4.indices = {1,3}; f1_4.local = 1;
    CanonicalFacet f1_5; f1_5.dim = 1; f1_5.indices = {2,3}; f1_5.local = 0;

    // triangles (normal facing outwards)
    CanonicalFacet f2_0; f2_0.dim = 2; f2_0.indices = {2,3,1}; f2_0.local = 0;
    CanonicalFacet f2_1; f2_1.dim = 2; f2_1.indices = {0,3,2}; f2_1.local = 1;
    CanonicalFacet f2_2; f2_2.dim = 2; f2_2.indices = {0,1,3}; f2_2.local = 2;
    CanonicalFacet f2_3; f2_3.dim = 2; f2_3.indices = {0,2,1}; f2_3.local = 3;

    // tetrahedron
    CanonicalFacet f3_0; f3_0.dim = 3; f3_0.indices = {0,1,2,3}; f3_0.local = 0;

    facets.push_back( f0 );
    facets.push_back( f1 );
    facets.push_back( f2 );
    facets.push_back( f3 );

    facets.push_back( f1_5 );
    facets.push_back( f1_4 );
    facets.push_back( f1_1 );
    facets.push_back( f1_2 );
    facets.push_back( f1_3 );
    facets.push_back( f1_0 );

    facets.push_back( f2_0 );
    facets.push_back( f2_1 );
    facets.push_back( f2_2 );
    facets.push_back( f2_3 );

    facets.push_back( f3_0 );
  }
  else
		avro_assert_not_reached;
}

short
find_orientation( const std::vector<index_t>& f , const std::vector<index_t>& g )
{
  index_t n = f.size();
  avro_assert( g.size() == n);

  // build the permutation matrix and determine sign via determinant
  matd<int> P(n,n); // starts as zeros
  for (index_t j = 0; j < n; j++)
  {
    index_t a = f[j];
    for (index_t i = 0; i < n; i++)
    {
      index_t b = g[i];
      if (a == b) {
        P(j,i) = 1;
      }
    }
  }
  return numerics::det(P);
}

template<>
void
Tessellation::_build( const Topology<Simplex>& topology ) {

  // get the canonical representation of the facets of the element
  const Simplex& element = topology.element();
  std::vector<CanonicalFacet> canonical;
  get_canonical_simplex_facets(element.number(),canonical);

  std::vector<std::map<MeshFacet,index_t>> facets(element.number());
  std::vector<std::vector<MeshFacet>> facets_(element.number());
  std::map<MeshFacet,index_t>::const_iterator it;

  std::vector<index_t> g;
  for (index_t k = 0; k < topology.nb(); k++)
  {
    for (index_t j = 0; j < canonical.size(); j++)
    {
      avro_assert( canonical[j].indices.size() == (canonical[j].dim+1) );

      if (canonical[j].dim == element.number()) continue;
      if (canonical[j].dim > 2) continue;

      // determine the indices of the facet
      MeshFacet f;
      f.dim = canonical[j].dim;
      f.indices.resize( canonical[j].indices.size() , 0 );
      for (index_t i = 0; i < f.indices.size(); i++)
        f.indices[i] = topology[k][ canonical[j].indices[i] ];

      // save the indices, then sort and determine positive/negative orientation
      g = f.indices;
      std::sort( f.indices.begin() , f.indices.end() );
      short orientation = find_orientation(g,f.indices);

      //print_inline(f.indices);

      // determine if this facet exists
      it = facets[f.dim].find(f);
      if (it == facets[f.dim].end()) {

        // this facet doesn't yet exist
        avro_assert_msg( facets_[f.dim].size() == facets[f.dim].size() ,
                          "|facets_| = %lu, |facets| = %lu",facets_.size(),facets.size());
        facets[f.dim].insert( {f,facets_[f.dim].size()} );
        f.orientation.push_back( orientation );
        f.local.push_back( canonical[j].local );
        f.parent.push_back( k );
        facets_[f.dim].push_back( f );
      }
      else {
        // this facet exists, add the parent data
        index_t id = it->second;
        avro_assert( id < facets_[f.dim].size() );
        facets_[f.dim][id].parent.push_back(k);
        facets_[f.dim][id].local.push_back( canonical[j].local );
        facets_[f.dim][id].orientation.push_back(orientation);
      }
    }
  }

  for (coord_t d = 0; d < facets_.size(); d++)
  for (index_t k = 0; k < facets_[d].size(); k++) {
    print_inline( facets_[d][k].indices );
  }

  get_primitives(topology,facets_);
}

void
triangulate( const Points& points , const index_t* v , index_t nv , std::vector<index_t>& triangles ) {
  avro_implement;
}

template<>
void
Tessellation::_build( const Topology<Polytope>& topology ) {

  // only linear polytope meshes are supported
  avro_assert( topology.element().order() == 1 );

  std::vector<index_t> edges;
  topology.get_edges(edges);

  // TODO: loop through the fields and make sure all fields are either zeroth or first order
  Topology<Simplex> simplices( topology.points() , 2 );

  std::vector<index_t> triangles;
  if (topology.number() == 2) {

    for (index_t k = 0; k < topology.nb(); k++) {
      triangles.clear();
      triangulate( topology.points() , topology(k) , topology.nv(k) , triangles );

      index_t nb_triangles = triangles.size()/3;
      for (index_t j = 0; j < nb_triangles; j++)
        simplices.add( triangles.data()+3*j , 3 );
    }
  }
  else if (topology.number() == 3) {

    std::vector<int> hrep;
    std::vector<index_t> vrep;
    for (index_t k = 0; k < topology.nb(); k++) {
      hrep.clear();
      topology.element().hrep( topology(k) , topology.nv(k) , hrep );

      for (index_t i = 0; i < hrep.size(); i++) {
        triangles.clear();
        vrep.clear();
        topology.element().vrep( topology(k) , topology.nv(k) , hrep[i] , vrep );

        avro_assert( vrep.size() > 2 );
        triangulate( topology.points() , vrep.data() , vrep.size() , triangles );

        index_t nb_triangles = triangles.size()/3;
        for (index_t j = 0; j < nb_triangles; j++)
          simplices.add( triangles.data()+3*j , 3 );
      }

    }

  }


}

void
Tessellation::build( const TopologyBase& topology ) {
  if (topology.type_name() == "simplex")
    _build( static_cast<const Topology<Simplex>&>(topology) );
  else
    avro_implement;
}

void
Tessellation::get_primitives( const TopologyBase& topology , const std::vector<std::vector<MeshFacet>>& facets ) {

}

} // graphics

} // avro