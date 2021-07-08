#include "common/error.h"

#include "geometry/entity.h"

#include "graphics/new/tessellation.h"
#include "graphics/gl.h"

#include "mesh/boundary.h"
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

  std::vector<std::map<MeshFacet,index_t>> facets(element.number()+1);
  std::vector<std::vector<MeshFacet>> facets_(element.number()+1);
  std::map<MeshFacet,index_t>::const_iterator it;

  std::vector<index_t> g;
  for (index_t k = 0; k < topology.nb(); k++)
  {
    for (index_t j = 0; j < canonical.size(); j++)
    {
      avro_assert( canonical[j].indices.size() == index_t(canonical[j].dim+1) );

      //if (canonical[j].dim == element.number()) continue;
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

  // convert the facets to primitives
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


static std::vector< std::vector< std::vector<index_t>> > canonical_tet_face = {
  { {} , {} , {} , {} },  // p = 0
  { {2,3,1} , {0,3,2} , {0,1,3} , {0,2,1} }, // p = 1
  { {2,3,1,5,6,4} , {0,3,2,4,7,8} , {0,1,3,5,8,9} , {0,2,1,6,9,7} } , // p = 2
  { {2,3,1,6,7,8,9,4,5,16} , {0,3,2,5,4,10,11,12,13,17} , {0,1,3,7,6,13,12,14,15,18} , {0,2,1,9,8,15,14,11,10,19} } // p = 3
};

static std::vector< std::vector< std::vector<index_t>> > canonical_tri_edge = {
  { {} , {} , {} , {} , {} , {} },  // p = 0
  { {2,1} , {0,2} , {1,0} }, // p = 1
  { {2,1,3} , {0,2,4} , {1,0,5} }, // p = 2
  { {2,1,4,3} , {0,2,6,5} , {1,0,8,7} } // p = 3
};

static std::vector< std::vector< std::vector<index_t>> > canonical_tet_edge = {
  { {} , {} , {} },  // p = 0
  { {2,3} , {1,3} , {1,2} , {0,2} , {0,3} , {0,1} }, // p = 1
  { {2,3,4} , {1,3,5} , {1,2,6} , {0,2,7} , {0,3,8} , {0,1,9} }, // p = 2
  { {2,3,4,5} , {1,3,7,6} , {1,2,8,9} , {0,2,11,10} , {0,3,12,13} , {0,1,14,15} } // p = 3
};


template<typename type>
void
Tessellation::get_primitives( const Topology<type>& topology , const std::vector<std::vector<MeshFacet>>& facets ) {

  // order and number of basis functions of the elements we will create
  coord_t order = topology.element().order();
  index_t nb_basis = nb_simplex_basis(2,order);

  // first add all the points
  points_ = std::make_shared<PointPrimitive>(topology.points());

  // count how many geometry entities there are
  std::set<Entity*> entities;
  for (index_t k = 0; k < topology.points().nb(); k++) {
    if (topology.points().entity(k) == nullptr) continue;
    entities.insert( topology.points().entity(k) );
  }
  printf("there are %lu entities\n",entities.size());
  bool geometryless = (entities.size() == 0);

  index_t nb_Nodes = 0;
  index_t nb_Edges = 0;
  index_t nb_Faces = 0;
  std::map<Entity*,index_t> entity2index;
  for (std::set<Entity*>::const_iterator it = entities.begin(); it != entities.end(); ++it) {
    Entity* entity = *it;
    index_t idx = 0;
    if (entity->number() == 0) idx = nb_Nodes++;
    if (entity->number() == 1) idx = nb_Edges++;
    if (entity->number() == 2) idx = nb_Faces++;
    entity2index.insert( {entity,idx} );
  }

  index_t nb_edge_primitives = nb_Edges + 1;
  index_t nb_triangle_primitives = nb_Faces + 1;

  // if no geometry was given, at least split up the interior/boundary triangles
  coord_t number = topology.number();
  if (nb_triangle_primitives == 1 && number > 2) nb_triangle_primitives = 2;
  if (nb_edge_primitives == 1 && number >= 2) nb_edge_primitives = 2;

  // allocate the primitives
  edges_.resize( nb_edge_primitives );
  for (index_t k = 0; k < nb_edge_primitives; k++)
    edges_[k] = std::make_shared<EdgePrimitive>(order);

  triangles_.resize( nb_triangle_primitives );
  for (index_t k = 0; k < nb_triangle_primitives; k++)
    triangles_[k] = std::make_shared<TrianglePrimitive>(order);

  // loop through triangles and bin them accordingly
  coord_t dim = 2;
  const std::vector<MeshFacet>& facets_d = facets[dim];
  printf("processing %lu %u-facets\n",facets_d.size(),dim);
  for (index_t k = 0; k < facets_d.size(); k++) {

    // retrieve the linear facet
    const MeshFacet& f = facets_d[k];

    // determine if this is a geometry facet
    index_t primitive_id = 0;
    if (number == 2) primitive_id = 0;
    else if (geometryless) {
      if (f.parent.size() == 1) primitive_id = 0; // interior
      else primitive_id = 1; // boundary
    }
    else {
      Entity* entity = BoundaryUtils::geometryFacet( topology.points() , f.indices.data() , f.indices.size() );
      if (entity == nullptr) primitive_id = nb_triangle_primitives -1; // interior at the end
      else primitive_id = entity2index[entity];
    }
    avro_assert( primitive_id < nb_triangle_primitives );

    // retrieve the high-order representation of the triangle
    if (number == 2) {
      // add the dof for the triangle
      index_t elem = f.parent[0];
      triangles_[primitive_id]->add( topology(elem) , topology.nv(elem) );
    }
    else if (number == 3) {
      // pick one of the parents - it doesn't matter which one (unless we want a specific side)
      index_t elem = f.parent[0];
      index_t face = f.local[0];

      // find the dof in the tetrahedron that maps to this triangle face
      std::vector<index_t> dof( nb_basis );
      for (index_t j = 0; j < nb_basis; j++) {
        dof[j] = topology(elem, canonical_tet_face[order][face][j] );
      }
      triangles_[primitive_id]->add( dof.data() , dof.size() );
    }
  }

  // loop through edges and bin them accordingly
  dim = 1;
  const std::vector<MeshFacet>& edges = facets[dim];
  printf("processing %lu %u-facets\n",edges.size(),dim);
  for (index_t k = 0; k < edges.size(); k++) {

    // retrieve the linear facet
    const MeshFacet& f = edges[k];

    // determine if this is a geometry facet
    index_t primitive_id = 0;
    if (number == 2 && geometryless) {
      if (f.parent.size() >= 2) primitive_id = 0; // interior
      else primitive_id = 1; // boundary
    }
    else if (number == 3) {
      // hard to tell if this is a geometry edge
      primitive_id = 0;
    }
    else {
      Entity* entity = BoundaryUtils::geometryFacet( topology.points() , f.indices.data() , f.indices.size() );
      if (entity == nullptr) primitive_id = nb_edge_primitives -1; // interior at the end
      else primitive_id = entity2index[entity];
    }
    avro_assert( primitive_id < nb_edge_primitives );

    // retrieve the high-order representation of the edge
    // pick one of the parents - it doesn't matter which one
    index_t elem = f.parent[0];
    index_t edge = f.local[0];

    // find the dof in the element that maps to this triangle face
    std::vector<index_t> dof( order+1 );
    for (index_t j = 0; j < dof.size(); j++) {
      if (number == 3)
        dof[j] = topology(elem, canonical_tet_edge[order][edge][j] );
      else if (number == 2)
        dof[j] = topology(elem, canonical_tri_edge[order][edge][j] );
      else avro_assert_not_reached;
    }
    edges_[primitive_id]->add( dof.data() , dof.size() );
  }
  edges_[0]->print();


  // generate the vertex array and primitive buffers
  GL_CALL( glGenVertexArrays( 1, &vertex_array_ ) );
  GL_CALL( glBindVertexArray(vertex_array_) );

  points_->write();
  for (index_t k = 0; k < triangles_.size(); k++) {
    triangles_[k]->write();
  }

  for (index_t k = 0; k < edges_.size(); k++) {
    edges_[k]->write();
  }

  const Fields& fields = topology.fields();

  // now create solution primitives to plot on the triangles
  solution_.resize( triangles_.size() );
  for (index_t k = 0; k < solution_.size(); k++) {

    solution_[k] = std::make_shared<FieldPrimitive>();

    printf("need to add %lu field data to field primitive\n",fields.nb());
    // TODO
  }

}

void
Tessellation::draw_triangles( ShaderProgram& shader ) {

  shader.use();

  // bind which attributes we want to draw
  GL_CALL( glBindVertexArray(vertex_array_) );

  GL_CALL( glBindBuffer( GL_ARRAY_BUFFER, points_->buffer()  ) );
  GL_CALL( glVertexAttribPointer( 0, 3, GL_FLOAT, GL_FALSE, 0, 0 ) ); // enable attribute in position 0 which holds coordinates
  GL_CALL( glEnableVertexAttribArray(0) );
  GL_CALL( glBindBuffer( GL_ARRAY_BUFFER , 0 ) );

  // bind the desired colormap
  glActiveTexture(GL_TEXTURE0 + 1);
  GLint colormap_location = glGetUniformLocation(shader.handle() , "colormap");
  glUniform1i(colormap_location, 1); // second sampler in fragment shader

  for (index_t k = 0; k < triangles_.size(); k++) {
    if (k == 1) continue;
    triangles_[k]->draw();
  }
}

void
Tessellation::draw_edges( ShaderProgram& shader ) {

  shader.use();

  // bind which attributes we want to draw
  GL_CALL( glBindVertexArray(vertex_array_) );

  GL_CALL( glBindBuffer( GL_ARRAY_BUFFER, points_->buffer()  ) );
  GL_CALL( glVertexAttribPointer( 0, 3, GL_FLOAT, GL_FALSE, 0, 0 ) ); // enable attribute in position 0 which holds coordinates
  GL_CALL( glEnableVertexAttribArray(0) );
  GL_CALL( glBindBuffer( GL_ARRAY_BUFFER , 0 ) );

  for (index_t k = 0; k < edges_.size(); k++) {
    if (k == 1) continue;
    edges_[k]->draw();
  }
}

} // graphics

} // avro
