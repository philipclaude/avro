#include "element/polytope.h"
#include "element/simplex.h"

#include "graphics/new/primitives.h"
#include "graphics/new/vao.h"

#include "mesh/points.h"
#include "mesh/topology.h"

extern "C" {
#define REAL avro::real_t
#define VOID void
#define ANSI_DECLARATORS
#include <triangle/triangle.h>
}


namespace avro
{

namespace graphics
{

static void
init_triangulateio( struct triangulateio& io ) {

  io.pointlist = NULL;
  io.pointattributelist = NULL;
  io.pointmarkerlist = NULL;
  io.numberofpoints = 0;
  io.numberofpointattributes = 0;

  io.trianglelist = NULL;
  io.triangleattributelist = NULL;
  io.trianglearealist = NULL;
  io.neighborlist = NULL;
  io.numberoftriangles = 0;
  io.numberofcorners = 0;
  io.numberoftriangleattributes = 0;

  io.segmentlist = NULL;
  io.segmentmarkerlist = NULL;
  io.numberofsegments = 0;

  io.holelist = NULL;
  io.numberofholes = 0;

  io.regionlist = NULL;
  io.numberofregions = 0;

  io.edgelist = NULL;
  io.edgemarkerlist = NULL;
  io.normlist = NULL;
  io.numberofedges = 0;
}


class Polygon {
public:
  Polygon( const Points& points , const index_t* v , index_t nv ) :
    points_(points),
    indices_(v),
    nb_vertices_(nv)
  {}

  void prepare() {
    coordinates_.clear();
    if (points_.dim() == 2) {
      for (index_t j = 0; j < nb_vertices_; j++) {
        for (coord_t d = 0; d < points_.dim(); d++)
          coordinates_.push_back( points_[ indices_[j] ][d] );
      }
    }
    else if (points_.dim() == 3) {

      // we need to find the coordinates of the points on the two-dimensional basis of the plane

      avro_implement;

    }
  }

  void triangulate_polygon() {

    // build up the list of coordinates
    prepare();

    // call triangle to triangulate the coordinates
    struct triangulateio input,output;
    init_triangulateio(input);
    init_triangulateio(output);

    input.numberofpoints = nb_vertices_;
    input.pointlist = coordinates_.data();

    ::triangulate("zQN", &input, &output, NULL);

    avro_assert( input.numberofpoints == output.numberofpoints );
    for (int k = 0; k < output.numberoftriangles; k++) {
      for (index_t j = 0; j < 3; j++)
        triangles_.push_back( indices_[output.trianglelist[3*k+j]] );
    }
    //print_inline(triangles_);
  }

  const std::vector<index_t>& triangles() const { return triangles_; }

private:
  const Points& points_;
  std::vector<index_t> triangles_;
  std::vector<real_t> coordinates_;
  const index_t* indices_;
  const index_t nb_vertices_;
};

class Triangulator {
public:
  Triangulator( const Points& points ) :
    points_(points)
  {}

  void add_polygon( const index_t* v , index_t nv , index_t parent ) {
    polygons_.push_back( std::make_shared<Polygon>(points_,v,nv) );
    parents_.push_back(parent);
  }

  void triangulate() {
    for (index_t k = 0; k < polygons_.size(); k++)
      polygons_[k]->triangulate_polygon();
  }

  index_t nb() const { return polygons_.size(); }
  index_t parent( index_t k ) const { return parents_[k]; }

  const Polygon& polygon(index_t k) const { return *polygons_[k].get(); }

private:
  const Points& points_;
  std::vector< std::shared_ptr<Polygon> > polygons_;
  std::vector<index_t> parents_;
};

template<>
void
VertexAttributeObject::_build( const Topology<Polytope>& topology ) {

  number_ = topology.number();
  order_  = topology.element().order();

  // only linear polytope meshes are supported
  avro_assert( topology.element().order() == 1 );

  std::vector<index_t> edges;
  topology.get_edges(edges);

  Triangulator triangulator( topology.points() );
  if (topology.number() == 2) {

    for (index_t k = 0; k < topology.nb(); k++) {
      triangulator.add_polygon( topology(k) , topology.nv(k) , k );
    }
  }
  else if (topology.number() == 3) {

    std::vector<int> hrep;
    std::vector<index_t> vrep;
    for (index_t k = 0; k < topology.nb(); k++) {

      // build up a list of all polytope facets

      // find the vertices on each facet (from one parent)

      // add the vertices to the triangulator

    }
  }

  // triangulate the polygons
  triangulator.triangulate();

  // retrieve the triangles

  // build up the list of MeshFacet objects
  std::vector<std::vector<MeshFacet>> facets(topology.number()+1);

  points_ = std::make_shared<PointPrimitive>(topology.points());

  triangles_.resize(1);
  triangles_[0] = std::make_shared<TrianglePrimitive>(order_);

  for (index_t k = 0; k < triangulator.nb(); k++) {

    const Polygon& polygon = triangulator.polygon(k);
    const std::vector<index_t>& triangles = polygon.triangles();
    for (index_t j = 0; j < triangles.size()/3; j++) {
      MeshFacet facet;

      facet.dim = 2;
      for (index_t i = 0; i < 3; i++)
        facet.indices.push_back(triangles[3*j+i]);

      facet.parent.push_back(k);
      facet.local.push_back(0);
      facet.orientation.push_back(0);

      facets[2].push_back(facet);

      triangles_[0]->add( triangles.data() + 3*j , 3);
    }
  }

  vec3 color = {0.8,0.8,0.2};
  triangles_[0]->set_color(color);


  // now create solution primitives to plot on the triangles
  coord_t number = number_;
  const Fields& fields = topology.fields();
  std::vector<std::string> field_names;
  fields.get_field_names(field_names);
  avro_assert_msg( field_names.size() == fields.nb() , "|names| = %lu, nb_fields = %lu" , field_names.size() , fields.nb() );

  if (fields.nb() > 0) {
    solution_.resize( triangles_.size() );
    for (index_t k = 0; k < solution_.size(); k++) {

      solution_[k] = std::make_shared<FieldPrimitive>();

      printf("need to add %lu field data to field primitive\n",fields.nb());
      for (index_t i = 0; i < fields.nb(); i++) {

        // retrieve the field
        const FieldHolder& fld = fields[field_names[i]];

        // determine the order of the solution
        index_t solution_order = fld.order();
        Simplex simplex(2,solution_order);
        avro_assert( solution_order == 0 ); // for now
        avro_assert( simplex.nb_basis() == 1 );

        // generate the interpolation data to evaluate the field on the primitive triangles
        Table<real_t> alpha( TableLayout_Rectangular , number );
        std::vector<index_t> parents;
        for (index_t j = 0; j < triangles_[0]->nb(); j++) {

          // evaluate the solution at the lagrange nodes of the reference simplex
          for (index_t n = 0; n < simplex.nb_basis(); n++) {

            // determine the reference coordinate in the parent element
            // using the face and orientation
            const real_t* x = simplex.reference().get_reference_coordinate(n);

            // store the interpolation data
            if (number == 2) {
              alpha.add( x , number );
            }
            else if (number == 3) {
              // determine the reference coordinates in the tetrahedron
              avro_implement;
            }
            parents.push_back( j );
          }
        }

        //alpha.print();

        // evaluate the field at the interpolation points
        for (index_t rank = 0; rank < fld.nb_rank(); rank++) {

          // create field data to hold the solution defined on these triangles
          // this will only by stored on the CPU until the activated field+rank is written to the GPU
          std::shared_ptr<FieldData> data = std::make_shared<FieldData>(solution_order);

          // evaluate the field at the reference coordinates of each element (parent)
          std::vector<real_t> values;
          fld.evaluate( rank, parents, alpha, values );
          avro_assert( values.size() == simplex.nb_basis() * triangles_[k]->nb() );

          // store the data and save it to the primitive
          for (index_t j = 0; j < triangles_[k]->nb(); j++)
            data->add( values.data() + simplex.nb_basis()*j , simplex.nb_basis() );

          // add the field data to the solution primitive
          solution_[k]->add( field_names[i] , rank , data );
          solution_[k]->set_active( field_names[i] , rank );

          //print_inline( data->data() );
        }
      }
    }
    avro_assert( solution_.size() == triangles_.size() );
  } // fields.nb > 0

}


} // graphics

} // avro
