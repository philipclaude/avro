#include "common/parallel_for.h"

#include "element/polytope.h"
#include "element/simplex.h"

#include "graphics/new/primitives.h"
#include "graphics/new/vao.h"

#include "mesh/points.h"
#include "mesh/topology.h"

#define TETLIBRARY
#include <tetgen1.5.0/tetgen.h>

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


class PolytopeDecomposition {
public:
  PolytopeDecomposition( coord_t dim , const Points& points , const index_t* v , index_t nv ) :
    dim_(dim),
    points_(points),
    indices_(v),
    nb_vertices_(nv)
  {}

  void prepare() {
    coordinates_.clear();
    for (index_t j = 0; j < nb_vertices_; j++) {
      for (coord_t d = 0; d < dim_; d++)
        coordinates_.push_back( points_[ indices_[j] ][d] );
    }
  }

  void triangulate() {

    // build up the list of coordinates
    prepare();

    // call triangle to triangulate the coordinates
    if (dim_ == 2) {
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
    }
    else if (dim_ == 3) {

      // compute the tetrahedralization of the polyhedron
      tetgenio input,output;

      // fill the vertices
      input.numberofpoints = coordinates_.size() / 3;

      input.pointlist = new REAL[input.numberofpoints*3];
      for (index_t i = 0; i < coordinates_.size(); i++)
        input.pointlist[i] = coordinates_[i];

      std::string switches = "fnnQ";
      tetrahedralize( (char*)switches.c_str() , &input , &output );

      avro_assert( input.numberofpoints == output.numberofpoints );
      for (int i = 0; i < output.numberoftrifaces; i++) {

        int t0 = output.adjtetlist[2*i  ];
        int t1 = output.adjtetlist[2*i+1];

        if (t0 == -1 || t1 == -1) {
          // boundary triangle
          for (coord_t d = 0; d < 3; d++)
            triangles_.push_back( indices_[output.trifacelist[3*i+d]] );
        }

      }
    }
  }

  const std::vector<index_t>& triangles() const { return triangles_; }

private:
  const coord_t dim_;
  const Points& points_;
  std::vector<index_t> triangles_;
  std::vector<real_t> coordinates_;
  const index_t* indices_;
  const index_t nb_vertices_;
};

class Triangulator {
public:
  typedef Triangulator thisclass;

  Triangulator( coord_t dim , const Points& points ) :
    dim_(dim),
    points_(points)
  {}

  void add_polytope( const index_t* v , index_t nv , index_t parent ) {
    polytopes_.push_back( std::make_shared<PolytopeDecomposition>(dim_,points_,v,nv) );
    parents_.push_back(parent);
  }

  void
  triangulate_polytope(index_t k) {
    polytopes_[k]->triangulate();
  }

  void triangulate() {
    ProcessCPU::parallel_for(
      parallel_for_member_callback( this , &thisclass::triangulate_polytope ),
      0,polytopes_.size() );
  }

  index_t nb() const { return polytopes_.size(); }
  index_t parent( index_t k ) const { return parents_[k]; }

  const PolytopeDecomposition& polytope(index_t k) const { return *polytopes_[k].get(); }

private:
  const coord_t dim_;
  const Points& points_;
  std::vector< std::shared_ptr<PolytopeDecomposition> > polytopes_;
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

  Triangulator triangulator( topology.number() , topology.points() );
  for (index_t k = 0; k < topology.nb(); k++) {
    triangulator.add_polytope( topology(k) , topology.nv(k) , k );
  }

  // triangulate the polygons
  triangulator.triangulate();

  // build the primitives
  points_ = std::make_shared<PointPrimitive>(topology.points());

  // edge primitives
  edges_.resize(1);
  edges_[0] = std::make_shared<EdgePrimitive>(order_);
  for (index_t k = 0; k < edges.size()/2; k++)
    edges_[0]->add( &edges[2*k] , 2 );

  // triangle primitives
  triangles_.resize(1);
  triangles_[0] = std::make_shared<TrianglePrimitive>(order_);

  std::vector<index_t> polytopes;
  for (index_t k = 0; k < triangulator.nb(); k++) {

    const PolytopeDecomposition& polytope = triangulator.polytope(k);
    const std::vector<index_t>& triangles = polytope.triangles();
    for (index_t j = 0; j < triangles.size()/3; j++) {
      triangles_[0]->add( triangles.data() + 3*j , 3);
      polytopes.push_back( triangulator.parent(k) );
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
        index_t nb_basis = 1;
        avro_assert( solution_order == 0 ); // for now

        // generate the interpolation data to evaluate the field on the primitive triangles
        Table<real_t> alpha( TableLayout_Rectangular , number );
        std::vector<index_t> parents;

        std::vector<real_t> s(number,1./(number+1));

        for (index_t j = 0; j < triangles_[k]->nb(); j++) {
          // for p = 0, we just use the middle
          alpha.add( s.data() , number );
          parents.push_back( polytopes[j] );
        }

        // evaluate the field at the interpolation points
        for (index_t rank = 0; rank < fld.nb_rank(); rank++) {

          // create field data to hold the solution defined on these triangles
          // this will only by stored on the CPU until the activated field+rank is written to the GPU
          std::shared_ptr<FieldData> data = std::make_shared<FieldData>(solution_order);

          // evaluate the field at the reference coordinates of each element (parent)
          std::vector<real_t> values;
          fld.evaluate( rank, parents, alpha, values );
          avro_assert( values.size() == nb_basis * triangles_[k]->nb() );

          // store the data and save it to the primitive
          for (index_t j = 0; j < triangles_[k]->nb(); j++)
            data->add( values.data() + nb_basis*j , nb_basis );

          // add the field data to the solution primitive
          solution_[k]->add( field_names[i] , rank , data );
          solution_[k]->set_active( field_names[i] , rank );
        }
      }
    }
    avro_assert( solution_.size() == triangles_.size() );

  } // fields.nb > 0


  // set the info for this vao
  std::vector<nlohmann::json> jtriangles;
  for (index_t k = 0; k < triangles_.size(); k++) {
    json jt;
    jt["name"] = "triangles" + std::to_string(k);
    jt["order"] = triangles_[k]->order();
    triangles_[k]->visible() = true;
    jtriangles.push_back(jt);
  }
  info_["triangles"] = jtriangles;

  std::vector<nlohmann::json> jedges;
  for (index_t k = 0; k < edges_.size(); k++) {
    json je;
    je["name"] = "edges" + std::to_string(k);
    je["order"] = edges_[k]->order();
    edges_[k]->visible() = true;
    jedges.push_back(je);
  }
  info_["edges"] = jedges;

  std::vector<nlohmann::json> jfields;
  for (index_t k = 0; k < fields.nb(); k++) {
    const FieldHolder& fld = fields[field_names[k]];
    std::vector<std::string> rank_names;
    std::vector<index_t> rank_index;
    for (index_t i = 0; i < fld.nb_rank(); i++) {
      rank_names.push_back( fld.get_name(i) );
      rank_index.push_back(i);
    }

    json jf;
    jf["name"]  = field_names[k];
    jf["ranks"] = rank_names;
    jf["rank_index"] = rank_index;
    jfields.push_back(jf);
  }
  info_["fields"]      = jfields;
  info_["field_names"] = field_names;

}


} // graphics

} // avro
