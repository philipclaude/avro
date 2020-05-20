#include "geometry/tessellation.h"
#include "geometry/body.h"
#include "geometry/entity.h"
#include "geometry/model.h"

#include "numerics/geometry.h"

#include <egads.h>

namespace avro
{

ModelTessellation::ModelTessellation( Model& model , TessellationParameters& params) :
  Mesh(0,0), model_(model)
{

  // need to set the number of the mesh (it's initialized to 0)
  for (index_t k=0;k<model.nb_bodies();k++)
  {
    if (model.body(k).number()>number_)
      number_ = model.body(k).number();
  }

  // set the dimension of the internal vertices to be filled
  set_number( number_ );
  points_.set_dim( number_+1 );
  points_.set_parameter_dim( number_ );
  internal_points_.set_dim( points_.dim() );
  internal_points_.set_parameter_dim( points_.udim() );

  volume_ = 0.;
  for (index_t k=0;k<model.nb_bodies();k++)
  {

    // save the original number of vertices
    index_t nb0 = points_.nb();

    // tessellate the body in the model
    std::shared_ptr<BodyTessellation> body_tess;
    if (params.type() == "simplex")
    {
      body_tess = std::make_shared<BodyTessellation>( points_ , model.body(k) , params );
    }
    else
      avro_implement;

    /*
    int sign = ( model.is_interior(k) ) ? -1 : 1;
    real db = body_tess.volume();
    //printf("body %d, sign = %d, volume = %g\n",int(k),sign,db);
    volume_ = volume_ +sign*db;
    */

    // copy the body tessellation into the model's tessellation
    // add the root topology of the body tessellation
    avro_assert( body_tess->nb_topologies()==1 );
    add( body_tess->topology_ptr(0) );

    // the vertex needs to receive information as to which body it's on
    // which is the number of children of the root topology
    for (index_t j=nb0;j<points_.nb();j++)
      points_.body(j) = k +1;

    // copy the interior/exterior flag
    //interior_.push_back( model.interior(k) );

    // retrieve the internal points from the body tessellation
    /*
    if (model.is_interior(k))
    {
      body_tess->make_internal_points();
      get_body_internal_points(bodyTess);
    }
    */
  }
}


Model&
ModelTessellation::model() const
{
  return model_;
}

void
ModelTessellation::get_body_internal_points( const BodyTessellation& body_tess )
{
  for (index_t k=0;k<body_tess.nb_internal();k++)
    internal_points_.create( body_tess.internal_point(k) );
  printf("retrieved internal points\n");
}

BodyTessellation::BodyTessellation( Points& model_points , Body& body , TessellationParameters& params ) :
  Mesh(body.number(),model_points.dim()),
  body_(body),
  params_(params), 
  model_points_(model_points)
{
  printf("number = %u, dim = %u\n",body.number(),points_.dim());
  body_.tessellate(*this);
}

#if 0
void
BodyTessellation::makeInternalPoints()
{
  // only closed bodies: solid bodies, closed shells (2-topology) or wirebodies (1-topology)
  if (!body_.closed())
  {
    printf("body is not closed: type = %s\n",body_.memberTypeName().c_str());
    return;
  }
  else printf("body is closed, making internal points\n");

  internalVertices_.setDimension( vertices_.dim() );

  // get all the n-simplices in the lower dimensional topologies
  // the assumption here is that these n-simplices bound the body
  Topology<Simplex> boundary(vertices_,number_);
  retrieveElements( number_ , boundary );

  Mesh<Simplex> volume( vertices_.dim() , number_+1 );
  if (number_==2)
  {
    TetGen tetgen(boundary,volume);
    tetgen.call("pQ");
  }
  else if (number_==1)
  {
    Triangle triangle(boundary,volume);
    triangle.call("pzYQ");
  }
  else
    avro_assert_not_reached;
  printf("called mesher\n");

  Topology<Simplex> simplices( volume.vertices() , number_+1 );
  volume.retrieveElements( number_+1 , simplices );

  // add the centroid of each volume simplex to the internal points list
  std::vector<real> centroid( vertices_.dim() );
  for (index_t k=0;k<simplices.nb();k++)
  {
    geometrics::centroid( simplices(k) , simplices.nv(k) , simplices.vertices() , centroid );
    internalVertices_.create( centroid );
  }
}

real
BodyTessellation::volume()
{
  Topology<Simplex> boundary(vertices_,number_);
  boundary.setSorted(false);
  retrieveElements(number_,boundary);
  return fabs(boundary.calculateBoundingVolume());
}
#endif


} // avro
