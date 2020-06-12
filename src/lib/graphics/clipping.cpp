#include "graphics/clipping.h"
#include "graphics/gl.h"
#include "graphics/math.h"

#include "library/ckf.h"

#include "numerics/geometry.h"

namespace avro
{

namespace graphics
{

ClippingPlane::ClippingPlane( coord_t dim , real_t length ) :
  normal_(dim),
  center_(dim),
  points_(dim),
  transformed_points_(dim),
  length_(length),
  sign_(1)
{
  initialize();
}

void
ClippingPlane::initialize()
{
  points_.clear(),
  transformed_points_.clear();
  transformation_ = mat4(1.0f);

  // the plane has topological dimension dim_-1
  index_t np = (index_t) std::pow( 2 , points_.dim()-1 );

  std::vector<index_t> dims(points_.dim()-1,2);
  CKF_Triangulation topology(dims);

  plot_points_.resize( 3*topology.points().nb() , 0.0 );
  plot_triangles_.resize( 3*topology.nb() , 0.0 );
  plot_colors_.resize( plot_points_.size() , 0.0 );

  avro_assert( topology.points().nb() == np );
  for (index_t k=0;k<topology.points().nb();k++)
  {
    std::vector<real_t> x( topology.points()[k] , topology.points()[k] + topology.points().dim() );
    x.push_back(0.0);
    points_.create( x.data() );
    for (coord_t d=0;d<points_.dim()-1;d++)
      points_[k][d] = length_*( points_[k][d] - 0.5 ); // ckf triangulation centered on 0.5, we want 0

    transformed_points_.create( points_[k] );
  }

  update();
}

void
ClippingPlane::update()
{
  // this function is specialized for clipping planes in 3d
  std::fill( center_.begin() , center_.end() , 0.0 );

  for (index_t k=0;k<points_.nb();k++)
  {
    vec4 p = { points_[k][0] , points_[k][1] , points_[k][2] , 1. };
    vec4 q = transformation_*p;

    for (index_t d=0;d<points_.dim();d++)
    {
      transformed_points_[k][d] = q[d];
      center_[d] += q[d];
    }
  }

  for (index_t d=0;d<points_.dim();d++)
    center_[d] /= points_.nb();

  // compute the normal and center
  real_t* x0 = transformed_points_[0];
  real_t* x1 = transformed_points_[1];
  real_t* x2 = transformed_points_[2];

  real_t u[3] = { x1[0] - x0[0] , x1[1] - x0[1] , x1[2] - x0[2] };
  real_t v[3] = { x2[0] - x0[0] , x2[1] - x0[1] , x2[2] - x0[2] };

  normal_[0] =   u[1]*v[2] - u[2]*v[1];
  normal_[1] = -(u[0]*v[2] - u[2]*v[0]);
  normal_[2] =   u[0]*v[1] - u[1]*v[0];

  numerics::normalize( normal_.data() , normal_.size() );

  //printf("--> clipping plane: normal = (%g,%g,%g) with center (%g,%g,%g)\n",normal_[0],normal_[1],normal_[2],center_[0],center_[1],center_[2]);
}

bool
ClippingPlane::visible( const Points& points , const index_t* v , index_t nv ) const
{
  for (index_t j=0;j<nv;j++)
  {
    // determine whether this point is on the correct side of the plane
    const real_t* x = points[v[j]];

    // compute the vector between this point and the plane center
    real_t dp = 0.0;
    for (coord_t d=0;d<transformed_points_.dim();d++)
      dp += (x[d] - center_[d])*normal_[d]*sign_;

    if (dp > 0.0) return false;
  }
  return true;
}

void
ClippingPlane::plot( GraphicsManager& manager0 , const real_t* focus , const mat4& transformation )
{
  OpenGL_Manager* manager = dynamic_cast<OpenGL_Manager*>(&manager0);

  real_t scale = 1./focus[3];
  for (index_t k=0;k<transformed_points_.nb();k++)
  {
    for (coord_t d=0;d<3;d++)
      plot_points_[3*k+d] = scale*(transformed_points_[k][d] - focus[d]);
    plot_colors_[3*k] = 255.0;
  }

  plot_triangles_[0] = 0; plot_triangles_[1] = 1; plot_triangles_[2] = 2;
  plot_triangles_[3] = 1; plot_triangles_[4] = 2; plot_triangles_[5] = 3;

  manager->write( "clip-plane" , 2 , plot_points_ , {} , plot_triangles_ , plot_colors_ );
  manager->select_shader( "clip-plane" , "wv" );
  DrawingParameters params;
  params.mvp = transformation;
  params.transparency = 0.25;
  manager->draw("clip-plane",2,params);
}

void
ClippingPlane::append_transformation( const mat4& mr , const mat4& mt )
{
  // append the rigid-body transformation: rotation (mr) and translation (mt) only
  transformation_ = mr*mt*transformation_;
}

} // graphics

} // avro
