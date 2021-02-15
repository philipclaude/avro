#include "graphics/clipping.h"
#include "graphics/gl.h"
#include "graphics/manager.h"
#include "graphics/math.h"

#include "library/ckf.h"

#include "numerics/geometry.h"
#include "numerics/quaternion.h"

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
  sign_(1),
  distance_(0),
  u_(dim),
  v_(dim)
{

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
  }

  initialize();
}

void
ClippingPlane::set_coordinates( const real_t* bounding_box )
{
  if (bounding_box==nullptr) return;

  // set the z-coordinate of the clipping plane to the middle of the bounding box
  real_t zm = 0.5*( bounding_box[5] + bounding_box[2] );

  // find the maximum length
  real_t lx = bounding_box[3] - bounding_box[0];
  real_t ly = bounding_box[4] - bounding_box[1];

  real_t lmax = lx;
  if (ly > lx) lmax = ly;

  for (index_t k=0;k<points_.nb();k++)
  {
    for (coord_t d=0;d<2;d++)
      points_[k][d] = (points_[k][d] -0.5)*lmax*2 -  bounding_box[d];
    points_[k][2] = zm;
  }
  initialize();
}

void
ClippingPlane::initialize()
{
  transformed_points_.clear();
  transformation_ = mat4(1.0f);

  for (index_t k=0;k<points_.nb();k++)
    transformed_points_.create( points_[k] );

  distance_ = 0.0;
  angles_[0] = 0.0;
  angles_[1] = 0.0;
  u_[0] = 0.;
  u_[1] = 1.;
  u_[2] = 0.;

  v_[0] = 0.;
  v_[1] = 0.;
  v_[2] = 1.;

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

  for (coord_t d=0;d<3;d++)
  {
    u_[d] = u[d];
    v_[d] = v[d];
  }
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
ClippingPlane::plot( GraphicsManager& manager , const real_t* focus , const mat4& transformation )
{
  //OpenGL_Manager* manager = dynamic_cast<OpenGL_Manager*>(&manager0);

  real_t scale = 1./focus[3];
  for (index_t k=0;k<transformed_points_.nb();k++)
  {
    for (coord_t d=0;d<3;d++)
      plot_points_[3*k+d] = scale*(transformed_points_[k][d] - focus[d]);
    plot_colors_[3*k] = 255.0;
  }

  plot_triangles_[0] = 0; plot_triangles_[1] = 1; plot_triangles_[2] = 2;
  plot_triangles_[3] = 1; plot_triangles_[4] = 2; plot_triangles_[5] = 3;

  manager.write( "clip-plane" , 2 , plot_points_ , {} , plot_triangles_ , plot_colors_ );
  DrawingParameters params;
  params.mvp = transformation;
  params.transparency = 0.25;
  //params.lighting = -1;
  manager.draw("clip-plane",2,params);
}

void
ClippingPlane::update( real_t distance , real_t* angles , int dir )
{
  sign_ = dir;
  initialize();

  // transform the points
  for (index_t k=0;k<4;k++)
  {
    for (coord_t d=0;d<3;d++)
      transformed_points_[k][d] -= normal_[d]*distance_;
  }

  mat4 M = mat4(1.0);
  mat4 T = mat4(1.0);

  // compound the rotation into M by first translating to the initial point
  for (coord_t d=0;d<3;d++)
    T[d][3] = -normal_[d]*distance_;
  mat4 M0 = T*M;

  // save the new distance
  distance_ = distance;

  // get the angles in radians
  angles[0] = angles[0]*M_PI/180.;
  angles[1] = angles[1]*M_PI/180.;

  mat4 Mr1 = mat4(0);
  mat4 Mr2 = mat4(0);
  Mr1[3][3] = 1;
  Mr2[3][3] = 1;

  // rotation about the first axis
  Quaternion q1(angles[0]-angles_[0],u_.data());
  mat3 R = q1.rotation_matrix();
  for (coord_t i=0;i<3;i++)
  for (coord_t j=0;j<3;j++)
    Mr1[i][j] = R[i][j];

  M = Mr1*M0;
  M0 = M;

  vec3 v0 = {v_[0],v_[1],v_[2]};
  vec3 v = R*v0;
  for (coord_t d=0;d<3;d++)
    v_[d] = v[d];

  // rotation about the second axis
  Quaternion q2(angles[1]-angles_[1],v_.data());
  R = q2.rotation_matrix();
  for (coord_t i=0;i<3;i++)
  for (coord_t j=0;j<3;j++)
    Mr2[i][j] = R[i][j];
  M = Mr2*M0;

  // also rotate the first axis
  vec3 u0 = {u_[0],u_[1],u_[2]};
  vec3 u = R*u0;
  for (coord_t d=0;d<3;d++)
    u_[d] = u[d];

  // save the angles
  angles_[0] = angles[0];
  angles_[1] = angles[1];

  // then translating back
  for (coord_t d=0;d<3;d++)
    T[d][3] = normal_[d]*distance_;
  M0 = T*M;

  // compute the normal as the cross product of the two axis vectors
  normal_[0] = u_[1]*v_[2] -u_[2]*v_[1];
  normal_[1] = u_[2]*v_[0] -u_[0]*v_[2];
  normal_[2] = v_[1]*u_[0] -v_[0]*u_[1];
  numerics::normalize(normal_.data(),normal_.size());

  // transform the points
  std::fill( center_.begin() , center_.end() , 0.0 );
  for (index_t k=0;k<4;k++)
  {
    vec4 p;
    for (coord_t d=0;d<3;d++)
      p[d] = transformed_points_[k][d];
    p[3] = 1.;
    vec4 x = M0*p;
    for (coord_t d=0;d<3;d++)
    {
      transformed_points_[k][d] = x[d] +normal_[d]*distance_;
      center_[d] += transformed_points_[k][d];
    }
  }

  for (coord_t d=0;d<3;d++)
    center_[d] /= transformed_points_.nb();

}

void
ClippingPlane::hide( GraphicsManager& manager )
{
  manager.remove("clip-plane");
}

void
ClippingPlane::append_transformation( const mat4& mr , const mat4& mt )
{
  // append the rigid-body transformation: rotation (mr) and translation (mt) only
  transformation_ = mr*mt*transformation_;
}

} // graphics

} // avro
