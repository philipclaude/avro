//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "geometry/egads/context.h"
#include "library/egads.h"
#include "numerics/geometry.h"

#include <array>

extern "C" int EG_isPlanar(ego object);

#define CROSS(a,b,c)      a[0] = (b[1]*c[2]) - (b[2]*c[1]);\
                          a[1] = (b[2]*c[0]) - (b[0]*c[2]);\
                          a[2] = (b[0]*c[1]) - (b[1]*c[0])
#define DOT(a,b)         (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])

namespace avro
{

namespace EGADS
{

Cube::Cube( const Context* context , coord_t number ) :
  Body(*context)
{
  set_object(object_);
}

Cube::Cube( const Context* context , const std::vector<real_t>& lens , const real_t* x0 ) :
  Cube(context,lens.size())
{
  std::vector<real_t> x(3,0);
  if (x0!=nullptr)
    x.assign(x0,x0+3);

  if (lens.size()==2)
  {
    real_t lx = lens[0];
    real_t ly = lens[1];

    real_t x1[3] = { x[0]     , x[1]     , 0.0 };
    real_t x2[3] = { x[0] +lx , x[1]     , 0.0 };
    real_t x3[3] = { x[0] +lx , x[1] +ly , 0.0 };
    real_t x4[3] = { x[0]     , x[1] +ly , 0.0 };

    Node node0(context,x1);
    Node node1(context,x2);
    Node node2(context,x3);
    Node node3(context,x4);

    Edge edge0(context,node0,node1);
    Edge edge1(context,node1,node2);
    Edge edge2(context,node2,node3);
    Edge edge3(context,node3,node0);

    EdgeLoop loop(context,CLOSED);
    loop.add(edge0,1);
    loop.add(edge1,1);
    loop.add(edge2,1);
    loop.add(edge3,1);
    loop.make();

    WireBody wire(context,loop);

    object_ = wire.object();
    set_object( wire.object() );
  }
  else if (lens.size()==3)
  {
    #ifndef AVRO_NO_ESP
    real_t data[6] = { x[0] , x[1] , x[2] , lens[0] , lens[1] , lens[2] };
    EGADS_ENSURE_SUCCESS( EG_makeSolidBody( *context_.get() , BOX , data , &object_ ) );
    #else
    printf("need full egads to make solid body\n");
    avro_assert_not_reached;
    #endif
    set_object( object_ );
  }
  this->build_hierarchy();
}

#ifndef AVRO_NO_ESP
void
SolidBody::make()
{
  EGADS_CHECK_SUCCESS( EG_makeSolidBody( *context_.get() , type_ , data_ , &object_ ) );
  set_object(object_);
  build_hierarchy();
}

Sphere::Sphere( const Context* context , real_t* x0 , real_t r ) :
  SolidBody(context,SPHERE)
{
  add_data(x0,3);
  add_data(&r,1);
  make();
}

Cone::Cone( const Context* context , real_t* apex , real_t *x0 , real_t r ) :
  SolidBody(context,CONE)
{
  add_data(apex,3);
  add_data(x0,3);
  add_data(&r,1);
  make();
}

Cylinder::Cylinder( const Context* context , real_t *x0 , real_t* x1 , real_t r ) :
  SolidBody(context,CYLINDER)
{
  add_data(x0,3);
  add_data(x1,3);
  add_data(&r,1);
  make();
}

Torus::Torus( const Context* context , real_t* x0 , real_t* dir , real_t R , real_t r ) :
  SolidBody(context,TORUS)
{
  add_data(x0,3);
  add_data(dir,3);
  add_data(&R,1);
  add_data(&r,1);
  make();
}

Plane::Plane( const Context* context , std::vector<real_t*>& x , real_t* uv ) :
  Object(*context)
{
  avro_assert( x.size()==3 );
  real_t u[3],v[3];
  for (coord_t d=0;d<3;d++)
  {
    u[d] = x[1][d] -x[0][d];
    v[d] = x[2][d] -x[0][d];
  }

  // get the normal to the simplex
  real_t n[3];
  CROSS(n,u,v);

  // get u x v
  CROSS(v,u,n);

  numerics::normalize(u,3);
  numerics::normalize(v,3);

  double data[9] = { x[0][0],  x[0][1],  x[0][2],
                    u[0], u[1], u[2],
                    v[0], v[1], v[2] };
  EGADS_CHECK_SUCCESS( EG_makeGeometry( *context->get() , SURFACE , PLANE , NULL , NULL , data , &object_ ) );

  // TODO use a better range
  if (uv==NULL)
  {
    range_[0] = 0;
    range_[1] = 1;
    range_[2] = 0;
    range_[3] = 1;
  }
  else
  {
    for (index_t j=0;j<4;j++)
      range_[j] = uv[j];
  }
}

// constructor from a set of points, end condition and tolerance
Spline::Spline( const Context* _context , const EndConditions end , const real_t tol , const std::vector<real_t>& points , bool _surface ) :
  Curve(_context),
  context_(_context), surface_(_surface), end_(end), tol_(tol)
{
  int sizes[2] = { (int)points.size()/3 , surface_ ? 1 : 0 };

  int status = EG_approximate( *context_->get() , (int)end , tol , sizes , &points[0] , &object_ );
  if (status==EGADS_SUCCESS) ok_ = true;
  else ok_ = false;
  EGADS_CHECK_SUCCESS( status ); // for printing message;
}

Spline::EndConditions
Spline::condition() const
{
  return end_;
}

real_t
Spline::tolerance() const
{
  return tol_;
}

const Context*
Spline::context() const
{
  return context_;
}

#define NEDGE 4

Face3D::Face3D( const Context* context , int style ) :
  Body(*context)
{
  // construct the nodes
  real_t xyz_[3] = {0.,0.,0.};
  Node n0(context,xyz_);

  xyz_[0] = 1.;
  Node n1(context,xyz_);

  xyz_[0] = 0.;
  xyz_[1] = 2.;
  Node n2(context,xyz_);

  xyz_[0] = 1.;
  Node n3(context,xyz_);

  // construct the curves
  double p[3] = {0.,0.,0.};
  double dx[3] = {1.,0.,2.};
  double dy[3];
  real_t r;

  p[0] = .5;
  dx[0] = dy[2] = -1.;
  p[1] = p[2] = dx[1] = dx[2] = dy[0] = dy[1] = 0.;
  r = 0.5;
  Circle circle0( context , p , dx , dy , r );

  p[0] = p[1] = p[2] = dx[2] = 0.;
  dx[0] = 0.;
  dx[1] = 2.;
  Line line1( context , p , dx );

  p[0]  = 1.;
  #if NEDGE == 4
  dx[0] = -1.;
  #endif
  Line line2( context , p , dx );

  p[0] = 0.5;
  p[1] = 2.;
  p[2] = 0.;

  dx[0] = dy[2] = -1.;
  dx[1] = dx[2] = dy[0] = dy[1] = 0.;
  r = .5;
  Circle circle3( context , p , dx , dy , r );

  // construct the edges
  Edge edge0(context,circle0,n0,n1);
  Edge edge1(context,n0,n2);

  #if NEDGE == 3
  Edge edge2(context,n1,n2);
  #else
  Edge edge2(context,circle3,n3,n2);
  Edge edge3(context,n1,n3);
  #endif

  // make a surface without a loop
  EdgeLoop loop0(context,CLOSED);
  loop0.add( edge0 , -1 );
  loop0.add( edge1 , 1 );
  loop0.add( edge2 , -1 );
  #if NEDGE==4
  loop0.add( edge3 , -1 );
  #endif


  // make the loop
  loop0.make();

  // make the face
  std::shared_ptr<Face> face = nullptr;
  if (EG_isPlanar(*loop0.object()) == 0)
  {
    face = std::make_shared<Face>(loop0,SFORWARD);
  }
  else
  {
    // get the surface that fits the loop
    Isocline surface0(loop0,style);
    loop0.delete_object();

    for (index_t k=0;k<NEDGE;k++)
      loop0.add_other( surface0 , k );

    // update the loop
    loop0.make( *surface0.object() );
    surface0.set_object( surface0.object() );

    face = std::make_shared<Face>( context , surface0 , loop0 , SFORWARD );
  }

  // make a face body
  FaceBody face_body(context,*face.get());
  set_object( face_body.object() );
  object_ = face_body.object();

  // retrieve the info about object class and model type
  //EG_getInfo( object_ , &data_.object_class , &data_.member_type , &data_.reference , &data_.previous , &data_.next );

  //number_ = 2;
  build_hierarchy();

  print();
}

NACA_Airfoil::NACA_Airfoil( const Context* context , real_t chord , real_t* digits0 , index_t np  ) :
  Body(*context)
{
  real_t digits[4] = {0,0,1,2};
  te_[0] = chord;
  te_[1] = te_[2] = 0.;
  if (digits0!=NULL)
  {
    for (index_t j=0;j<4;j++)
      digits[j] = digits0[j];
  }

  printf("get naca airfoil points!\n");
  avro_implement;

  // interpolate the airfoil
  std::vector<real_t> points;
  Spline spline( context , Spline::Slope , 1e-7 , points , false );
  Edge edge(context,spline);

  EdgeLoop loop(context,CLOSED);
  loop.add(edge,1);
  loop.make();

  WireBody body(context,loop);
  set_object( body.object() );

  build_hierarchy();
}

void
NACA_Airfoil::add_wake( const Context* context , const real_t length , real_t* dir0 )
{
  real_t dir[3] = {0,0,0};
  if (dir0==NULL)
  {
    dir[0] = 1.;
    dir[1] = 0.;
    dir[2] = 0.;
  }

  // get the node representing the trailing edge
  ego te = get_trailing_edge();

  real_t xyzWake[3];
  for (coord_t d=0;d<3;d++)
    xyzWake[d] = te_[d] +length*dir[d];

  // create the wake edge
  Node teNode(context,&te);
  Node wakeNode( context , xyzWake );
  Edge wake(context,teNode,wakeNode);

  // retrieve the airfoil edge
  //                   loop       edge
  ego* airfoil_ego = static_cast<EGADS::Object*>(&this->child(0)->child(0))->object();
  Edge airfoil(context,airfoil_ego);

  // create an open loop with both airfoil and wake edges
  EdgeLoop loop(context,OPEN);
  loop.add( airfoil , -1 );
  loop.add( wake , -1 );
  loop.make();

  // make a wirebody from the single loop
  WireBody body(context,loop);
  set_object(body.object());

  // build the entity hierarchy
  build_hierarchy();
}

ego
NACA_Airfoil::get_trailing_edge()
{
  //     loop         edge     node
  return *static_cast<EGADS::Object*>(&this->child(0)->child(0).child(0))->object();
}

/*
Circle::Circle( const Context* context , const real_t* xc , const real_t R )
{
  real_t dx[3] = {1.,0.,0.};
  real_t dy[3] = {0.,1.,0.};
  Circlele circle( context , xc , dx , dy , R  );
  Edge edge(context,circle);

  te_[0] = R;
  te_[1] = te_[2] = 0.;

  EdgeLoop loop(CLOSED);
  loop.add(edge,1);
  loop.make(context);

  WireBody wire(context,loop);
  set_object(wire.object());

  build_hierarchy();
}
*/

Smiley::Smiley( const Context* context , real_t* x0 , real_t rf , real_t rm , real_t hm , real_t tm , real_t re , real_t de , real_t te ) :
  Body(*context)
{
  Circle head(context,x0,rf);

  real_t x1[3] = { x0[0] +  rm    *cos(tm) , x0[1] -  rm    *sin(tm) , 0.0 };
  real_t x2[3] = { x0[0] + (rm+hm)*cos(tm) , x0[1] - (rm+hm)*sin(tm) , 0.0 };
  real_t x3[3] = { x0[0] - (rm+hm)*cos(tm) , x0[1] - (rm+hm)*sin(tm) , 0.0 };
  real_t x4[3] = { x0[0] -  rm    *cos(tm) , x0[1] -  rm    *sin(tm) , 0.0 };

  Node node0(context,x1);
  Node node1(context,x2);
  Node node2(context,x3);
  Node node3(context,x4);

  Edge edge0(context,node0,node1);
  Edge edge1(context,node1,node2);
  Edge edge2(context,node2,node3);
  Edge edge3(context,node3,node0);

  real_t xr[3] = { x0[0] + de*cos(te) , x0[1] + de*sin(te) , 0.0 };
  Circle eyeR( context , xr , re );
  Edge eyeR_edge( context , eyeR );

  real_t xl[3] = { x0[0] - de*cos(te) , x0[1] + de*sin(te) , 0.0 };
  Circle eyeL( context , xl , re );
  Edge eyeL_edge( context , eyeL );

  EdgeLoop loop(context,CLOSED);
  loop.add(edge0,1);
  loop.add(edge1,1);
  loop.add(edge2,1);
  loop.add(edge3,1);
  //loop.add(eyeR_edge,1);
  //loop.add(eyeL_edge,-1);
  loop.make();

  Face face(loop,SFORWARD);
  FaceBody body(context,face);
  object_ = body.object();

  set_object(object_);
  build_hierarchy();
}


#endif

} // EGADS

} // avro
