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

#include <array>

namespace avro
{

namespace EGADS
{

Cube::Cube( const Context* context , coord_t number ) :
  Body(*context,&object_)
{}

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

    Body::object_ = wire.object();
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
  }
  this->build_hierarchy();
}

#ifndef AVRO_NO_ESP
void
SolidBody::make()
{
  EGADS_CHECK_SUCCESS( EG_makeSolidBody( *context_.get() , type_ , data_ , object_ ) );
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

#if 0
Plane::Plane( const Context* context , std::vector<real_t*>& x , real_t* uv )
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

  geometrics::normalize(u,3);
  geometrics::normalize(v,3);

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
Spline::Spline( const Context* _context , const endConditions end , const real_t tol , const std::vector<real_t>& points , bool _surface ) :
  context_(_context), surface_(_surface), end_(end), tol_(tol)
{
  int sizes[2] = { (int)points.size()/3 , surface_ ? 1 : 0 };

  int status = EG_approximate( *context_->get() , (int)end , tol , sizes , &points[0] , &object_ );
  if (status==EGADS_SUCCESS) ok_ = true;
  else ok_ = false;
  EGADS_CHECK_SUCCESS( status ); // for printing message;
}

Spline::Spline( const Context* _context , const endConditions end , const real_t tol , Topology<Simplex>& topology ) :
  context_(_context), end_(end), tol_(tol)
{

  Vertices& x = topology.vertices();

  if (topology.number()==1)
  {
    surface_ = false;

    // assume the topology edges are contiguous
    std::vector<real_t> points( 3*(topology.nb()+1) , 0. );
    for (index_t d=0;d<x.dim();d++)
      points[d] = x[ topology(0,0) ][d];
    for (index_t k=0;k<topology.nb();k++)
    {
      for (index_t d=0;d<x.dim();d++)
        points[3*(k+1)+d] = x[ topology(k,1) ][d];
    }

    int sizes[2] = { (int)points.size()/3 , 0 };

    int status = EG_approximate( *context_->get() , (int)end , tol , sizes , &points[0] , &object_ );
    if (status==EGADS_SUCCESS) ok_ = true;
    else ok_ = false;
    EGADS_CHECK_SUCCESS( status ); // for printing message;

  }
  else if (topology.number()==2)
  {
    surface_ = true;

    #if 1
    // fit the entire triangulation
    int ntri = topology.nb();

    int len = int(x.nb()); // not necessary to have all of them, but as a first pass ok
    std::vector<real_t> points(3*len,0.);

    for (index_t k=0;k<x.nb();k++)
    for (index_t d=0;d<x.dim();d++)
      points[3*k+d] = x[k][d];

    std::vector<int> tris( 3*ntri );
    for (index_t k=0;k<topology.nb();k++)
    for (index_t j=0;j<3;j++)
      tris[3*k+j] = topology(k,j)+1; // 1-bias in egads fitTriangles
    #else
    // fit only one triangle, resulting in a planar geometry
    int ntri = 1;
    int len = 3;
    std::vector<real_t> points(3*len,0);

    // use the first triangle coordinates
    for (index_t j=0;j<len;j++)
    for (coord_t d=0;d<x.dim();d++)
    {
      points[3*j+d] = x[topology(0,j)][d];
    }

    std::vector<int> tris = {1,2,3};
    #endif

    int status = EG_fitTriangles( *context_->get() , len , points.data() , ntri , tris.data() , NULL , tol_ , &object_ );
    if (status==EGADS_SUCCESS) ok_ = true;
    else ok_ = false;
    EGADS_CHECK_SUCCESS( status ); // for printing message;

    status = EG_getRange( object_ , range_ , &periodic_ );
    if (status==EGADS_SUCCESS) ok_ = true;
    else ok_ = false;
    EGADS_CHECK_SUCCESS( status ); // for printing message;
  }
  else
    avro_assert_not_reached;

}

Spline::endConditions
Spline::condition() const
{
  return end_;
}

real
Spline::tolerance() const
{
  return tol_;
}

const Context*
Spline::context() const
{
  return context_;
}

Face3D::Face3D( const Context* context , int style ) :
  Body(context,"makeFace3D")
{

  // construct the nodes
  real_t xyz_[3] = {0.,0.,0.};
  Node n0(context_,xyz_);

  xyz_[0] = 1.;
  Node n1(context_,xyz_);

  xyz_[0] = 0.;
  xyz_[1] = 2.;
  Node n2(context_,xyz_);

  xyz_[0] = 1.;
  Node n3(context_,xyz_);

  // construct the curves
  double p[3] = {0.,0.,0.};
  double dx[3] = {1.,0.,2.};
  double dy[3];
  real_t r;

  p[0] = .5;
  dx[0] = dy[2] = -1.;
  p[1] = p[2] = dx[1] = dx[2] = dy[0] = dy[1] = 0.;
  r = 0.5;
  Circlele circle0( context_ , p , dx , dy , r );

  p[0] = p[1] = p[2] = dx[2] = 0.;
  dx[0] = 0.;
  dx[1] = 2.;
  Line line1( context_ , p , dx );

  p[0]  = 1.;
  #if NEDGE==3
  dx[0] = -1.;
  #endif
  Line line2( context_ , p , dx );

  p[0] = 0.5;
  p[1] = 2.;
  p[2] = 0.;

  dx[0] = dy[2] = -1.;
  dx[1] = dx[2] = dy[0] = dy[1] = 0.;
  r = .5;
  Circlele circle3( context_ , p , dx , dy , r );

  // construct the edges
  Edge edge0(context_,circle0,n0,n1);
  Edge edge1(context_,n0,n2);

  #if NEDGE == 3
  Edge edge2(context_,n1,n2);
  #else
  Edge edge2(context_,circle3,n3,n2);
  Edge edge3(context_,n1,n3);
  #endif

  // make a surface without a loop
  EdgeLoop loop0(CLOSED);
  loop0.add( edge0 , -1 );
  loop0.add( edge1 , 1 );
  loop0.add( edge2 , -1 );
  #if NEDGE==4
  loop0.add( edge3 , -1 );
  #endif

  ego faceObject;

  // make the loop
  loop0.make( context_ );
  if (EG_isPlanar(loop0.object()) == 0) {

    Face face0(loop0,SFORWARD);
    faceObject = face0.object();

  }
  else
  {
    // get the surface that fits the loop
    EGADSIsocline surface0(loop0,style);
    loop0.deleteObject();

    for (index_t k=0;k<NEDGE;k++)
      loop0.addOther( surface0 , k );

    // update the loop
    loop0.make( context_ , surface0.object() );
    surface0.setObject( surface0.object() );

    Face face0( context_ , surface0 , loop0 , SFORWARD );
    faceObject = face0.object();

  }

  // make a face body
  EGADSBody body(context_,faceObject);
  setObject( body.object() );

  // retrieve the info about object class and model type
  EG_getInfo( object_ , &objectClass_ , &memberType_ , &top_ , &prev_ , &next_ );
  determine_number();

  buildHierarchy();

}

void
getPoints( Topology<Simplex>& topology , std::vector<real_t>& points )
{
  points.clear();
  std::vector<index_t> indices;
  for (index_t k=0;k<topology.nb();k++)
  {
    indices.push_back( topology(k,0) );
    if (k==topology.nb()-1) indices.push_back( topology(k,1) );
  }

  coord_t dlim = (topology.vertices().dim()<3) ? topology.vertices().dim() : 3;
  for (index_t k=0;k<indices.size();k++)
  {
    for (coord_t d=0;d<dlim;d++)
      points.push_back(topology.vertices()[indices[k]][d] );
    for (coord_t d=dlim;d<3;d++)
      points.push_back( 0. );
  }
}

Airfoil::Airfoil( const Context* context , const std::string& type , real_t chord , real_t* digits0 , index_t np  )
{
#if 0
  real_t digits[4] = {0,0,1,2};
  te_[0] = chord;
  te_[1] = te_[2] = 0.;
  if (digits0!=NULL)
  {
    for (index_t j=0;j<4;j++)
      digits[j] = digits0[j];
  }

  if (type=="naca")
  {
    //mesh = smart_new(NACA)(np,chord,digits);
    //topology = &mesh->topology(0);
  }
  else
  {
    avro_assert_not_reached;
  }
#else
  avro_implement;
#endif

  // interpolate the airfoil
  std::vector<real_t> points;
  Spline spline( context , Spline::Slope , 1e-7 , points , false );
  Edge edge(context,spline);

  EdgeLoop loop(CLOSED);
  loop.add(edge,1);
  loop.make();

  EGADSWireBody body(context,loop);
  setObject( body.object() );

  buildHierarchy();
}

void
Airfoil::addWake( const Context* context , const real_t length , real_t* dir0 )
{
  real_t dir[3] = {0,0,0};
  if (dir0==NULL)
  {
    dir[0] = 1.;
    dir[1] = 0.;
    dir[2] = 0.;
  }

  // get the node representing the trailing edge
  ego te = getTrailingEdge();

  real_t xyzWake[3];
  for (coord_t d=0;d<3;d++)
    xyzWake[d] = te_[d] +length*dir[d];

  // create the wake edge
  Node teNode(&te);
  Node wakeNode( context , xyzWake );
  Edge wake(context,teNode,wakeNode);

  // retrieve the airfoil edge
  //                   loop       edge
  Edge airfoil( entity(0)->child(0)->pobject() );

  // create an open loop with both airfoil and wake edges
  EdgeLoop loop(OPEN);
  loop.add( airfoil , -1 );
  loop.add( wake , -1 );
  loop.make();

  // make a wirebody from the single loop
  WireBody body(context,loop);
  setObject(body.object());

  // build the entity hierarchy
  buildHierarchy();
}

ego
Airfoil::getTrailingEdge()
{
  //     loop         edge     node
  return entity(0)->child(0)->child(0)->object();
}

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
  setObject(wire.object());

  buildHierarchy();
}
#endif

#endif

} // EGADS

} // avro
