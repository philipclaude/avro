//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_LIBRARY_EGADS_H_
#define avro_LIB_LIBRARY_EGADS_H_

#include "common/tools.h"

#include "geometry/egads/body.h"
#include "geometry/egads/context.h"
#include "geometry/egads/object.h"

#include "numerics/geometry.h"

namespace avro
{

namespace EGADS
{

class Context;

class Cube : public Body
{
public:
  Cube( const Context* context , coord_t dim );
  Cube( const Context* context , const std::vector<real_t>& lens , const real_t* x0=nullptr );

private:
  ego object_;
};

class Face3D : public Body
{
public:
  Face3D( const Context* context , int style );
};

class SolidBody : public Body
{
public:
  SolidBody( const Context* context , const int _type ) :
    Body(*context,nullptr),
    type_(_type) , count_(0)
  {}

  void make();

  void add_data( real_t* x , index_t nx )
  {
    for (index_t k=0;k<nx;k++)
      data_[count_++] = x[k];
  }
private:
  int type_;
  real_t data_[8];
  index_t count_;
};

class Sphere : public SolidBody
{
public:
  Sphere( const Context* context , real_t* x0 , real_t r );
};

class Cone : public SolidBody
{
public:
  Cone( const Context* context , real_t* apex , real_t *x0 , real_t r );
};

class Cylinder : public SolidBody
{
public:
  Cylinder( const Context* context , real_t *x0 , real_t* x1 , real_t r );
};

class Torus : public SolidBody
{
public:
  Torus( const Context* context , real_t* x0 , real_t* dir , real_t R , real_t r );
};

class NACA_Airfoil : public Body
{
public:
  NACA_Airfoil( const Context* context , real_t chord=1., real_t* digits0=NULL, index_t np=100 );

  void add_wake( const Context* context , const real_t length , real_t* dir=NULL );

protected:
  ego get_trailing_edge();
  real_t te_[3];
};


class Square : public Body
{
public:
  Square( const Context* context , const real_t* xc , const real_t lx , const real_t ly , const real_t theta=0. );
private:
  ego object_;
};

class Plane : public Object
{
public:

  // default constructor, basis is cartesian ex, ey about (0,0,0)
  Plane( const Context* context  ) :
    Object(*context)
  {
    std::vector<real_t> x = {0.,0.,0.};
    std::vector<real_t> e0 = {1.,0.,0.};
    std::vector<real_t> e1 = {0.,1.,0.};
    Plane(context,x.data(),e0.data(),e1.data());
  }

  // constructor from a basis e0, e1 about x0
  Plane( const Context* context , double *x , double *e0 , double *e1 ) :
    Object(*context)
  {
    double data[9] = { x[0],  x[1],  x[2],
                      e0[0], e0[1], e0[2],
                      e1[0], e1[1], e1[2] };
    #ifndef AVRO_NO_ESP
    EGADS_CHECK_SUCCESS( EG_makeGeometry( *context->get() , SURFACE , PLANE , NULL , NULL , data , &object_ ) );
    #else
    UNUSED(data);
    printf("need full EGADS to make geometry\n");
    avro_assert_not_reached;
    #endif
    construct(&object_);
  }

  Plane( const Context* context , std::vector<real_t*>& x , real_t* uv=NULL );

  real_t* range() { return range_; }

private:
  real_t range_[4];
  ego object_;
};

class Node : public Object
{
public:
  // constructor from cartesian coordinates
  Node( const Context* context , double *xyz ) :
    Object(*context)
  {
    for (coord_t d=0;d<3;d++)
      xyz_[d] = xyz[d];
    #ifndef AVRO_NO_ESP
    EGADS_CHECK_SUCCESS( EG_makeTopology( *context->get() , NULL , NODE , 0 , xyz , 0 , NULL , NULL , &object_ ) );
    #else
    printf("need full EGADS to make topology\n");
    avro_assert_not_reached;
    #endif
    construct(&object_);
  }

  Node( const Context* context , ego* node ) :
    Object(*context),
    object_(*node)
  {
    xyz_[0] = xyz_[1] = xyz_[2] = 0.; // or from EG_getGeometry ?
    construct(node);
  }

  real_t operator[] ( const coord_t d ) const { return xyz_[d]; }
  const real_t* x() const { return xyz_; }
  real_t* x() { return xyz_; }
private:
  real_t xyz_[3];
  ego object_;
};

class Curve : public Object
{
public:
  Curve( const Context* context ) :
    Object(*context)
  {}
protected:
  ego object_;
};

class Line : public Curve
{
public:
  // constructor from a cartesian coordinates and non-unit direction vector
  Line( const Context* context , const double* p0 , const double *dx ) :
    Curve(context)
  {
    double data[6] = { p0[0] , p0[1] , p0[2] , dx[0] , dx[1] , dx[2] };
    #ifndef AVRO_NO_ESP
    EGADS_CHECK_SUCCESS( EG_makeGeometry(*context->get(), CURVE , LINE , NULL , NULL , data , &object_ ) );
    #else
    UNUSED(data);
    printf("need full EGADS to make geometry\n");
    avro_assert_not_reached;
    #endif
    construct(&object_);
  }

  // constructor from two nodes
  Line( const Context* context , const Node& n0 , const Node& n1 ) :
    Curve(context)
  {
    double data[6] = { n0[0] , n0[1] , n0[2] , n1[0]-n0[0] , n1[1]-n0[1] , n1[2]-n0[2] };
    #ifndef AVRO_NO_ESP
    EGADS_CHECK_SUCCESS( EG_makeGeometry( *context->get() , CURVE , LINE , NULL , NULL , data , &object_ ) );
    #else
    UNUSED(data);
    printf("need full EGADS to make geometry\n");
    avro_assert_not_reached;
    #endif
    construct(&object_);
  }
};

class Circle : public Curve
{
public:
  Circle( const Context* context , const real_t* x , const real_t* dx , const real_t* dy , const real_t& r ) :
    Curve(context)
  {
    double prv[] = { x[0] , x[1] , x[2] , 1 , 0 , 0 , 0 , 1 , 0 , r };
    if (dx!=NULL)
    {
      avro_assert( dy!=NULL );
      prv[3] = dx[0]; prv[4] = dx[1]; prv[5] = dx[2];
      prv[6] = dy[0]; prv[7] = dy[1]; prv[8] = dy[2];
    }
    c_[0]  = x[0];  c_[1]  = x[1];  c_[2]  = x[2];
    dx_[0] = dx[0]; dx_[1] = dx[1]; dx_[2] = dx[2];
    dy_[0] = dy[0]; dy_[1] = dy[2]; dy_[2] = dy[2];
    r_ = r;
    EGADS_CHECK_SUCCESS( EG_makeGeometry( *context->get() , CURVE , CIRCLE , NULL , NULL , prv , &object_ ) );
    construct(&object_);
  }

  Circle( const Context* context , const real_t* x , const real_t r ) :
    Curve(context)
  {
    real_t dx[3] = {1.,0.,0.};
    real_t dy[3] = {0.,1.,0.};
    double prv[] = { x[0] , x[1] , x[2] , 1 , 0 , 0 , 0 , 1 , 0 , r };
    prv[3] = dx[0]; prv[4] = dx[1]; prv[5] = dx[2];
    prv[6] = dy[0]; prv[7] = dy[1]; prv[8] = dy[2];
    c_[0]  = x[0];  c_[1]  = x[1];  c_[2]  = x[2];
    dx_[0] = dx[0]; dx_[1] = dx[1]; dx_[2] = dx[2];
    dy_[0] = dy[0]; dy_[1] = dy[2]; dy_[2] = dy[2];
    r_ = r;
    EGADS_CHECK_SUCCESS( EG_makeGeometry( *context->get() , CURVE , CIRCLE , NULL , NULL , prv , &object_ ) );
    construct(&object_);
  }

  void point( real_t* p ) const
  {
    real_t d = std::sqrt( dx_[0]*dx_[0] +dx_[1]*dx_[1] +dx_[2]*dx_[2] );

    p[0] = c_[0] +dx_[0]*r_/d;
    p[1] = c_[1];
    p[2] = c_[2];
  }

private:
  real_t c_[3];
  real_t dx_[3];
  real_t dy_[3];
  real_t r_;
};

class Spline : public Curve
{
public:
  enum EndConditions { Natural , Slope , SlopeQuadratic };

  // constructor from a set of points, end condition and tolerance
  Spline( const Context* _context , const EndConditions end , const real_t tol , const std::vector<real_t>& points , bool _surface=false );

  EndConditions condition() const;
  real_t tolerance() const;
  const Context* context() const;

  real_t* range() { return range_; }

  bool ok() const { return ok_; }

private:
  const Context* context_;
  bool surface_;
  EndConditions end_;
  real_t tol_;
  real_t range_[4];
  bool ok_;
};

class Edge : public Object
{
public:

  // constructor from two nodes, creates an interim line in passing
  Edge( const Context* context , Node& n0 , Node& n1 ) :
    Object(*context)
  {
    real_t dummy[3];
    Line line(context,n0,n1);
    EGADS_CHECK_SUCCESS( EG_invEvaluate( *line.object() , n0.x() , &range_[0] , dummy ) );
    EGADS_CHECK_SUCCESS( EG_invEvaluate( *line.object() , n1.x() , &range_[1] , dummy ) );
    ego nodes[2] = {*n0.object(),*n1.object()};

    #ifndef AVRO_NO_ESP
    EGADS_CHECK_SUCCESS( EG_makeTopology( *context->get() , *line.object() , EDGE , TWONODE , range_ , 2 , nodes , NULL , &object_ ) );
    #else
    UNUSED(nodes);
    printf("need full EGADS to make topology\n");
    avro_assert_not_reached;
    #endif
    construct(&object_);
  }

  // constructor from either an arc or a line
  Edge( const Context* context , Curve& curve, Node& n0 , Node& n1 ) :
    Object(*context)
  {
    real_t dummy[3];
    EGADS_CHECK_SUCCESS( EG_invEvaluate( *curve.object() , n0.x() , &range_[0] , dummy ) );
    EGADS_CHECK_SUCCESS( EG_invEvaluate( *curve.object() , n1.x() , &range_[1] , dummy ) );
    if (range_[1]<range_[0]) range_[1] = 2.*range_[0];
    ego nodes[2] = {*n0.object(),*n1.object()};

    #ifndef AVRO_NO_ESP
    EGADS_CHECK_SUCCESS( EG_makeTopology( *context->get() , *curve.object() , EDGE , TWONODE , range_ , 2 , nodes , NULL , &object_ ) );
    #else
    UNUSED(nodes);
    printf("need full EGADS to make topology\n");
    avro_assert_not_reached;
    #endif
    construct(&object_);
  }

  // constructor from a full circle [0,2pi]
  Edge( const Context* context , Circle& circle ) :
    Object(*context)
  {
    int periodic;
    EGADS_CHECK_SUCCESS( EG_getRange( *circle.object() , &range_[0] , &periodic ) );
    avro_assert( periodic==1 );
    real_t p[3];
    circle.point(p);
    Node n0(context,p);
    #ifndef AVRO_NO_ESP
    EGADS_CHECK_SUCCESS( EG_makeTopology( *context->get() , *circle.object() , EDGE , ONENODE , &range_[0] , 1 , n0.object() , NULL , &object_ ) );
    #else
    printf("need full EGADS to make topology\n");
    avro_assert_not_reached;
    #endif
    construct(&object_);
  }

  Edge( const Context* context , Circle& circle , Node& n0 , Node& n1 ) :
    Object(*context)
  {
    // constructor for an arc
    real_t xyz[3];
    real_t d;
    EGADS_ENSURE_SUCCESS( EG_invEvaluate( *circle.object() , n0.x() , &range_[0] , xyz ) );
    d = numerics::distance( n0.x() , xyz , 3 );
    avro_assert( d < 1e-3 );
    EGADS_ENSURE_SUCCESS( EG_invEvaluate( *circle.object() , n1.x() , &range_[1] , xyz ) );
    d = numerics::distance( n1.x() , xyz , 3 );
    avro_assert( d < 1e-3 );
    ego nodes[2] = {*n0.object(),*n1.object()};
    #ifndef AVRO_NO_ESP
    EGADS_CHECK_SUCCESS( EG_makeTopology( *context->get() , *circle.object() , EDGE , TWONODE , &range_[0] , 2 , nodes , NULL , &object_ ) );
    #else
    printf("need full EGADS to make topology\n");
    avro_assert_not_reached;
    #endif
    construct(&object_);
  }

  // constructor from a spline
  Edge( const Context* context , Spline& spline , real_t* x0=NULL ) :
    Object(*context),
    object_(*spline.object())
  {

    range_[0] = 0.;
    range_[1] = 1.;

    ego node = NULL;
    if (x0!=NULL)
    {
      Node node0( context , x0 );
      node = *node0.object();
    }
    else
    {
      std::vector<real_t> X(3);
      std::vector<real_t> U(range_,range_+2);
      spline.evaluate(U,X);
      Node node0( context , X.data() );
      node = *node0.object();
    }
    #ifndef AVRO_NO_ESP
    EGADS_CHECK_SUCCESS( EG_makeTopology( *context->get() , *spline.object() , EDGE , ONENODE , range_ , 1 , &node , NULL , &object_ ) );
    #else
    printf("need full EGADS to make topology\n");
    avro_assert_not_reached;
    #endif
    construct(&object_);
  }

  Edge( const Context *context , ego* object ) :
    Object(*context),
    object_(*object)
  {
    construct(&object_);
  }

  real_t range( const coord_t d ) const { return range_[d]; }

private:
  real_t range_[2];
  ego    object_;
};

class Surface : public Object
{
public:

  Surface( const Context* context ) :
    Object(*context)
  {}

  void other_curve( ego& input , ego& output , double tol=0. )
  {
    EGADS_CHECK_SUCCESS( EG_otherCurve( object_ , input , tol , &output ) );
  }

protected:
  ego object_;
};

class EdgeLoop : public Object
{
public:
  EdgeLoop( const Context* context , int member_type=-1 ) :
    Object(*context),
    context_(context)
  {
    data_.member_type = member_type;
  }

  void make( ego ref=NULL )
  {
    //avro_assert( egos_.size() == senses_.size() );
    #ifndef AVRO_NO_ESP
    EGADS_ENSURE_SUCCESS( EG_makeTopology(*context_->get(),ref,LOOP,data_.member_type,NULL, senses_.size()  , &egos_[0] , &senses_[0] , &object_ ) );
    #else
    printf("need full EGADS to make topology\n");
    avro_assert_not_reached;
    #endif
    construct(&object_);
  }

  void add( Edge& obj , const int sense )
  {
    egos_.push_back( *obj.object() );
    senses_.push_back(sense);
  }

  void add_other( Surface& surface , const index_t k )
  {
    ego other;
    surface.other_curve( egos_[k] , other );
    egos_.push_back( other );
    //senses_.push_back( senses_[k] );
  }

  std::vector<int>& senses() { return senses_; }

  ego* operator[] ( const index_t k ) { return &egos_[k]; }
  int operator() ( const index_t k ) const { return senses_[k]; }

  const Context* context() const { return context_; }

private:
  const Context* context_;
  ego object_;
  std::vector<ego>  egos_;
  std::vector<int>  senses_;
};

class Isocline : public Surface
{
public:
  Isocline( EdgeLoop& loop , int style ) :
    Surface(loop.context())
  {
    EGADS_CHECK_SUCCESS( EG_isoCline(*loop.object(), style , 0. , &object_ ) );
    construct(&object_);
  }
};

class Face : public Object
{
public:
  Face( EdgeLoop& loop , int memberType ) :
    Object( *loop.context() )
  {
    // the easy way for planar faces
    EGADS_CHECK_SUCCESS( EG_makeFace( *loop.object() , memberType , NULL , &object_ ) );
    construct(&object_);
  }

  Face( const Context* context , Spline& surface , int memberType ) :
    Object(*context)
  {
    EGADS_CHECK_SUCCESS( EG_makeFace( *surface.object() , memberType , surface.range() , &object_ ) );
    construct(&object_);
  }

  Face( const Context* context , Plane& plane , int memberType ) :
    Object(*context)
  {
    EGADS_CHECK_SUCCESS( EG_makeFace( *plane.object() , memberType , plane.range() , &object_ ) );
    construct(&object_);
  }

  Face( const Context* context , Isocline& surface , EdgeLoop& loop , int memberType ) :
    Object(*context)
  {
    // the hard way for non-planar faces
    #ifndef AVRO_NO_ESP
    EGADS_CHECK_SUCCESS( EG_makeTopology( *context->get() , *surface.object() , FACE , memberType , NULL , 1 , loop.object() , loop.senses().data() , &object_ ) );
    #else
    printf("need full EGADS to make topology\n");
    avro_assert_not_reached;
    #endif
    construct(&object_);
  }

private:
  ego object_;
};

class Shell : public Object
{
public:
  Shell( const Context* context , ego* face , index_t nface ) :
    Object(*context,&object_)
  {
    #ifndef AVRO_NO_ESP
    EGADS_CHECK_SUCCESS( EG_makeTopology(*context->get(), NULL, SHELL ,CLOSED, NULL, nface , face , NULL, &object_ ));
    #else
    printf("need full EGADS to make topology\n");
    avro_assert_not_reached;
    #endif
  }
private:
  ego object_;
};

class SheetBody : public Body
{
public:
  SheetBody( const Context* context , Shell& shell ) :
    Body(*context)
  {
    #ifndef AVRO_NO_ESP
    EGADS_CHECK_SUCCESS( EG_makeTopology( *context->get() , NULL , BODY, SHEETBODY , NULL , 1 , shell.object() , NULL, &object_) );
    #else
    printf("need full EGADS to make topology\n");
    avro_assert_not_reached;
    #endif
    set_object( object_ );
  }
private:
  ego object_;
};

class FaceBody : public Body
{
public:
  FaceBody( const Context* context , Face& face ) :
    Body(*context)
  {
    // face body from a single face
    #ifndef AVRO_NO_ESP
    EGADS_CHECK_SUCCESS( EG_makeTopology( *context->get() , NULL , BODY, FACEBODY , NULL , 1 , face.object() , NULL, &object_) );
    #else
    printf("need full EGADS to make topology\n");
    avro_assert_not_reached;
    #endif
    set_object( object_ );
  }
private:
  ego object_;
};

class WireBody : public Body
{
public:

  WireBody( const Context* context , Edge edge ) :
    Body(*context)
  {
    EdgeLoop loop(context,CLOSED);
    loop.add(edge,1);
    loop.make();
    EGADS_CHECK_SUCCESS( EG_makeTopology(*context->get(),NULL, BODY , WIREBODY , NULL , 1 , loop.object() , NULL , &object_ ) );
    set_object( object_ );
  }

  WireBody( const Context* context , ego* edge , index_t nedge ) :
    Body(*context)
  {
    EdgeLoop loop(context,CLOSED);
    for (index_t k=0;k<nedge;k++)
    {
      Edge e(context,&edge[k]);
      loop.add(e,1);
    }
    loop.make();
    #ifndef AVRO_NO_ESP
    EGADS_CHECK_SUCCESS( EG_makeTopology(*context->get(),NULL, BODY , WIREBODY , NULL , 1 , loop.object() , NULL , &object_ ) );
    #else
    printf("need full EGADS to make topology\n");
    avro_assert_not_reached;
    #endif
    set_object( object_ );
  }

  WireBody( const Context* context , EdgeLoop& loop ) :
    Body(*context)
  {
    #ifndef AVRO_NO_ESP
    EGADS_CHECK_SUCCESS( EG_makeTopology(*context->get(),NULL, BODY , WIREBODY , NULL , 1 , loop.object() , NULL , &object_ ) );
    #else
    printf("need full EGADS to make topology\n");
    avro_assert_not_reached;
    #endif
    set_object( object_ );
  }
private:
  ego object_;

};

class Smiley : public Body
{
public:
  Smiley( const Context* context , real_t* x0 , real_t rf , real_t rm , real_t hm , real_t tm , real_t re , real_t de , real_t te );

private:
  ego object_;
};

} // EGADS

} // avro

#endif
