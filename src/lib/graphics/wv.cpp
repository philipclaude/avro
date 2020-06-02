#include "graphics/application.h"
#include "graphics/wv.h"

#include <wsss.h>
#include <wsserver.h>

wvContext* __context__;
void* __web_server__;

#define WV_CHECK_STATUS(X) { if (X<0) { printf("error: wv returned code %d in file %s (line %d)\n",X,__FILE__,__LINE__); avro_assert(X>=0); } }

void
browserMessage( void* wsi , char* text , /*@unused*/ int )
{
  printf("browser message!\n");
  avro::graphics::Application<avro::graphics::Web_Interface> *app
    = (avro::graphics::Application<avro::graphics::Web_Interface>*)(__web_server__);
  app->receive(text);
}

namespace avro
{

namespace graphics
{

WV_Manager::WV_Manager()
{
  // set defaults
  bias_ = 0;
  fov_  = 30;
  znear_ = 1.;
  zfar_  = 100.;
  eye_[0] = eye_[1] = 0; eye_[2] = 7.;
  center_[0] = center_[1] = center_[2] = 0;
  up_[0] = up_[2] = 0; up_[1] = 1;

  // create context
  context_ = wv_createContext(0,fov_,znear_,zfar_,eye_,center_,up_);
}

void
WV_Manager::write( Primitive& primitive )
{
  if (primitive.points().size()==0) return;

  int status = 0;
  wvData items[5];
  float color[3] = {0,0,0}; // black

  char label[256];
  sprintf(label,"%p",(void*)&primitive);

  // send the vertex coordinates (everyone needs to do this)
  printf("primitive poitns = %lu\n",primitive.points().size());
  status = wv_setData( WV_REAL64 , primitive.points().size() , (void*)primitive.points().data() , WV_VERTICES , &items[0] );
  WV_CHECK_STATUS( status );

  //wv_adjustVerts( items[0] , focus ); // TODO retrieve focus...but from where??

  if (primitive.number()==0)
  {
    // send the vertex colours
    status = wv_setData( WV_REAL32 , 1 , (void*) color , WV_COLORS , &items[1] );

    if (index_.find(&primitive)!=index_.end())
    {
      status = wv_modGPrim( context_ , index_[&primitive] , 2 , items );
      WV_CHECK_STATUS( status );
      index_[&primitive] = status;
    }
    else
    {
      // create the gprim, but first check it really doesn't exist
      int idx = wv_indexGPrim(context_,label);
      if (idx<0)
      {
        idx = wv_addGPrim(context_,label , WV_POINT , WV_ON , 2 , items );
        index_.insert( {&primitive,idx} );
      }

      WV_CHECK_STATUS( idx );
      context_->gPrims[index_[&primitive]].pSize = 10;
    }
  }

  if (primitive.number()==1)
  {
    // add the line connectivity
    status = wv_setData( WV_INT32 , primitive.edges().size() , (void*) primitive.edges().data() , WV_INDICES , &items[1] );
    WV_CHECK_STATUS( status );

    // add the line colour
    status = wv_setData( WV_REAL32 , 1 , (void*) color , WV_COLORS , &items[2] );
    WV_CHECK_STATUS( status );

    if (index_.find(&primitive)!=index_.end())
    {
      status = wv_modGPrim( context_ , index_[&primitive] , 3 , items );
      WV_CHECK_STATUS( status );
      index_[&primitive] = status;
    }
    else
    {
      // create the gprim, but first check it really doesn't exist
      int idx = wv_indexGPrim(context_,label);
      if (idx<0)
      {
        idx = wv_addGPrim(context_,label , WV_LINE , WV_ON , 3 , items );
        index_.insert( {&primitive,idx} );
      }
      WV_CHECK_STATUS( idx );
    }
  }

  if (primitive.number()>=2)
  {
    // add the triangle connectivity
    std::vector<int> triangles( primitive.triangles().begin() , primitive.triangles().end() );
    status = wv_setData( WV_INT32 , triangles.size() , (void*)triangles.data() , WV_INDICES , &items[1] );
    WV_CHECK_STATUS( status );

    // add the triangle colors
    std::vector<float> colors( primitive.colors().begin() , primitive.colors().end() );
    status = wv_setData( WV_REAL32 , primitive.points().size(), (void*)colors.data() , WV_COLORS , &items[2] );
    WV_CHECK_STATUS( status );

    // add the  edges
    std::vector<int> edges( primitive.edges().begin() , primitive.edges().end() );
    status = wv_setData( WV_INT32 , edges.size() , (void*)edges.data() , WV_LINDICES , &items[3] );
    WV_CHECK_STATUS( status );

    std::vector<float> normals(primitive.normals().begin(),primitive.normals().end());
    status = wv_setData( WV_REAL32 , primitive.points().size() , (void*)normals.data() , WV_NORMALS , &items[4] );
    WV_CHECK_STATUS( status );

    if (index_.find(&primitive)!=index_.end())
    {
      status = wv_modGPrim( context_ , index_[&primitive] , 5 , items );
      WV_CHECK_STATUS( status );
      index_[&primitive] = status;
    }
    else
    {
      // create the gprim, but first check it really doesn't exist
      int idx = wv_indexGPrim(context_,label);
      if (idx<0)
      {
        idx = wv_addGPrim(context_,label , WV_TRIANGLE, WV_ON|WV_SHADING , 5 , items );
        index_.insert( {&primitive,idx} );
      }
      WV_CHECK_STATUS( idx );
    }
  }
}

} // graphics

} // avro
