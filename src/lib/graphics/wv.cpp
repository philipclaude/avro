//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
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

  index_t nb_points = primitive.points().size()/3;

  // send the vertex coordinates
  // no need to ask wv to focus them (with wv_adjustVerts) since we've already normalized the data
  status = wv_setData( WV_REAL64 , nb_points , (void*)primitive.points().data() , WV_VERTICES , &items[0] );
  WV_CHECK_STATUS( status );

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
    std::vector<int> edges( primitive.edges().begin() , primitive.edges().end() );
    status = wv_setData( WV_INT32 , edges.size() , (void*) edges.data() , WV_INDICES , &items[1] );
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
    for (index_t k=0;k<colors.size();k++)
      colors[k] /= 256.0;
    status = wv_setData( WV_REAL32 , nb_points , (void*)colors.data() , WV_COLORS , &items[2] );
    WV_CHECK_STATUS( status );

    // add the  edges
    std::vector<int> edges( primitive.edges().begin() , primitive.edges().end() );
    status = wv_setData( WV_INT32 , edges.size() , (void*)edges.data() , WV_LINDICES , &items[3] );
    WV_CHECK_STATUS( status );

    std::vector<float> normals(primitive.normals().begin(),primitive.normals().end());
    status = wv_setData( WV_REAL32 , nb_points , (void*)normals.data() , WV_NORMALS , &items[4] );
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

void
WV_Manager::write( const std::string& name , coord_t number , const std::vector<real_t>& points0 , const std::vector<index_t>& edges0 , const std::vector<index_t>& triangles0 , const std::vector<real_t>& colors0 )
{
  int status = 0;
  wvData items[4];

  index_t nb_points = points0.size()/3;

  // send the vertex coordinates
  // no need to ask wv to focus them (with wv_adjustVerts) since we've already normalized the data
  std::vector<real_t> points(points0.begin(),points0.end());
  status = wv_setData( WV_REAL64 , nb_points , (void*)points.data() , WV_VERTICES , &items[0] );
  WV_CHECK_STATUS( status );

  if (number>=2)
  {
    // add the triangle connectivity
    std::vector<int> triangles( triangles0.begin() , triangles0.end() );
    status = wv_setData( WV_INT32 , triangles.size() , (void*)triangles.data() , WV_INDICES , &items[1] );
    WV_CHECK_STATUS( status );

    // add the triangle colors
    std::vector<float> colors( colors0.begin() , colors0.end() );
    //for (index_t k=0;k<colors.size();k++)
    //  colors[k] /= 256.0;
    status = wv_setData( WV_REAL32 , nb_points , (void*)colors.data() , WV_COLORS , &items[2] );
    WV_CHECK_STATUS( status );

    // add the  edges
    if (edges0.size()>0)
    {
      std::vector<int> edges( edges0.begin() , edges0.end() );
      status = wv_setData( WV_INT32 , edges.size() , (void*)edges.data() , WV_LINDICES , &items[3] );
      WV_CHECK_STATUS( status );
    }

    if (aux_index_.find(name)!=aux_index_.end())
    {
      status = wv_modGPrim( context_ , int(aux_index_[name]) , 4 , items );
      WV_CHECK_STATUS( status );
      aux_index_[name] = status;
    }
    else
    {
      // create the gprim, but first check it really doesn't exist
      int idx = wv_indexGPrim(context_, const_cast<char*>(name.c_str()) );
      if (idx<0)
      {
        idx = wv_addGPrim(context_, const_cast<char*>(name.c_str()) , WV_TRIANGLE, WV_ON|WV_SHADING|WV_TRANSPARENT , 4 , items );
        aux_index_.insert( {name,idx} );
      }
      WV_CHECK_STATUS( idx );
    }
  }
}

void
WV_Manager::remove( const std::string& name )
{
  avro_assert( aux_index_.find(name)!=aux_index_.end() );
  index_t idx = aux_index_[name];
  aux_index_.erase(name);
  wv_removeGPrim(context_, int(idx) );
}

} // graphics

} // avro
