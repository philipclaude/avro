//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef AVRO_LIB_GRAPHICS_POSTSCRIPT_H_
#define AVRO_LIB_GRAPHICS_POSTSCRIPT_H_

#include "graphics/math.h"

#include <string>

namespace avro
{

namespace graphics
{

class BSPTriangle;

class PostScriptWriter {

public:
  PostScriptWriter( const std::string& filename );

  void begin( index_t width , index_t height );
  void write( const std::vector<BSPTriangle*>& triangles , const mat4& view_matrix , const mat4& projection_matrix , vec3 color = {0.5,0.5,0.5} );
  void end();

private:
  std::string filename_;
  mat4 screen_matrix_;
  FILE* fid_;

};

} // graphics

} // avro

#endif
