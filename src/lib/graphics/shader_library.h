//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef AVRO_LIB_GRAPHICS_SHADER_LIBRARY_H_
#define AVRO_LIB_GRAPHICS_SHADER_LIBRARY_H_

#include "graphics/shader.h"

#include <map>
#include <string>

namespace avro
{

namespace graphics
{

class Shaders {

public:
  Shaders( int pmin , int pmax , int qmin , int qmax ) {
    generate(pmin,pmax,qmin,qmax);
  }

  void generate( int pmin , int pmax , int qmin , int qmax) {

    #if AVRO_HEADLESS_GRAPHICS == 0
    std::string version = "410";
    #else
    std::string version = "330";
    #endif
    for (int p = pmin; p <= pmax; p++) {

      for (int q = qmin; q <= qmax; q++) {

        // generate shaders with and without tessellation shaders
        std::string base = "-p" + std::to_string(p) + "-q" + std::to_string(q);

        std::vector<std::string> macros = {"#version " + version, "#define SOLUTION_ORDER " + std::to_string(p),
                                           "#define GEOMETRY_ORDER " + std::to_string(q) ,
                                           "#define WITH_TESSELLATION 1" , "#define WITH_FIELD 1"};

        #if AVRO_HEADLESS_GRAPHICS == 0
        std::shared_ptr<ShaderProgram> t0 = std::make_shared<ShaderProgram>("triangles",true,macros);
        std::shared_ptr<ShaderProgram> e0 = std::make_shared<ShaderProgram>("edges",true,macros);
        shaders_.insert( {"triangles" + base + "-tess=on" , t0 } );
        shaders_.insert( {"edges" + base + "-tess=on" , e0 } );
        #endif

        macros[3] = "#define WITH_TESSELLATION 0";
        std::shared_ptr<ShaderProgram> t1 = std::make_shared<ShaderProgram>("triangles",false,macros);
        std::shared_ptr<ShaderProgram> e1 = std::make_shared<ShaderProgram>("edges",false,macros);
        shaders_.insert( {"triangles" + base + "-tess=off" , t1 } );
        shaders_.insert( {"edges" + base + "-tess=off" , e1 } );
      }
    }
    std::vector<std::string> macros = {"#version " + version};
    shaders_.insert( {"points" , std::make_shared<ShaderProgram>("points",false,macros) } );
  }

  ShaderProgram& get( const std::string& type , int p , int q , bool with_tess=true ) {
    std::string name = type + "-p" + std::to_string(p) + "-q" + std::to_string(q);
    if (with_tess) name += "-tess=on";
    else (name += "-tess=off");
    return get(name);
  }




  ShaderProgram& get( const std::string& name ) {
    if (shaders_.find(name) == shaders_.end()) {
      printf("could not load shader %s\n",name.c_str());
      avro_assert_not_reached;
    }
    return *shaders_.at(name).get();
  }

private:
  std::map<std::string,std::shared_ptr<ShaderProgram>> shaders_;
};

extern std::shared_ptr<Shaders> __shaders__;

} // graphics

} // avro

#endif
