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
  Shaders( coord_t pmin , coord_t pmax , coord_t qmin , coord_t qmax ) {
    generate(pmin,pmax,qmin,qmax);
  }

  void generate( coord_t pmin , coord_t pmax , coord_t qmin , coord_t qmax) {

    for (coord_t p = pmin; p <= pmax; p++) {

      for (coord_t q = qmin; q <= qmax; q++) {

        // generate shaders with and without tessellation shaders

        printf("generating shader for p = %u, q = %u\n",p,q);

        std::string base = "-p" + std::to_string(p) + "-q" + std::to_string(q);

        std::vector<std::string> macros = {"#define SOLUTION_ORDER " + std::to_string(p),
                                           "#define GEOMETRY_ORDER " + std::to_string(q) };

        std::shared_ptr<ShaderProgram> t0 = std::make_shared<ShaderProgram>("triangles",true,macros);
        std::shared_ptr<ShaderProgram> e0 = std::make_shared<ShaderProgram>("edges",true,macros);

        std::shared_ptr<ShaderProgram> t1 = std::make_shared<ShaderProgram>("triangles",false,macros);
        std::shared_ptr<ShaderProgram> e1 = std::make_shared<ShaderProgram>("edges",false,macros);

        shaders_.insert( {"triangles" + base + "-tess=on" , t0 } );
        shaders_.insert( {"edges" + base + "-tess=on" , e0 } );

        shaders_.insert( {"triangles" + base + "-tess=off" , t1 } );
        shaders_.insert( {"edges" + base + "-tess=off" , e1 } );
      }
    }
    shaders_.insert( {"points" , std::make_shared<ShaderProgram>("points") } );
  }

  ShaderProgram& get( const std::string& type , const coord_t p , const coord_t q , bool with_tess=true ) {
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
