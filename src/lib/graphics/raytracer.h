#ifndef AVRO_LIB_GRAPHICS_RAYTRACER_H_
#define AVRO_LIB_GRAPHICS_RAYTRACER_H_

#include "graphics/window.h"

#include <memory>
#include <vector>

namespace avro
{

namespace graphics
{

class SceneObject;

class BoundaryVolumeHierarchy {

public:

  void build( const std::vector<SceneObject*>& objects );

private:
  std::vector<SceneObject*> objects_;
};

struct Ray {
  vec3 origin;
  vec3 direction;
};

struct Material {
  int type; // 0 = diffuse, 1 = reflective, 2 = refractive
  vec3 kd;  // diffusive color
  vec3 ks;  // specular color
  real_t eta; // index of refraction
};

struct Intersection {

  real_t t;
  real_t normal;
  SceneObject* object;
};

struct AABB {
  real_t length;
  vec3 center;
  vec3 min;
  vec3 max;
};

class Canvas {

public:
  Canvas( int width , int height ) :
    width_(width),
    height_(height),
    pixel_(width_*height_)
  {}

  vec3& operator() (int i , int j ) { return pixel_[i*height_+j]; }

  void to_bytes();

  void to_ppm( const std::string& filename ) const;
  void to_framebuffer( gl_index framebuffer ) const;

private:
  int width_, height_;
  std::vector<vec3> pixel_;
  std::vector<char> bytes_;
};

class SceneObject {
public:

  SceneObject( const Material& material );
  virtual ~SceneObject() {}

  virtual void intersect( const Ray& ray , Intersection& ixn ) const = 0;
  const Material& material() const { return material_; }
  const AABB& box() const { return box_; }

protected:
  const Material& material_;
  AABB box_;
};

class RayTracer {

public:
  RayTracer();

  void build_bvh();
  void get_color( const Ray& ray );

  void trace( index_t k );
  void draw( Canvas& canvas );

private:
  Window window_;
};

} // graphics

} // avro

#endif
