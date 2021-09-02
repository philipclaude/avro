#ifndef AVRO_LIB_GRAPHICS_RAYTRACER_H_
#define AVRO_LIB_GRAPHICS_RAYTRACER_H_

#include "avro_params.h"

#include "element/simplex.h"

#include "graphics/window.h"

#include <memory>
#include <vector>

namespace avro
{

template<typename type> class Topology;

namespace graphics
{

class ShaderProgram;
struct Material;
struct Ray;
struct Intersection;

struct AABB {
  real_t length;
  vec3 center;
  vec3 min;
  vec3 max;
};


class Object {
public:

  Object( const Material& material );
  virtual ~Object() {}

  virtual void intersect( const Ray& ray , Intersection& ixn ) const = 0;
  const Material& material() const { return material_; }
  const AABB& box() const { return box_; }

protected:
  const Material& material_;
  AABB box_;
};

class BVH_Node : public Object {

private:
  std::shared_ptr<Object> left_;
  std::shared_ptr<Object> right_;

  AABB box;
};

class BoundaryVolumeHierarchy {

public:

  void build( const std::vector<Object*>& objects );

private:
  std::vector<Object*> objects_;
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
  Object* object;
};

class Canvas {

public:
  Canvas( int width , int height );

  void init_gl();
  void draw_gl();

  vec3& operator() (int i , int j ) { return pixel_[i*height_+j]; }

  void to_ppm( const std::string& filename ) const;
  void to_framebuffer( gl_index framebuffer ) const;

  void convert();

private:
  int width_, height_;
  std::vector<vec3> pixel_;
  std::vector<gl_float> data_;
  gl_index pixel_buffer_;
  gl_index pixel_texture_;
  gl_index vertex_array_;

  std::shared_ptr<ShaderProgram> shader_;
};

class Sphere : public Object {
public:
  Sphere( real_t radius , const Material& material ) :
    Object(material),
    radius_(radius)
  {}

private:
  real_t radius_;
};

class Triangle : public Object {
public:
  Triangle( const Topology<Simplex>& topology , index_t k );

private:
  vec3 vertex_[3];
};

class RayTracer {

public:
  RayTracer(int width, int height);

  void add_objects( const TopologyBase& topology );
  void add_object( const Sphere& sphere );

  void build_bvh();
  void get_color( const Ray& ray );

  void trace( index_t k );
  void render();
  void draw();

  Window& window() { return window_; }

private:

  BoundaryVolumeHierarchy bvh_;
  Window window_;
  Canvas canvas_;
};

} // graphics

} // avro

#endif
