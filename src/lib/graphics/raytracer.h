#ifndef AVRO_LIB_GRAPHICS_RAYTRACER_H_
#define AVRO_LIB_GRAPHICS_RAYTRACER_H_

#include "avro_params.h"

#include "element/simplex.h"

#include "graphics/window.h"

#include <algorithm>
#include <memory>
#include <string>
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
struct Light;

struct AABB {
  vec3 min;
  vec3 max;
  bool intersect( const Ray& ray , Intersection& ixn , real_t tmin , real_t tmax ) const;
};

class Object {
public:

  Object( const Material& material ) :
    material_(material)
  {}
  virtual ~Object() {}

  virtual bool intersect( const Ray& ray , Intersection& hit , real_t tmin , real_t tmax ) const = 0;
  virtual bool bounding_box( AABB& output ) const = 0;

  const Material& material() const { return material_; }
  std::string& name() { return name_; }

protected:
  const Material& material_;
  AABB box_;
  std::string name_;
};

class BVH_Node : public Object {

public:
  BVH_Node( const Material& material , const std::vector<std::shared_ptr<Object>>& objects , index_t start , index_t end );

  virtual bool intersect( const Ray& ray , Intersection& hit , real_t tmin , real_t tmax ) const override;
  virtual bool bounding_box( AABB& output ) const override;

private:
  std::shared_ptr<Object> left_;
  std::shared_ptr<Object> right_;
};

struct Ray {
  vec3 origin;
  vec3 direction;
};

struct Material {
  int type = 0; // 0 = diffuse, 1 = reflective, 2 = refractive
  vec3 ka;
  vec3 kd;  // diffusive color
  vec3 ks;  // specular color
  real_t eta; // index of refraction
  float shine = 32.0;

  vec3 shade( const Ray& ray , const Intersection& hit , const Light& light ) const;
  bool scatter( const Ray& ray , const Intersection& hit , Ray& scatter_ray ) const;
};

struct RefractiveMaterial : Material {

  vec3 scatter() const;
};

struct Intersection {
  real_t t;
  vec3 point;
  vec3 normal;
  const Object* object;
  int boundary = -1;
};

struct Light {
  vec3 La;
  vec3 Ld;
  vec3 Ls;
  vec3 position;
};

class Canvas {

public:
  Canvas( int width , int height );

  void init_gl();
  void draw_gl();

  vec3& operator() (int i , int j ) { return pixel_[i*width_+j]; }

  void to_ppm( const std::string& filename ) const;
  void to_framebuffer( gl_index framebuffer ) const;

  void convert();

  int width() const { return width_; }
  int height() const { return height_; }

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
  Sphere( vec3 center , real_t radius , const Material& material ) :
    Object(material),
    center_(center),
    radius_(radius)
  {}

  bool intersect( const Ray& ray , Intersection& ixn , real_t tmin , real_t tmax ) const;

  vec3 normal( const vec3& point ) const;
  bool bounding_box( AABB& output_box ) const {
    avro_implement;
    return false;
  }

private:
  vec3 center_;
  real_t radius_;
};

class Triangle : public Object {
public:
  Triangle( const Topology<Simplex>& topology , index_t k , const Material& material );

  virtual bool intersect( const Ray& ray , Intersection& ixn , real_t tmin , real_t tmax ) const override;
  virtual bool bounding_box( AABB& output_box ) const override {
    output_box = box_;
    return true;
  }

private:
  vec3 vertex_[3];
};

class CurvilinearTriangle : public Object {
public:
  CurvilinearTriangle( const Topology<Simplex>& topology , index_t k , const Material& material );

  virtual bool intersect( const Ray& ray , Intersection& ixn , real_t tmin , real_t tmax ) const override;
  virtual bool bounding_box( AABB& output_box ) const override {
    avro_implement;
    return false;
  }

private:
  std::vector<vec3> vertex_;

};

class Scene {
public:

  bool intersect( const Ray& ray , Intersection& ixn ) const;

  void add( const Topology<Simplex>& topology , const Material& material , bool use_bvh = true );
  void add( const Sphere& sphere ) {
    objects_.push_back(&sphere);
  }

private:
  std::vector<const Object*> objects_;
  //std::vector< std::shared_ptr<Triangle> > triangles_;
  std::vector< std::shared_ptr<Object> > items_;

  std::vector< std::shared_ptr<BVH_Node> > nodes_;
};

class RayTracer {

public:
  RayTracer(int width, int height);

  vec3 pixel2world( real_t u , real_t v ) const;

  void get_color( const Ray& ray , vec3& color , int depth = 10 );
  void trace( index_t k );
  void render();

  Window& window() { return window_; }
  Scene& scene() { return scene_; }

  index_t nb_pixels() const { return window_.width() * window_.height(); }
  void pixel( index_t k , index_t& i , index_t& j ) const {
    i = k / window_.width();
    j = k % window_.width();
  }

  std::vector<Light>& lights() { return lights_; }

private:
  mat3 basis_;
  Scene scene_;
  Window window_;
  Canvas canvas_;
  vec3 ca_; // ambient light color
  std::vector<Light> lights_;
};

} // graphics

} // avro

#endif
