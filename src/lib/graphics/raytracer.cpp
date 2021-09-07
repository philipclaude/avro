#include "common/parallel_for.h"

#include "graphics/gl.h"
#include "graphics/raytracer.h"
#include "graphics/shader.h"

#include "mesh/points.h"
#include "mesh/topology.h"

#include "numerics/vec.hpp"

#include <time.h>

namespace avro
{

namespace graphics
{

Canvas::Canvas( int width , int height ) :
  width_(width),
  height_(height),
  pixel_(width_*height_)
{}

void
Canvas::init_gl() {

  // only do the following if opengl is supported
  GL_CALL( glGenVertexArrays( 1, &vertex_array_ ) );
  GL_CALL( glBindVertexArray(vertex_array_) );

  // bind the colormap values to a buffer
  GL_CALL( glGenBuffers( 1 , &pixel_buffer_ ) );
  GL_CALL( glGenTextures( 1 , &pixel_texture_ ) );

  // initialize the shader
  std::vector<std::string> macros = {"#version 330"};
  shader_ = std::make_shared<ShaderProgram>("raytracer",false,macros);
  shader_->use();

  shader_->setUniform("u_width",width_);
  shader_->setUniform("u_height",height_);

  // bind the desired texture
  glActiveTexture(GL_TEXTURE0 + 0);
  GLint pixels_location = glGetUniformLocation(shader_->handle() , "pixels");
  glUniform1i(pixels_location, 0); // first sampler in fragment shader
}

void
Canvas::draw_gl() {

  // clear the screen
  glClearColor(1,1,1,0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

  // write the pixel data to the pixel buffer
  GL_CALL( glBindBuffer( GL_TEXTURE_BUFFER , pixel_buffer_) );
  GL_CALL( glBufferData( GL_TEXTURE_BUFFER , sizeof(gl_float)*data_.size() , data_.data() , GL_STATIC_DRAW) );

  // bind the pixel buffer to the pixel texture
  GL_CALL( glActiveTexture( GL_TEXTURE0 + 0 ) );
  GL_CALL( glBindTexture( GL_TEXTURE_BUFFER , pixel_texture_ ) );
  GL_CALL( glTexBuffer( GL_TEXTURE_BUFFER , GL_RGB32F , pixel_buffer_ ) );

  // nothing is actually drawn here, we just rely on the interpolation to give (u,v) coordinates to then look up the texture value
  GL_CALL( glDrawArrays( GL_POINTS , 0 , 1 ) );
}

void
Canvas::convert() {

  data_.resize( width_*height_*3 , 0.0 );

  for (index_t i = 0; i < height_; i++)
  for (index_t j = 0; j < width_; j++)
  for (index_t d = 0; d < 3; d++)
    data_[3*(i*width_+j)+d] = gl_float( (*this)(i,j)[d] );
}

vec3
Sphere::normal( const vec3& point ) const {
  return { (point[0] - center_[0])/radius_ , (point[1] - center_[1])/radius_ , (point[2] - center_[2])/radius_ };
}

bool
Sphere::intersect( const Ray& ray , Intersection& hit , real_t tmin , real_t tmax ) const {

  // initialize to no intersection
  hit.t      = -1;
  hit.object = nullptr;

  // compute the coefficients in the quadratic equation for the intersection
  vec3 offset = ray.origin - center_;
  real_t b    = glm::dot( ray.direction , offset );
  real_t c    = glm::dot( offset , offset ) - radius_*radius_;
  real_t d    = b*b - c;
  if (d < 0.0) return false;

  // compute both possible intersections
  real_t t1 = -b - std::sqrt(d);
  real_t t2 = -b + std::sqrt(d);

  // check if t1 is in the admissible range (it is smaller), otherwise check t2)
  if      (t1 > tmin && t1 < tmax) hit.t = t1;
  else if (t2 > tmin && t2 < tmax) hit.t = t2;
  else return false;

  // save the intersection information
  hit.point  = ray.origin + hit.t * ray.direction;
  hit.normal = normal(hit.point);
  hit.object = this;

  return true;
}

bool
Scene::intersect( const Ray& ray , Intersection& hit ) const {

  real_t tmin = 1e-3;
  real_t tmax = 1e20;

  // default to no intersection
  hit.t = tmax;

  // check every object in the scene
  bool intersected = false;
  for (index_t k = 0; k < objects_.size(); k++) {

    // check the intersection with this object
    Intersection hit_k;
    if (objects_[k]->intersect( ray , hit_k , tmin , tmax) ) {

      // there is an intersection
      intersected = true;

      // is the intersection closer than the previous t value?
      if (hit_k.t < hit.t) {
        hit = hit_k;
      }
    }
  }

  return intersected;
}

RayTracer::RayTracer( int width , int height ) :
  window_(width,height),
  canvas_(width,height)
{
  window_.init();
  canvas_.init_gl();

  window_.camera().set_lookat( {0.,0.,0.} );
  //window_.camera().set_eye( {0.,2.,10.} );
  window_.camera().set_eye( {10.,5.,10.} );
  window_.camera().set_fov( M_PI/6.0 );
}

vec3
Material::shade( const Ray& ray , const Intersection& hit , const Light& light ) const {

  // direction to light
  vec3 l = glm::normalize( light.position - hit.point );

  // fraction of light reflected
  float cos_theta = std::max( 0.0f , glm::dot( l , hit.normal ) );

  // compute diffuse contribution
  vec3 cd = light.Ld * kd * cos_theta;

  // compute h = (l + v) / | l + v | (v = opposite ray direction)
  vec3 h = glm::normalize( l - ray.direction );

  // compute the fraction of light that is reflected due to the specular highlight
  vec3 cs;
  if (shine > 0.0) {
    float specular = std::max( 0.0f , std::pow( glm::dot(hit.normal,h) , shine ) );
    cs = light.Ls * specular;
  }

  return cd + cs;
}

vec3
reflect( const vec3& x , const vec3& y ) {
  float dp = 2.0 * glm::dot(x,y);
  return glm::normalize( x - dp*y );
}

vec3
refract( const vec3& v , const vec3& n , float n1_over_n2 ) {

  float dt = glm::dot(v,n);
  float discriminant = 1.0 - n1_over_n2*n1_over_n2*( 1.0 - dt*dt );
  if (discriminant < 0.0) {
    // total internal reflection
    return reflect(v,n);
  }

  vec3 r1 = n1_over_n2 * (v - dt * n);
  vec3 r2 = std::sqrt(discriminant) * n;

  return glm::normalize( r1 - r2 );
}

vec3
REFRACT( const vec3& v , const vec3& n , float eta ) {
  vec3 outward_normal;
  float n1_over_n2 = 0.0;

  if (glm::dot(v,n) > 0.0) {
    // exiting material: (n1/n2) = (material eta)
    n1_over_n2 = eta;
    outward_normal = -1.0 * n;
  }
  else {
    // entering material: (n1/n2) = 1/(material eta)
    n1_over_n2 = 1./eta;
    outward_normal = n;
  }
  return refract( v , outward_normal , n1_over_n2 );
}

bool
Material::scatter( const Ray& ray , const Intersection& hit , Ray& scatter_ray ) const {

  if (type == 0) return false;
  else if (type == 1) {
    // reflection ray
    scatter_ray.origin    = hit.point;
    scatter_ray.direction = reflect( ray.direction , hit.normal );
  }
  else if (type == 2) {
    // refraction ray
    scatter_ray.origin    = hit.point;
    scatter_ray.direction = REFRACT( ray.direction , hit.normal , eta );
  }
  else
    avro_implement;

  return true;
}

void
RayTracer::get_color( const Ray& ray , vec3& color , int depth ) {

  Intersection hit;

  // determine which objects in the scene are hit
  bool intersected = scene_.intersect( ray , hit );

  if (!intersected || depth == 0) {
    // no intersection, use the background color
    real_t t = 0.5*( ray.direction[1] ) + 0.2;
    vec3 cA = {1.,1.,1.};
    vec3 cB = {0.5,0.7,1.0};
    color = (1.0 - t)*cA + t*cB;
    return;
  }

  // check if we want to render the boundary (of a triangle)
  if (hit.boundary > 0) {
    color = {0.,0.,0.};
    return;
  }

  // initialize to color from ambient lighting
  color = hit.object->material().ka * ca_ * 0.5;

  // loop through all the lights
  for (index_t k = 0; k < lights_.size(); k++) {

    // check if we are in the shadow by intersecting a shadow ray with the objects in the scene
    Ray shadow;
    shadow.origin = hit.point;
    shadow.direction = glm::normalize( lights_[k].position - hit.point );
    Intersection hit_k;
    bool in_shadow = scene_.intersect( shadow , hit_k );
    if (in_shadow) continue;

    vec3 color_k;
    color_k = hit.object->material().shade(ray, hit, lights_[k]);

    color = color + color_k;
  }

  // determine if a ray is scattered
  Ray scatter_ray;
  bool scattered = hit.object->material().scatter( ray , hit , scatter_ray );

  if (scattered) {
    vec3 scatter_color;
    get_color( scatter_ray , scatter_color , depth-1 );
    color = color * 0.1 + scatter_color * 0.8;
  }
}

vec3
RayTracer::pixel2world( real_t u , real_t v ) const {

  // retrieve some parameters
  const real_t fov = window_.camera().fov();
  const real_t d = glm::norm( window_.camera().lookat() - window_.camera().eye() );
  const real_t a = real_t(canvas_.width())/real_t(canvas_.height());
  const real_t h = 2*d*tan(fov/2.0);
  const real_t w = a*h;

  // pixel coordiantes in camera space
  real_t pu = -0.5*w + w*u;
  real_t pv = -0.5*h + h*v;
  real_t pw = -d;
  vec3 q = {pu,pv,pw};

  // pixel coordinates in world space
  return basis_ * q + window_.camera().eye();
}

void
RayTracer::trace( index_t k ) {

  // get the pixel indices in [0,width] x [0,height]
  index_t i, j;
  pixel(k,i,j);

  // initialize the color
  canvas_(i,j) = {0.,0.,0.};

  index_t nb_samples_ = 8; // TODO make this a user parameter
  for (index_t s = 0; s < nb_samples_; s++) {

    // get the pixel coordinates in [0,1] x [0,1]
    real_t u = (j + random_within(0.,1.))/window_.width();
    real_t v = (i + random_within(0.,1.))/window_.height();

    // compute the world coordinates of the pixel
    vec3 p = pixel2world(u,v);

    // create a ray passing through the pixel
    Ray ray;
    ray.origin    = window_.camera().eye();
    ray.direction = glm::normalize( p - window_.camera().eye() );

    // determine which objects are intersected by the ray
    vec3 color;
    get_color( ray , color );

    canvas_(i,j) = canvas_(i,j) + color * (1./nb_samples_);
  }
}

void
RayTracer::render() {

  // compute the ambient light color
  ca_ = {0.,0.,0.};
  for (index_t k = 0; k < lights_.size(); k++) {
    ca_ = ca_ + lights_[k].La;
  }
  ca_ = ca_ * (1./lights_.size());

  // retrieve the view parameters
  const vec3& center = window_.camera().lookat();
  const vec3& eye = window_.camera().eye();
  vec3 up = {0.,1.,0.};

  // calculate the orthonormal basis for the camera
  vec3 g = center - eye;
  vec3 w = g * real_t(-1./glm::norm(g));
  vec3 u = glm::cross(up,w);
  u = u * (1./glm::norm(u));
  vec3 v = glm::cross(w,u);

  // save the orthonormal basis
  for (coord_t d = 0; d < 3; d++) {
    basis_(d,0) = u[d];
    basis_(d,1) = v[d];
    basis_(d,2) = w[d];
  }

  // cast a ray through each pixel
  clock_t t0 = clock();
  #if 0 // ray trace in serial
  for (index_t k = 0; k < nb_pixels(); k++) {
    trace(k);
  }
  index_t nb_thread = 1;
  #else // option to ray trace in parallel
  typedef RayTracer thisclass;
  ProcessCPU::parallel_for(
    parallel_for_member_callback( this , &thisclass::trace ), 0,nb_pixels()
  );
  index_t nb_thread = ProcessCPU::maximum_concurrent_threads();
  #endif
  printf("--> render time: %g sec.\n",(clock()-t0)/real_t(CLOCKS_PER_SEC*nb_thread) );

  // determine if we want to render to the OpenGL framebuffer, or to an image
  canvas_.convert();
  canvas_.draw_gl();
  glfwSwapBuffers(window_.window());

  #if 0
  // write a ppm file
  const real_t f = 255.99;
  FILE* fid = fopen("test.ppm","w");
  fprintf(fid,"P3\n%d %d\n255\n",canvas_.width(),canvas_.height());
  for (int i = canvas_.height()-1; i >= 0; i--) {
    for (index_t j = 0; j < canvas_.width(); j++) {
      real_t r = f*std::min( 1.0f , canvas_(i,j)[0] );
      real_t g = f*std::min( 1.0f , canvas_(i,j)[1] );
      real_t b = f*std::min( 1.0f , canvas_(i,j)[2] );
      fprintf(fid,"%d %d %d\n",int(r),int(g),int(b));
    }
  }
  fclose(fid);
  #endif
}

Triangle::Triangle( const Topology<Simplex>& triangles , index_t k , const Material& material ) :
  Object(material)
{
  for (index_t j = 0; j < 3; j++)
  for (index_t d = 0; d < 3; d++)
    vertex_[j][d] = triangles.points()[ triangles(k,j) ][d];

  box_.max = vertex_[0];
  box_.min = vertex_[0];
  for (index_t j = 1; j < 3; j++) {
    for (coord_t d = 0; d < 3; d++) {
      if (vertex_[j][d] < box_.min[d]) box_.min[d] = vertex_[j][d];
      if (vertex_[j][d] > box_.max[d]) box_.max[d] = vertex_[j][d];
    }
  }
}

bool
Triangle::intersect( const Ray& ray , Intersection& hit , real_t tmin , real_t tmax ) const {

  // setup the system of equations to determine the intersection point
  mat3 M;
  vec3 B;
  for (coord_t d = 0; d < 3; d++) {
    M(d,0) = vertex_[0][d] - vertex_[2][d];
    M(d,1) = vertex_[1][d] - vertex_[2][d];
    M(d,2) = -ray.direction(d);
    B(d)   = ray.origin(d) - vertex_[2][d];
  }

  // solve the system
  mat3 Minv = glm::inverse(M);
  vec3 C = Minv * B;

  // compute the barycentric coordinates
  float alpha = C(0);
  if (alpha < 0. || alpha > 1.) return false;
  float beta = C(1);
  if (beta < 0. || beta > 1.) return false;
  float gamma = 1.0 - alpha - beta;
  if (gamma < 0. || gamma > 1.) return false;

  // check for an admissible intersection
  float t = C(2);
  if (t > tmin && t < tmax) {

    hit.point = ray.origin + t * ray.direction;
    vec3 u,v;
    for (coord_t d = 0; d < 3; d++) {
      u(d) = vertex_[1][d] - vertex_[0][d];
      v(d) = vertex_[2][d] - vertex_[0][d];
    }
    hit.normal = glm::normalize( glm::cross(u,v) );
    hit.t = t;
    hit.object = this;

    // TODO calculate a tolerance that gives a constant thickness for the mesh edges
    real_t tol = 1e-2;
    hit.boundary = -1;
    if (alpha < tol || beta < tol || gamma < tol) hit.boundary = 1;

    return true;
  }
  return false;
}

CurvilinearTriangle::CurvilinearTriangle( const Topology<Simplex>& triangles , index_t k , const Material& material ) :
  Object(material)
{
  avro_implement;
}

bool
CurvilinearTriangle::intersect( const Ray& ray , Intersection& hit , real_t tmin , real_t tmax ) const {
  avro_implement;
}

void
Scene::add( const Topology<Simplex>& triangles , const Material& material , bool use_bvh ) {

  // we will need to extract triangles if this is a tetrahedral mesh
  if (triangles.number() != 2) avro_implement;
  avro_assert( triangles.points().dim() == 3 );

  index_t order = triangles.element().order();

  index_t k0 = items_.size();
  for (index_t k = 0; k < triangles.nb(); k++) {

    std::shared_ptr<Object> triangle;
    if (order == 1) triangle = std::make_shared<Triangle>(triangles,k,material);
    else triangle = std::make_shared<CurvilinearTriangle>(triangles,k,material);

    items_.push_back(triangle);
    triangle->name() = "triangle " + std::to_string(k);
  }

  // compute BVH of all triangles
  if (use_bvh) {
    std::shared_ptr<BVH_Node> bvh = std::make_shared<BVH_Node>(material,items_,k0,items_.size());
    objects_.push_back(bvh.get());
    nodes_.push_back(bvh);
  }
  else {
    for (index_t k = k0; k < items_.size(); k++)
      objects_.push_back( items_[k].get() );
  }
}

bool
box_compare(const std::shared_ptr<Object> a , const std::shared_ptr<Object> b, int axis) {
  AABB box_a;
  AABB box_b;
  if (!a->bounding_box(box_a) || !b->bounding_box(box_b))
    avro_assert_not_reached;
  return box_a.min[axis] < box_b.min[axis];
}

bool
box_x_compare (const std::shared_ptr<Object> a, const std::shared_ptr<Object> b) {
  return box_compare(a, b, 0);
}

bool
box_y_compare (const std::shared_ptr<Object> a, const std::shared_ptr<Object> b) {
  return box_compare(a, b, 1);
}

bool
box_z_compare (const std::shared_ptr<Object> a, const std::shared_ptr<Object> b) {
  return box_compare(a, b, 2);
}

AABB
surrounding_box(const AABB& box0, const AABB& box1) {
  vec3 small = { std::fmin(box0.min[0], box1.min[0]),
                 std::fmin(box0.min[1], box1.min[1]),
                 std::fmin(box0.min[2], box1.min[2]) };

  vec3 big  = { std::fmax(box0.max[0], box1.max[0]),
                std::fmax(box0.max[1], box1.max[1]),
                std::fmax(box0.max[2], box1.max[2]) };

  AABB aabb;
  aabb.min = small;
  aabb.max = big;
  return aabb;
}

BVH_Node::BVH_Node( const Material& material , const std::vector<std::shared_ptr<Object>>& src_objects , index_t start , index_t end ) :
  Object(material)
{
  int axis = int(random_within(0,2));
  std::vector<std::shared_ptr<Object>> objects = src_objects;
  auto comparator = (axis == 0) ? box_x_compare
                  : (axis == 1) ? box_y_compare
                                : box_z_compare;

  index_t object_span = end - start;
  if (object_span == 1 ) {
    left_ = right_ = objects[start];
  }
  else if (object_span == 2) {
    if (comparator(objects[start],objects[start+1])) {
      left_  = objects[start];
      right_ = objects[start+1];
    }
    else {
      left_  = objects[start+1];
      right_ = objects[start];
    }
  }
  else {

    std::sort( objects.begin() + start , objects.begin() + end , comparator );

    index_t mid = start + object_span/2;
    left_  = std::make_shared<BVH_Node>(material,objects, start, mid);
    right_ = std::make_shared<BVH_Node>(material,objects, mid, end);

    left_->name()  = "bvh node";
    right_->name() = "bvh node";
  }
  //printf("left = %s, right = %s\n",left_->name().c_str(),right_->name().c_str());

  // compute the bounding box of the left and right nodes
  AABB box_left, box_right;
  if (!left_->bounding_box(box_left) || !right_->bounding_box(box_right)) {
    avro_assert_not_reached;
  }

  // compute the bounding box of this node by accumulating the AABBs of the children
  box_ = surrounding_box(box_left,box_right);
}

bool
BVH_Node::bounding_box( AABB& output_box ) const {
  output_box = box_;
  return true;
}

bool
BVH_Node::intersect( const Ray& ray , Intersection& hit , real_t tmin , real_t tmax ) const {

  // if the ray doesn't hit this box, then it won't intersect the children boxes
  if (!box_.intersect(ray, hit, tmin, tmax))
    return false;

  // recursively determine if any of the children AABBs are intersected
  bool hit_left  = left_->intersect(ray, hit , tmin, tmax );
  bool hit_right = right_->intersect(ray, hit, tmin, hit_left ? hit.t : tmax );
  return hit_left || hit_right;
}

bool
AABB::intersect( const Ray& ray , Intersection& hit , real_t tmin , real_t tmax ) const {

  // loop through the box planes and determine if there is an intersection
  for (int d = 0; d < 3; d++) {
    real_t t0 = std::fmin((min[d] - ray.origin[d]) / ray.direction[d],
                          (max[d] - ray.origin[d]) / ray.direction[d]);
    real_t t1 = std::fmax((min[d] - ray.origin[d]) / ray.direction[d],
                          (max[d] - ray.origin[d]) / ray.direction[d]);
    tmin = std::fmax(t0, tmin);
    tmax = std::fmin(t1, tmax);
    if (tmax <= tmin)
      return false;
  }
  return true;
}

} // graphics

} // avro
