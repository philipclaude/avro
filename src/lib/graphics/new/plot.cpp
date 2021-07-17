#include "graphics/new/plot.h"
#include "graphics/new/primitives.h"
#include "graphics/new/vao.h"

#include "mesh/topology.h"

namespace avro
{

namespace graphics
{

ClipPlane::ClipPlane(	const Plot& plot) :
	plot_(plot),
	visible_(false),
	direction_(1),
	style_(0),
	dimension_(2), // z direction
	distance_(0) {
	transform_matrix_ = glm::identity();
}

void
ClipPlane::initialize( const vec3& c , float length ) {

	center_ = c;
	length_scale_ = length;

	// initialize with normal in the z-direction
	vec3 u = {1,0,0};
	vec3 v = {0,1,0};
	for (coord_t d = 0; d < 3; d++) {
		coordinates_[3*0+d] =   length_scale_*u(d) - length_scale_*v(d);
		coordinates_[3*1+d] =   length_scale_*u(d) + length_scale_*v(d);
		coordinates_[3*2+d] = - length_scale_*u(d) + length_scale_*v(d);
		coordinates_[3*3+d] = - length_scale_*u(d) - length_scale_*v(d);
	}

	// write the buffers to gl manager
	write();
}

void
ClipPlane::get( vec3& c , vec3& n ) const {
  n.zero();
  n(dimension_) = direction_;
  c = center_;
  c(dimension_) += distance_;
}

void
ClipPlane::update() {
  // translation from the origin to the clip center
  vec3 t = center_;
  t(dimension_) += distance_;

  // axis of rotation, depending on the clip normal direction
  vec3 axis;
  if (dimension_ == 0) axis = {0,1,0};
  if (dimension_ == 1) axis = {1,0,0};
  if (dimension_ == 2) axis = {0,0,1};

  // the transformation of the clip in world space is a rotation, followed by a translation
  transform_matrix_ = glm::translate( glm::identity() , t ) * glm::rotate(glm::identity(),M_PI/2.0,axis);
}

Plot::Plot( const TopologyBase& topology ) :
	topology_(topology),
	clip_(*this),
	length_scale_(1.0),
	hidden_(false) {
  model_matrix_ = glm::identity();
}

void
Plot::build() {
  if (topology_.number() == 4) {
    // build a spacetime topology and add the vao for each boundary (this is the topology we pass into build)

    // we discard the spacetime topology, but hold on to the actual 4d topology
    // in case the final vao (the slice) needs to be recomputed from the 4d topology
    avro_implement;
  }
  else {
    vao_.push_back( std::make_shared<VertexAttributeObject>() );
    active_vao_ = vao_[0].get();
    active_vao_->build(topology_);
  	vao_labels_.push_back("group 0");
  }
  compute_center();

  clip_.initialize( center_ , length_scale_ );
}

void
Plot::compute_center() {

  vec3 xmin = { 1e20, 1e20, 1e20};
  vec3 xmax = {-1e20,-1e20,-1e20};

  // compute the center
  center_.zero();
  const std::vector<gl_float>& coordinates = active_vao_->points().coordinates();
  for (index_t k = 0; k < coordinates.size()/3; k++)
  for (coord_t d = 0; d < 3; d++) {
    gl_float x = coordinates[3*k+d];
    center_(d) += x;

    if (x < xmin(d)) xmin(d) = x;
    if (x > xmax(d)) xmax(d) = x;
  }

  length_scale_ = -1;
  for (coord_t d = 0; d < 3; d++) {
    center_(d) *= (3./coordinates.size());
    float L = (xmax(d) - xmin(d));
    if (L > length_scale_) length_scale_ = L;
  }

  // store the center translation and its inverse
  center_translation_ =  glm::translate( glm::identity() , center_ );
  inverse_center_translation_ = glm::inverse(center_translation_);
}

void
Plot::transform_center(const mat4& m) {
  center_translation_ = m * center_translation_;
  for (coord_t d = 0; d < 3; d++)
    center_(d) = center_translation_(d,3);
  inverse_center_translation_ = glm::inverse(center_translation_);
}

void
Plot::transform( const mat4& m , bool centered) {
  // m is the transformation relative to the trackball
  if (centered) {
    model_matrix_ = center_translation_ * m * inverse_center_translation_ * model_matrix_;
  }
  else {
    model_matrix_ = m * model_matrix_;
  }
}

} // graphics

} // avro
