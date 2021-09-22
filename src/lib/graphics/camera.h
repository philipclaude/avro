//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef AVRO_LIB_GRAPHICS_CAMERA_H_
#define AVRO_LIB_GRAPHICS_CAMERA_H_

#include "graphics/math.h"

namespace avro
{

namespace graphics
{

class Camera {

public:
	Camera(float fov , index_t width , index_t height);

	float fov() const { return fov_; }
	const vec3& eye() const { return eye_; }
  const vec3& lookat() const { return lookat_; }

	void set_fov( float fov ) { fov_ = fov; }
  void set_eye( const vec3& eye ) { eye_ = eye; compute_view(); }
  void set_lookat( const vec3& lookat ) { lookat_ = lookat; compute_view(); }

  void compute_view();
  void compute_projection( index_t width , index_t height );

  const mat4& view_matrix() const { return view_matrix_; }
  const mat4& projection_matrix() const { return projection_matrix_; }

private:
  float fov_;
  index_t width_;
  index_t height_;

	vec3 eye_;
	vec3 up_;
	vec3 lookat_;

	mat4 projection_matrix_;
	mat4 view_matrix_;
};

} // graphics

} // avro

#endif
