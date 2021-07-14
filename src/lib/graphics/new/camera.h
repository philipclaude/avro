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

	vec3& eye() { return eye_; }

  void set_eye( const vec3& eye ) { eye_ = eye; compute_view(); }
  void set_lookat( const vec3& lookat ) { lookat_ = lookat; compute_view(); }

  void compute_view();

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
