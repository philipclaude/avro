#ifndef AVRO_LIB_GRAPHICS_PLOT_H_
#define AVRO_LIB_GRAPHICS_PLOT_H_

#include "graphics/gl.h"
#include "graphics/math.h"

#include <vector>

namespace avro
{

class TopologyBase;

namespace graphics
{

class Plot;
class VertexAttributeObject;

class ClipPlane {
public:

	ClipPlane(const Plot& plot);
	~ClipPlane();

	void initialize( const vec3& c , float length );
	void write();
	void draw(const mat4& view , const mat4& projection) const;
	void get( vec3& c , vec3& n ) const;
	void update();

	bool  visible() const { return visible_; }
	bool& visible()       { return visible_; }
	int& style() { return style_; }
	int style() const { return style_; }
	float length_scale() const { return length_scale_; }
	int& dimension() { return dimension_; }
	float& distance() { return distance_; }
	void flip() { direction_ *= -1; }

private:
	const Plot& plot_;
	bool visible_;
	vec3 center_;
	float length_scale_;
	mat4 transform_matrix_; // transformation with respect to the initial position
	short direction_;
	int style_; // 0 for none, 1 for pixel, 2 for primitive
	int dimension_;
	float distance_;

	gl_index point_buffer_;
	gl_index index_buffer_;
	gl_index vertex_array_;

	gl_float coordinates_[12];
	const gl_index indices_[6] = {0,1,2,0,2,3};
};

class Plot {

public:

	Plot( const TopologyBase& topology , bool with_clip=true );

	void build();

  void compute_center();
  void transform_center(const mat4& m);
	void transform( const mat4& m , bool centered);

	const ClipPlane& clip() const { return *clip_.get(); }
	      ClipPlane& clip()       { return *clip_.get(); }

  const mat4& model_matrix() const { return model_matrix_; }
  const vec3& center() const { return center_; }
  float length_scale() const { return length_scale_; }

  index_t nb_vao() const { return vao_.size(); }
  VertexAttributeObject& vao( index_t k ) { return *vao_[k].get(); }

  const VertexAttributeObject& active_vao() const { return *active_vao_; }
  VertexAttributeObject& active_vao() { return *active_vao_; }
	void set_active( index_t k ) { active_vao_ = vao_[k].get(); }

	const std::vector<std::string>& vao_labels() const { return vao_labels_; }

	bool& hidden() { return hidden_; }

private:

	const TopologyBase& topology_;
	VertexAttributeObject* active_vao_;
	mat4 model_matrix_;

	std::shared_ptr<ClipPlane> clip_;
	vec3 center_;
	mat4 center_translation_;
	mat4 inverse_center_translation_;

  std::vector< std::shared_ptr<VertexAttributeObject> > vao_;

  float length_scale_;

	std::vector<std::string> vao_labels_;
	bool hidden_;
};

} // graphics

} // avro

#endif
