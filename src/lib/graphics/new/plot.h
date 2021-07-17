#ifndef AVRO_LIB_GRAPHICS_PLOT_H_
#define AVRO_LIB_GRAPHICS_PLOT_H_

#include "graphics/math.h"
#include "graphics/new/primitives.h"
#include "graphics/new/vao.h"

#include "mesh/topology.h"

#include <vector>

namespace avro
{

class TopologyBase;

namespace graphics
{


class ClipPlane {
public:
	ClipPlane() :
		visible_(false),
		modifying_(false),
		direction_(1),
		style_(0)
	{
		u_ = {1,0,0};
		v_ = {0,1,0};
		coordinates_.resize( 12 , 0.0 );
		indices_ = {0,1,2,0,2,3};
		transform_matrix_ = glm::identity();
	}

	~ClipPlane();

	bool visible() const   { return visible_; }
	bool modifying() const { return modifying_; }

	bool& visible()   { return visible_; }
	bool& modifying() { return modifying_; }

	void initialize( const vec3& c , float length ) {
		center_ = c;
		length_ = length;
		compute();
	}

	void compute() {
		float L = length_;
		for (coord_t d = 0; d < 3; d++) {
			coordinates_[3*0+d] = center_(d) + L*u_(d) - L*v_(d);
			coordinates_[3*1+d] = center_(d) + L*u_(d) + L*v_(d);
			coordinates_[3*2+d] = center_(d) - L*u_(d) + L*v_(d);
			coordinates_[3*3+d] = center_(d) - L*u_(d) - L*v_(d);
		}
		write();
	}

	void write();
	void draw(const mat4& view , const mat4& projection) const;

	void transform( const mat4& m , bool centered ) {
		if (centered) {
			mat4 Tc = glm::translate( glm::identity() , center_ );
			mat4 Tc_inv = glm::inverse( Tc );
			model_matrix_ = Tc * m * Tc_inv * model_matrix_;
			if (modifying_) transform_matrix_ = Tc * m * Tc_inv * transform_matrix_;
		}
		else {
			model_matrix_ = m * model_matrix_;
			if (modifying_) transform_matrix_ = m * transform_matrix_;
		}
	}

	void get( vec3& c , vec3& n ) const {
		c = glm::to_vec3( transform_matrix_ * glm::to_vec4(center_,1.0) );
		vec3 u = glm::to_vec3( transform_matrix_ * glm::to_vec4(u_,0.0) );
		vec3 v = glm::to_vec3( transform_matrix_ * glm::to_vec4(v_,0.0) );
		if (direction_ < 0) n = glm::cross(v,u);
		else n = glm::cross(u,v);
	}

	int& style() { return style_; }
	int style() const { return style_; }

private:
	bool visible_;
	bool modifying_;
	vec3 center_;
	float length_;
	vec3 u_; // easier to maintain tangent vectors instead of normal
	vec3 v_;
	mat4 model_matrix_; // transformation used for rendering
	mat4 transform_matrix_; // transformation of the clipping plane in physical space
	short direction_;
	int style_; // 0 for none, 1 for pixel, 2 for primitive

	gl_index point_buffer_;
	gl_index index_buffer_;
	gl_index vertex_array_;

	std::vector<gl_float> coordinates_;
	std::vector<gl_index> indices_;
};

class Plot {

public:

	Plot( const TopologyBase& topology ) :
		topology_(topology),
    length_scale_(1.0)
	{
    model_matrix_ = glm::identity();
  }

	void build() {

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

  void compute_center() {
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
		center_.print();

    length_scale_ = -1;
    for (coord_t d = 0; d < 3; d++) {
      center_(d) *= (3./coordinates.size());
      float L = (xmax(d) - xmin(d));
      if (L > length_scale_) length_scale_ = L;
    }
		center_.print();

    // store the center translation and its inverse
    center_translation_ =  glm::translate( glm::identity() , center_ );
    inverse_center_translation_ = glm::inverse(center_translation_);
  }

  void transform_center(const mat4& m) {
    center_translation_ = m * center_translation_;
    for (coord_t d = 0; d < 3; d++)
      center_(d) = center_translation_(d,3);
    inverse_center_translation_ = glm::inverse(center_translation_);
  }

	void transform( const mat4& m , bool centered) {
		// m is the transformation relative to the trackball
		if (!clip_.modifying()) {
			if (centered) {
				model_matrix_ = center_translation_ * m * inverse_center_translation_ * model_matrix_;
			}
			else {
        model_matrix_ = m * model_matrix_;
      }

			transform_clip(m,false);
		}
		else {
			transform_clip(m,true);
		}
  }

	void transform_clip( const mat4& m , bool centered ) {
		clip_.transform(m,centered);
	}

	const ClipPlane& clip() const { return clip_; }
	      ClipPlane& clip()       { return clip_; }


  const mat4& model_matrix() const { return model_matrix_; }
  const vec3& center() const { return center_; }
  float length_scale() const { return length_scale_; }

  index_t nb_vao() const { return vao_.size(); }
  VertexAttributeObject& vao( index_t k ) { return *vao_[k].get(); }

  void set_triangle_shader( ShaderProgram* program ) { triangle_shader_ = program; }
  void set_edge_shader( ShaderProgram* program ) { edge_shader_ = program; }
  void set_point_shader( ShaderProgram* program ) { point_shader_ = program; }

  ShaderProgram* triangle_shader() { return triangle_shader_; }
  ShaderProgram* edge_shader() { return edge_shader_; }
  ShaderProgram* point_shader() { return point_shader_; }

  const VertexAttributeObject& active_vao() const { return *active_vao_; }
  VertexAttributeObject& active_vao() { return *active_vao_; }
	void set_active( index_t k ) { active_vao_ = vao_[k].get(); }

	const std::vector<std::string>& vao_labels() const { return vao_labels_; }

	bool& hidden() { return hidden_; }

private:

	const TopologyBase& topology_;
	VertexAttributeObject* active_vao_;
	mat4 model_matrix_;

	ClipPlane clip_;
	vec3 center_;
	mat4 center_translation_;
	mat4 inverse_center_translation_;

  std::vector< std::shared_ptr<VertexAttributeObject> > vao_;

  ShaderProgram* triangle_shader_;
  ShaderProgram* edge_shader_;
  ShaderProgram* point_shader_;
  float length_scale_;

	std::vector<std::string> vao_labels_;
	bool hidden_;
};

} // graphics

} // avro

#endif
