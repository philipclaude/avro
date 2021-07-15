#ifndef AVRO_LIB_GRAPHICS_PLOT_H_
#define AVRO_LIB_GRAPHICS_PLOT_H_

#include "graphics/math.h"
#include "graphics/new/vao.h"

#include "mesh/topology.h"

#include <vector>

namespace avro
{

class TopologyBase;

namespace graphics
{

class Plot {

private:
	class ClippingPlane {
	public:
		ClippingPlane() {}

		bool visible() const { return visible_; }
		bool modifying() const { return modifying_; }

		void draw();

	private:
		bool visible_;
		bool modifying_;
		vec3 center_;
		vec3 normal_;
		vec3 model_matrix_;

		gl_index buffer_;
	};

public:

	Plot( const TopologyBase& topology ) :
		topology_(topology)
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
    }

		// compute the center

		// store the center translation

		// store the inverse center translation
	}

	void transform( const mat4& m , bool centered) {
		// m is the transformation relative to the trackball
		if (!clip_.modifying()) {
			if (centered) {
				model_matrix_ = inverse_center_translation_ * m * center_translation_ * model_matrix_;
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
		// translate clip to origin

		// apply transformation

		// translate backwards

		// compound the total transformation to the model matrix
	}

	void draw() {
		active_vao_->draw();

		if (clip_.visible()) {
			// draw clipping plane
		}
	}

  const mat4& model_matrix() const { return model_matrix_; }

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

private:

	const TopologyBase& topology_;
	VertexAttributeObject* active_vao_;
	mat4 model_matrix_;

	ClippingPlane clip_;
	vec3 center_;
	mat4 center_translation_;
	mat4 inverse_center_translation_;

  std::vector< std::shared_ptr<VertexAttributeObject> > vao_;

  ShaderProgram* triangle_shader_;
  ShaderProgram* edge_shader_;
  ShaderProgram* point_shader_;
};

} // graphics

} // avro

#endif
