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
	{}

	void build() {
		vao_.build(topology_);

		// compute the center

		// store the center translation

		// store the inverse center translation
	}

	void transform( const mat4& m , bool centered) {
		// m is the transformation relative to the trackball
		if (!clip_.modifying()) {
			if (centered) {
				mat4 m1 = inverse_center_translation_ * m * center_translation_;
				vao_.apply_transformation(m1);
			}
			else
				vao_.apply_transformation(m);

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
		vao_.draw();

		if (clip_.visible()) {
			// draw clipping plane
		}
	}

private:

	const TopologyBase& topology_;
	VertexAttributeObject vao_;
	mat4 model_matrix_;

	ClippingPlane clip_;
	vec3 center_;
	mat4 center_translation_;
	mat4 inverse_center_translation_;
};

} // graphics

} // avro

#endif
