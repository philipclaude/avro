#ifndef AVRO_LIB_GRAPHICS_PLOT_H_
#define AVRO_LIB_GRAPHICS_PLOT_H_

#include "graphics/math.h"
#include "graphics/new/vao.h"

#include <vector>

namespace avro
{

class TopologyBase;

namespace graphics
{

class VertexAttributeObject;

class Plot {

public:
	Plot( const TopologyBase& topology );

	void add( const TopologyBase& topology );
	/* {
		coord_t number = topology.number();
		coord_t order  = topology.shape().order();
		std::shared_ptr<VertexAttributeObject> vao = std::make_shared<inverse>(number,order);
		vao->build(topology);
		vertex_array_objects_.push_back(vao);
	}*/

	void apply_transformation( const mat4& m , int picked=-1 ) {
    if (picked < 0) {
      for (index_t k = 0; k < vao_.size(); k++) {
        vao_[k]->apply_transformation(m);
      }
    }
    else {
      vao_[picked]->apply_transformation(m);
    }
  }

private:
	std::vector< std::shared_ptr<VertexAttributeObject> > vao_;
};

} // graphics

} // avro

#endif
