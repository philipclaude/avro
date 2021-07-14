#ifndef AVRO_LIB_GRAPHICS_APPLICATION_H_
#define AVRO_LIB_GRAPHICS_APPLICATION_H_

#include "avro_params.h"

#include "graphics/new/plot.h"
#include "graphics/new/vao.h"
#include "graphics/new/window.h"

#include "graphics/gl.h"
#include "graphics/math.h"

#include <memory>
#include <vector>

namespace avro
{

class TopologyBase;

namespace graphics
{

class Plot;


class OpenGL_Application {

public:
	OpenGL_Application();

	void run();

	void add( const TopologyBase& topology );

private:
	Window window_;
	std::vector<std::shared_ptr<Plot>> plot_;

};

class WebGL_Application {

	void run();

	void add( const TopologyBase& topology );

private:
	//WebGL_Manager manager_;
	std::vector<std::shared_ptr<Plot>> plot_;
};

} // graphics

} // avro

#endif
