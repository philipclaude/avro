//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef AVRO_LIB_GRAPHICS_APPLICATION_H_
#define AVRO_LIB_GRAPHICS_APPLICATION_H_

#include "avro_params.h"

#include "graphics/gui.h"
#include "graphics/plot.h"
#include "graphics/vao.h"
#include "graphics/window.h"
#include "graphics/webglpp.h"

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

class ViewerBase {
public:
	virtual void add( const TopologyBase& topology ) = 0;
	virtual void run(bool quit) = 0;
	virtual ~ViewerBase() {};
};

class OpenGL_Application : public ViewerBase {

public:
	OpenGL_Application();

	void run(bool quit = false);
	void add( const TopologyBase& topology );

	const Window& window() const { return window_; }
	Window& window() { return window_; }

private:
	Window window_;
	std::shared_ptr<GUI> gui_;
	std::vector<std::shared_ptr<Plot>> plot_;

};

class WebGL_Application : public ViewerBase {

public:
	void run(bool quit = false);
	void run_thread(bool quit = false);
	void add( const TopologyBase& topology );

private:
	WebGL_Manager manager_;
	std::vector<std::shared_ptr<Plot>> plot_;
};

class Viewer {
public:
	Viewer( bool web=false );

	void run(bool quit = false);
	void add( const TopologyBase& topology );

private:
	std::shared_ptr<ViewerBase> app_;
};

} // graphics

} // avro

#endif
