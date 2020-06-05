//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "graphics/controls.h"
#include "graphics/scene.h"

namespace avro
{

namespace graphics
{

void
SceneGraph::update_matrices( const Controls& controls )
{
  normal_matrix_ = controls.normal();
  mvp_matrix_ = controls.model_view_projection();
}

} // graphics

} // avro
