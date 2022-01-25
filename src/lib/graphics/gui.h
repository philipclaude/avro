//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef AVRO_LIB_GRAPHICS_GUI_H_
#define AVRO_LIB_GRAPHICS_GUI_H_

#include <imgui/imgui.h>

namespace avro
{

namespace graphics
{

class Window;

class GUI {

public:
  GUI( Window& window );

  void draw();
  void end_draw();

  bool redraw() {
    if (count_ == 0) {
      count_ = 5;
      return false;
    }
    else {
      count_--;
      return true;
    }
  }

private:
  Window& window_;
  //ImGuiIO context_;

  int count_;
};

} // graphics

} // avro

#endif
