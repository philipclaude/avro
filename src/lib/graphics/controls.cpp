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

#if USE_GLM
#define GLM_ENABLE_EXPERIMENTAL
#include "glm/gtx/string_cast.hpp"
#endif

#include "json/json.hpp"

#include <fstream>
#include <iomanip>

namespace avro
{

namespace graphics
{

#define NEW_CONTROLS 1

Controls::Controls( float fov , int width , int height , float znear , float zfar ) :
  width_(width),
  height_(height),
  offleft_(0),
  offtop_(0),
  scale_(1.0),
  enabled_(true)
{
  perspective_ = glm::perspective(fov,(float)width_/(float)height_,znear,zfar);
  perspective_ = perspective_ * glm::lookAt( eye_ , center_ , up_ );

  model_view_ = glm::identity();
  ui_matrix_  = model_view_;

  startx_ = -1;
  starty_ = -1;
  cursorx_ = -1;
  cursory_ = -1;

  modifier = 0;
  dragging = false;
}

bool
Controls::update()
{
  ui_matrix_  = model_view_;
  model_view_ = glm::identity();

  bool updated = false;
  if (!enabled_) return updated;

  rotation_ = glm::identity();
  translation_ = glm::identity();

  // mouse-movement
  if ( dragging )
  {
    if (modifier==4) // control is down (rotate)
    {
      float anglex =  (starty_ -cursory_)/100.0;
      float angley = -(startx_ -cursorx_)/100.0;

      if (anglex!=0.0f || angley!=0.0f)
      {
        model_view_ = glm::rotate( model_view_ , anglex , {1,0,0} );
        model_view_ = glm::rotate( model_view_ , angley , {0,1,0} );
        updated = true;

        rotation_ = model_view_;
      }
    }
    else if (modifier==3) // alt-shift is down (rotate)
    {
      float anglex =  (starty_ -cursory_)/10.0;
      float angley = -(startx_ -cursorx_)/10.0;

      if (anglex!=0.0f || angley!=0.0f)
      {
        model_view_ = glm::rotate( model_view_ , anglex , {1,0,0} );
        model_view_ = glm::rotate( model_view_ , angley , {0,1,0} );
        updated = true;

        rotation_ = model_view_;
      }
    }
    else if (modifier==2) // alt is down (spin)
    {
      float xf = 1.0f*startx_ - width_*1.0f/2.0;
      float yf = 1.0f*starty_ - height_*1.0f/2.0;

      if (xf!=0.0 || yf!=0.0)
      {
        float theta = std::atan2(yf,xf);
        xf = 1.0f*cursorx_ - width_*1.0f/2.;
        yf = 1.0f*cursory_ - height_*1.0f/2.;

        if (xf!=0.0 || yf!=0.0)
        {
          float dtheta = std::atan2(yf,xf) - theta;
          if (dtheta < 1.5708)
          {
            float anglez = 128*dtheta/3.141592654;
            model_view_  = glm::rotate( model_view_ , anglez , {0,0,1} );
            updated = true;

            rotation_ = model_view_;
          }
        }

      }
    }
    else if (modifier==1) // shift is down (zoom)
    {
      if (cursory_!=starty_)
      {
        float scale = std::exp( (cursory_ - starty_)/512.0 );
        model_view_ = glm::scale( model_view_ , {scale,scale,scale} );
        scale_ *= scale;
        updated = true;
      }

    }
    else // no modifier (translate)
    {
      float transx = (cursorx_ - startx_) / 256.0;
      float transy = (cursory_ - starty_) / 256.0;
      if (transx!=0.0 || transy!=0.0)
      {
        model_view_ = glm::translate( model_view_ , {transx,transy,0.0f} );
        updated = true;
        translation_ = model_view_;
      }
    }

    startx_ = cursorx_;
    starty_ = cursory_;
  }

  return updated;
}

void
Controls::calculate_view()
{
  // accumulate the complete transformation
  model_view_ = model_view_*ui_matrix_;

  // construct normal matrix from model-view
  normal_ = glm::inverse( model_view_ );
  normal_ = glm::transpose( glm::scale( normal_ , {scale_,scale_,scale_} ) );

  // construct model-view-projection from original perspective and model-view
  model_view_projection_ = glm::scale( perspective_ , {1.,1.,1./scale_} ) * model_view_;
}

void
Controls::mouse_up()
{
  dragging = false;
  modifier = 0;
}

void
Controls::mouse_down(int button, int action, int mods,int x,int y)
{
  startx_ = x;
  starty_ = y;

  startx_ -= offleft_ +1;
  starty_  = height_ - starty_ + offtop_ +1;

  dragging = true;

  modifier = 0;
  #ifdef AVRO_WITH_GL
  if (mods == GLFW_MOD_SHIFT) modifier = 1;
  if (mods == GLFW_MOD_ALT) modifier = 2;
  if (mods == GLFW_MOD_CONTROL) modifier = 4;
  #endif
}

void
Controls::mouse_move(int x, int y)
{
  cursorx_ = x;
  cursory_ = y;

  cursorx_ -= offleft_ + 1;
  cursory_  = height_ - cursory_ + offtop_ + 1;

  #if NEW_CONTROLS
  if (dragging) {

  }
  #endif
}

void
Controls::key_down(int key)
{
  if (key == 77 /*GLFW_KEY_M*/)
  {
    std::cout << glm::to_string(model_view_projection_) << std::endl;
    std::cout << glm::to_string(ui_matrix_) << std::endl;

    nlohmann::json J;
    std::vector<real_t> mvp(16);
    std::vector<real_t> ui(16);
    int k = 0;
    for (coord_t i = 0; i < 4; i++)
    for (coord_t j = 0; j < 4; j++)
    {
      mvp[k] = model_view_projection_[i][j];
      ui[k++]  = ui_matrix_[i][j];
    }
    J["mvp"] = mvp;
    J["ui"] = ui;
    J["scale"] = scale_;
    std::ofstream output("view.json");
    output << std::setw(4) << J << std::endl;
    output.close();

    // TODO write perspective, but this should always be the same

  }
}

void
Controls::load( const std::string& filename )
{
  std::fstream file;
  file.open(filename);
  avro_assert( file.is_open() );
  std::stringstream s;
  s << file.rdbuf();

  nlohmann::json J = nlohmann::json::parse( s.str().c_str() );

  std::vector<real_t> mvp = J.at("mvp");
  std::vector<real_t> ui = J.at("ui");
  scale_ = J.at("scale");

  int k = 0;
  for (coord_t i = 0; i < 4; i++)
  for (coord_t j = 0; j < 4; j++)
  {
    model_view_projection_[i][j] = mvp[k];
    ui_matrix_[i][j] = ui[k++];
  }

  std::cout << glm::to_string(model_view_projection_) << std::endl;
  std::cout << glm::to_string(ui_matrix_) << std::endl;

}

void
Controls::mouse_wheel(double deltax ,double deltay)
{
  if (deltay > 0)
  {
    model_view_ = glm::scale( model_view_ , {1.1,1.1,1.1} );
    scale_ *= 1.1;
  }
  else if (deltay < 0)
  {
    model_view_ = glm::scale( model_view_ , {0.9,0.9,0.9} );
    scale_ *= 0.9;
  }
}

void
Controls::reset()
{
  model_view_ = glm::identity();
  ui_matrix_  = glm::identity();
  scale_ = 1.0;

  startx_ = -1;
  starty_ = -1;
  cursorx_ = -1;
  cursory_ = -1;

  modifier = 0;
  dragging = false;

  calculate_view();
}


} // graphics

} // avro
