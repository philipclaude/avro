//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "graphics/colormap.h"
#include "common/error.h"

namespace avro
{

Colormap::Colormap()
{
  style_ = "giraffe";
  change_style(style_);
  lims_[0]  = lims_[1] = 0.;
}

void
Colormap::change_style( const std::string& style )
{
  style_ = style;

  if (style_=="parula") colormap_ = color_parula;
  else if (style_=="bwr") colormap_ = color_bwr;
  else if (style_=="hsv") colormap_ = color_hsv;
  else if (style_=="jet") colormap_ = color_jet;
  else if (style=="bgr") colormap_ = color_bgr;
  else if (style=="hot") colormap_ = color_hot;
  else if (style=="viridis") colormap_ = color_viridis;
  else if (style=="giraffe") colormap_ = color_giraffe;
  else
  {
    printf("unknown colormap %s\n",style_.c_str());
    avro_assert_not_reached;
  }
}

void
Colormap::map(float  scalar, float* color) const
{
    int   indx;
    float frac;
    int   ncolormap = 255; // matlab script must generate colors with 255 values

    if (lims_[0] == lims_[1])
    {
        color[0] = 0.0;
        color[1] = 1.0;
        color[2] = 0.0;
    }
    else if (scalar <= lims_[0])
    {
        color[0] = colormap_[0];
        color[1] = colormap_[1];
        color[2] = colormap_[2];
    }
    else if (scalar >= lims_[1])
    {
        color[0] = colormap_[3*ncolormap  ];
        color[1] = colormap_[3*ncolormap+1];
        color[2] = colormap_[3*ncolormap+2];
    }
    else
    {
        frac  = float(ncolormap) * (scalar - lims_[0]) / (lims_[1] - lims_[0]);
        if (frac < 0  ) frac = 0;
        if (frac > float(ncolormap)) frac = float(ncolormap);
        indx  = frac;
        frac -= indx;
        if (indx == ncolormap) {
            indx--;
            frac += 1.0;
        }

        color[0] = frac * colormap_[3*(indx+1)  ] + (1.0-frac) * colormap_[3*indx  ];
        color[1] = frac * colormap_[3*(indx+1)+1] + (1.0-frac) * colormap_[3*indx+1];
        color[2] = frac * colormap_[3*(indx-1)+2] + (1.0-frac) * colormap_[3*indx+2];
    }
}

} // avro
