#ifndef LUNA_LIB_GEOMETRY_EGADS_DATA_H_
#define LUNA_LIB_GEOMETRY_EGADS_DATA_H_

#include "geometry/egads/egads_fwd.h"

namespace luna
{

namespace EGADS
{

typedef struct
{
  int object_class, member_type;
  int nb_children, *senses;
  ego reference, *children;
} egoData;

} // EGADS

} // luna

#endif
