#ifndef LUNA_LIB_GEOMETRY_EGADS_CONTEXT_H_
#define LUNA_LIB_GEOMETRY_EGADS_CONTEXT_H_

#include "geometry/egads/egads_fwd.h"

namespace luna
{

namespace EGADS
{

class Context
{
public:
  Context();
  Context( ego* context );
  ~Context();

  ego* get();
  const ego* get() const;

  void print() const;

private:
  ego* context_;
  bool mine_;
};

} // EGADS

} // luna

#endif
