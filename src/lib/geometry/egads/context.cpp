#include "geometry/egads/context.h"

#include <egads.h>

namespace luna
{

namespace EGADS
{

Context::Context() :
  mine_(true)
{
	context_ = (ego*) malloc( sizeof(ego) );
	EG_open( context_ );
}

Context::Context( ego* _context ) :
  context_( _context ),
  mine_(false)
{}

ego*
Context::get()
{
	return context_;
}

const ego*
Context::get() const
{
  return context_;
}

void
Context::print() const
{
  int major,minor;
  const char *occrev;
  EG_revision(&major,&minor,&occrev);
  printf("creating geometry context with EGADS v%d.%d and OpenCASCADE %s\n",major,minor,occrev);
}

Context::~Context()
{
  if (mine_)
  {
	   EG_close(*context_);
	    free(context_);
  }
}

} // EGADS

} // luna
