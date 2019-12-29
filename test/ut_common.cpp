#include "common/process.h"
#include "common/tools.h"

//#include "graphics/plotter.h"

typedef luma::real_t REAL;
#include <tetgen1.5.0/predicates.h>

#include <triangle/predicates.h>

using namespace luma;

void
ut_pre(int argc, char** argv)
{
  UNUSED(argc);
  UNUSED(argv);

  int plot_on = 0;
  if (argc>1) plot_on = atoi( argv[1] );

  printInfo();

  // initialize the distributed/shared memory task/thread managers
  #ifdef luma_WITH_MPI
  ProcessMPI::initialize();
  #endif
  ProcessCPU::initialize();
  ProcessGPU::initialize();

  // initialize the predicates
  exactinit(1,0,0,10,10,10);
  exactinit();

  // create a graphics plotter instance and save globally
  //Server server;
  //Plotter* plotter = NULL;
  //plotter = new Plotter(&server);
  //__plotter__ = (void*) plotter;

  if (!plot_on)
  {
    //plotter->on() = false;
  }
}

void
ut_post()
{
  //GET_PLOTTER(plotter);
  //if (plotter==NULL) return;
  //delete plotter;
}
