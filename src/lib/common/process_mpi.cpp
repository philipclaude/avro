//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "common/process.h"

#ifdef avro_WITH_MPI
#include <mpi.h>
#include "common/mpi.hpp"
#endif

namespace avro
{

namespace ProcessMPI
{

// task manager used to control distributed memory applications
class TaskManager
{
public:
  virtual index_t nb_processes() = 0;
  virtual index_t rank() = 0;

  virtual void shape_begin( Thread& thread ) = 0;
  virtual void worker_do( Thread& thread ) = 0;
  virtual void shape_end( Thread& thread) = 0;
  virtual void barrier() = 0;

protected:
  virtual ~TaskManager() {}

};

std::shared_ptr<TaskManager> task_manager_;

// overrides the task distribution and runs everything on the same processor
class SerialTaskManager : public TaskManager
{
public:
  index_t nb_processes() { return 1; }
  index_t rank() { return 0; }

  void shape_begin( Thread& thread ) { thread.run(); }
  void worker_do( Thread& thread ) { thread.run(); }
  void shape_end( Thread& thread) { thread.run(); }
  void barrier() {}

};

#ifdef avro_WITH_MPI
class MPITaskManager : public TaskManager
{
public:
  MPITaskManager() :
    comm_(mpi::communicator::world),
    instance_(0,NULL)
  {
  }

  index_t nb_processes()
  {
    return mpi::size(comm_);
  }

  index_t rank()
  {
    return mpi::rank();
  }

  void shape_begin( Thread& thread )
  {
    if (rank()!=0) return;
    thread.run();
  }

  void worker_do( Thread& thread )
  {
    if (rank()==0) return;
    thread.run(rank());
  }

  void shape_end( Thread& thread)
  {
    mpi::barrier(comm_);
    if (rank()!=0) return;
    thread.run();
  }

  void barrier()
  {
    printf("waiting..\n");
    mpi::barrier(comm_);
  }

private:
  mpi::communicator comm_;
  mpi::instance instance_;

};
#endif

index_t nb_processes()
  { return task_manager_->nb_processes(); }

index_t rank()
  { return task_manager_->rank(); }

void initialize()
{
  #ifdef avro_WITH_MPI
  task_manager_ = std::make_shared<MPITaskManager>();
  #else
  task_manager_ = std::make_shared<SerialTaskManager>();
  #endif
}

void terminate()
{
  //delete task_manager_;
}

void
shape_begin( Thread& thread )
{
  task_manager_->shape_begin(thread);
}

void
worker_do( Thread& thread )
{
  task_manager_->worker_do(thread);
}

void
shape_end( Thread& thread )
{
  task_manager_->shape_end(thread);
}

void
barrier()
{
  task_manager_->barrier();
}

} // ProcessMPI

} // avro
