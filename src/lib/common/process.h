//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_COMMON_PROCESS_H_
#define avro_COMMON_PROCESS_H_

#include "common/thread.h"
#include "common/types.h"

#include <memory>

namespace avro
{

namespace ProcessCPU
{

  void initialize();
  void terminate();

  index_t maximum_concurrent_threads();
  void run_threads(ThreadGroup& threads );

  void enterCriticalSection();
  void leaveCriticalSection();

  index_t nb_cores();
  void set_thread_manager( std::shared_ptr<ThreadManager> manager );
  bool is_running_threads();

  void enable_FPE(bool flag);
  bool FPE_enabled();

  void enable_multithreading(bool flag);
  bool multithreading_enabled();

  void set_max_threads(index_t num_threads);
  index_t max_threads();

} // ProcessCPU

namespace ProcessGPU
{

  void initialize();
  void terminate();

  index_t maximum_concurrent_threads();
  template<typename type> void add_kernel_arg( type* arg , index_t n );
  template<typename type> void retrieve_kernel_value( index_t k , type* value , index_t n );
  void run_threads(std::string& src , std::string& name , index_t n );
  void synchronize();

  void enterCriticalSection();
  void leaveCriticalSection();

  index_t nb_cores();
  void set_thread_manager( std::shared_ptr<ThreadManager> manager );
  bool is_running_threads();

  void enable_FPE(bool flag);
  bool FPE_enabled();

  void enable_multithreading(bool flag);
  bool multithreading_enabled();

  void set_max_threads(index_t num_threads);
  index_t max_threads();

} // ProcessGPU

namespace ProcessMPI
{

  void initialize();
  void terminate();

  index_t nb_processes();
  index_t rank();

  void barrier();

  //  functions
  void main_begin( Thread& task );
  void main_end( Thread& task );

  // worker functions
  void worker_do( Thread& task );

} // ProcessMPI

} // avro

#endif
