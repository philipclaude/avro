//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_COMMON_PARALLEL_FOR_H_
#define avro_COMMON_PARALLEL_FOR_H_

#include "common/process.h"
#include "common/thread.h"
#include "avro_types.h"

#include <algorithm>

namespace avro
{

template<class T>
class ParallelForMemberCallback
{
public:

  // function pointer to member callback
  typedef void (T::*fptr)(index_t);

  // store the object and a function pointer
  // to the callback each thread will run
  ParallelForMemberCallback(T* object, fptr f) :
    object_(object), f_(f)
  {}

  // call the member function
  void operator() ( index_t i )
    { (*object_.*f_)(i); }

  private:
    T* object_;
    fptr f_;
};

// helper function to create ParallelForMemberCallback
template<class T> ParallelForMemberCallback<T>
parallel_for_member_callback(T* obj, void (T::* fun)(index_t))
  { return ParallelForMemberCallback<T>(obj,fun); }

template<class T>
class TaskMemberCallback
{
public:
  typedef void (T::*fptr)(index_t);

  TaskMemberCallback( T* object , fptr f ) :
    object_(object), f_(f)
  {}

  void operator() (index_t k)
    { (*object_.*f_)(k); }

private:
  T* object_;
  fptr f_;
};

template<class T> TaskMemberCallback<T>
task_member( T* obj , void (T::* fun)(index_t) )
  { return TaskMemberCallback<T>(obj,fun); }

namespace ProcessCPU
{
#include "common/parallel_for.hpp"
}

namespace ProcessGPU
{

template<typename type>
inline void
parallel_add_kernel_arg( type* arg , index_t n = 1 )
{
  add_kernel_arg(arg,n);
}

template<typename type>
inline void
parallel_retrieve_value( index_t k , type* value , index_t N )
{
  retrieve_kernel_value(k,value,N);
}

inline void
parallel_for( std::string& source , std::string& name , index_t n )
{
  run_threads(source,name,n);
}

} // ProcessGPU

namespace ProcessMPI
{

template<class Func>
inline void
main_begin( const Func& func )
{
  Thread_var task = std::make_shared<Task<Func>>(func);
  main_begin(*task);
}

template<class Func>
inline void
main_end( const Func& func )
{
  Thread_var task = std::make_shared<Task<Func>>(func);
  main_end(*task);
}

template<class Func>
inline void
worker_do( const Func& func )
{
  Thread_var task = std::make_shared<Task<Func>>(func);
  worker_do(*task);
}

template<typename T>
inline void
send( int dest , int tag , const T& value )
{
  send(dest,tag,value);
}

template<typename T>
inline void
send( int dest , int tag , T* values )
{
  send(dest,tag,values);
}

template<typename T>
inline int
receive( int source , int tag , T& value )
{
  return receive(source,tag,value);
}

template<typename T>
inline int
receive( int source , int tag , T* values , int n )
{
  return receive(source,tag,values,n);
}

}

} // avro

#endif
