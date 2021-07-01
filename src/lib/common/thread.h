//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_COMMON_THREAD_H_
#define avro_COMMON_THREAD_H_

#include "common/error.h"
#include "avro_types.h"

#include <memory>
#include <string>
#include <vector>

#define avro_THREAD_LOCAL __thread

namespace avro
{

class Thread
{
public:
  Thread() : id_(0)
  {}

  virtual void run() = 0;
  virtual void run( const index_t k ) {}

  index_t id() const { return id_; }

  static Thread* current();

protected:
  virtual ~Thread() {}

private:
  void set_id( index_t _id ) { id_ = _id; }
  static void setCurrent( Thread* _thread );
  index_t id_;

  friend class ThreadManager;
};

typedef std::shared_ptr<Thread> Thread_var;
typedef std::vector<Thread_var> ThreadGroup;

class ThreadManager
{
public:

  virtual std::string name() const = 0;

  // run threads from a group of threads defining a common task
  virtual void
  runThreads( ThreadGroup& threads )
  {
    index_t max_threads = maxConcurrentThreads();
    runConcurrentThreads(threads,max_threads);
  }

  // run threads for a kernel stored as a string in source
  virtual void
  runThreads( std::string& source , std::string& name , index_t n )
    { runConcurrentThreads(source,name,n); }

  // default thread runner from kernel source
  virtual void
  runConcurrentThreads( std::string& source , std::string& name , index_t n )
    { avro_assert_not_reached; }

  // default argument adder for kernels compiled from a string
  virtual void
  addArgument( void* arg , size_t sz , index_t n )
    { avro_assert_not_reached ; }

  // value retrieval upon terimnination of all threads
  virtual void
  retrieveValue( index_t k , index_t n , size_t sz , void* arg )
    { avro_assert_not_reached; }

  // thread synchornization
  virtual void synchronize()
    { avro_assert_not_reached; }

  // these must be defined by the derived thread managers
  virtual index_t maxConcurrentThreads() = 0;
  virtual void enterCriticalSection() = 0;
  virtual void leaveCriticalSection() = 0;

  virtual ~ThreadManager() {}

protected:
  virtual void
  runConcurrentThreads( ThreadGroup& threads , index_t max_threads ) = 0;
  static void setThreadId( Thread* thread , index_t _id )
    { thread->set_id(_id); }
  static void setCurrentThread( Thread* thread )
    { Thread::setCurrent(thread); }

};

class SerialThreadManager : public ThreadManager
{
public:
  SerialThreadManager() {}
  void runThreads( ThreadGroup& threads )
  {
    for (index_t i=0;i<threads.size();i++)
      threads[i]->run();
  }

  std::string name() const { return "serial"; }

  index_t maxConcurrentThreads() { return 1; }

  // nothing to do upon entry/exit of critical sections
  void enterCriticalSection() {}
  void leaveCriticalSection() {}

  void runConcurrentThreads( ThreadGroup& threads , index_t max_threads )
  { avro_implement; }
};

template<class Func>
class ParallelForThread : public Thread
{
public:
  ParallelForThread( const Func& _func , index_t _from , index_t _to , index_t _step = 1 ) :
    func_(_func),
    from_(_from),
    to_(_to),
    step_(_step)
  {}

  virtual void run()
  {
    for (index_t k=from_;k<to_;k+=step_)
        const_cast<Func&>(func_)(k);
  }

  virtual ~ParallelForThread() {}

private:
  const Func& func_;
  index_t from_;
  index_t to_;
  index_t step_;
};

template<class Func>
class Task : public Thread
{
public:
  Task( const Func& _func ) :
    func_(_func)
  {}

  virtual void run()
  {
    const_cast<Func&>(func_)(0);
  }

  virtual void run( const index_t k )
  {
    const_cast<Func&>(func_)(k);
  }

  virtual ~Task() {}

private:
  const Func& func_;
};

} // avro

#endif
