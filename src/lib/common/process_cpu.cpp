// avro: Adaptive Voronoi Remesher
// Copyright 2017-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "common/process.h"
#include "common/thread.h"

#ifdef avro_CPU_THREAD_MANAGER_PTHREAD
#include <pthread.h>
#endif

#ifdef avro_CPU_THREAD_MANAGER_CPP
#include <thread>
#endif

#ifdef avro_CPU_THREAD_MANAGER_OPENMP
#include <omp.h>
#endif

#ifdef avro_CPU_THREAD_MANAGER_EMP
#include <emp.h>
#endif

#include <memory>
#include <unistd.h>

avro_THREAD_LOCAL avro::Thread* avro_current_thread_ = 0;

namespace avro
{

void
Thread::setCurrent(Thread* thread)
  { avro_current_thread_ = thread; }

namespace ProcessCPU
{
  std::shared_ptr<ThreadManager> thread_manager_;
  int running_threads_invocations_ = 0;

  #ifdef avro_CPU_THREAD_MANAGER_PTHREAD
  class PThreadManager : public ThreadManager
  {
  public:
    PThreadManager()
    {
      pthread_mutex_init(&mutex_,0);
      pthread_attr_init(&attr_);
      pthread_attr_setdetachstate(&attr_,PTHREAD_CREATE_JOINABLE);
    }

    index_t maxConcurrentThreads() { return ProcessCPU::nb_cores(); }

    void enterCriticalSection()
    {
      pthread_mutex_lock(&mutex_);
    }

    void leaveCriticalSection()
    {
      pthread_mutex_unlock(&mutex_);
    }

    static void* run_thread(void* thread_in)
    {
      Thread* thread = reinterpret_cast<Thread*>(thread_in);
      setCurrentThread(thread);
      thread->run();
      return NULL;
    }

    void runConcurrentThreads( ThreadGroup& threads0 , index_t max_threads )
    {
      printf("running with pthreads...\n");
      threads_.resize( threads0.size() );
      for (index_t i=0;i<threads0.size();i++)
      {
        Thread* T = threads0[i].get();
        setThreadId(T,i);
        pthread_create(&threads_[i],&attr_,&run_thread,T);
      }

      for (index_t i=0;i<threads_.size();i++)
        pthread_join(threads_[i],NULL);
    }

  virtual ~PThreadManager()
  {
    pthread_attr_destroy(&attr_);
    pthread_mutex_destroy(&mutex_);
  }

  private:
    pthread_attr_t attr_;
    pthread_mutex_t mutex_;
    std::vector<pthread_t> threads_;

  };
  #endif

  #ifdef avro_CPU_THREAD_MANAGER_CPP
  class CppThreadManager : public ThreadManager
  {
  public:
    CppThreadManager() {}

    index_t maxConcurrentThreads()
      { return std::thread::hardware_concurrency(); }

    void enterCriticalSection()
      { avro_implement; }

    void leaveCriticalSection()
      { avro_implement; }

    static void* run_thread( void* thread_in )
    {
      Thread* thread = reinterpret_cast<Thread*>(thread_in);
      setCurrentThread(thread);
      thread->run();
      return NULL;
    }

    void runConcurrentThreads( ThreadGroup& threads0 , index_t max_threads )
    {
      threads_.resize( threads0.size() );
      for (index_t i=0;i<threads0.size();i++)
      {
        Thread* T = threads0[i].get();
        setThreadId(T,i);
        threads_[i] = std::thread( &run_thread , T );
      }

      for (index_t i=0;i<threads_.size();i++)
        threads_[i].join();
    }
  private:
    std::vector<std::thread> threads_;

  };
  #endif

  #ifdef avro_CPU_THREAD_MANAGER_OPENMP
  class OpenMPThreadManager : public ThreadManager
  {
  public:
    OpenMPThreadManager()
    {
      omp_init_lock(&lock_);
    }

    index_t maxConcurrentThreads()
    {
      return omp_get_max_threads();
    }
    void enterCriticalSection()
    {
      omp_set_lock(&lock_);
    }

    void leaveCriticalSection()
    {
      omp_unset_lock(&lock_);
    }

    void runConcurrentThreads( ThreadGroup& threads , index_t max_threads )
    {
      #pragma omp parallel for schedule(dynamic) // pcaplan is schedule dynamic not supported with intel?
      for (index_t i=0;i<threads.size();i++)
      {
        setThreadId(threads[i].get(),i);
        setCurrentThread(threads[i].get());
        threads[i]->run();
      }
    }

    virtual ~OpenMPThreadManager()
    {
      omp_destroy_lock(&lock_);
    }

  private:
    omp_lock_t lock_;

  };
  #endif

  #ifdef avro_CPU_THREAD_MANAGER_EMP
  class EMPThreadManager : public ThreadManager
  {
  public:
    EMPThreadManager()
    {
      EMP_Init(NULL);
      lock_ = EMP_LockCreate();
    }

    index_t maxConcurrentThreads() { return ProcessCPU::nb_cores(); }

    void enterCriticalSection()
    {
      EMP_LockSet(lock_);
    }

    void leaveCriticalSection()
    {
      EMP_LockRelease(lock_);
    }

    static void run_thread( void* thread_in )
    {
      Thread* thread = reinterpret_cast<Thread*>(thread_in);
      setCurrentThread(thread);
      thread->run();
    }

    void runConcurrentThreads( ThreadGroup& threads0 , index_t max_threads )
    {

      // setup the threads
      threads_.resize( threads0.size() );
      for (index_t i=0;i<threads0.size();i++)
      {
        Thread* T = threads0[i].get();
        setThreadId(T,i);
        threads_[i] = EMP_ThreadCreate( &run_thread , T );
      }

      // synchronize
      for (index_t i=0;i<threads_.size();i++)
      {
        if (threads_[i]!=NULL) EMP_ThreadWait(threads_[i]);
      }

      // cleanup
      for (index_t i=0;i<threads_.size();i++)
      {
        if (threads_[i]!=NULL) EMP_ThreadDestroy(threads_[i]);
      }
    }

    virtual ~EMPThreadManager()
    {
      EMP_LockDestroy(lock_);
      EMP_Done(NULL);
    }

  private:
    void* lock_;
    std::vector<void*> threads_;
  };
  #endif

void set_thread_manager( std::shared_ptr<ThreadManager> manager )
{
  thread_manager_ = manager;
}

void initialize()
{

  // set the appropriate thread manager
  #if defined(avro_CPU_THREAD_MANAGER_OPENMP)
    printf("enabling openmp threads\n");
    set_thread_manager( std::make_shared<OpenMPThreadManager>() );
  #elif defined(avro_CPU_THREAD_MANAGER_CPP)
    printf("enabling cpp threads\n");
    set_thread_manager(  std::make_shared<CppThreadManager>() );
  #elif defined(avro_CPU_THREAD_MANAGER_PTHREAD)
    printf("enabling pthreads\n");
    set_thread_manager( std::make_shared<PThreadManager>() );
  #elif defined(avro_CPU_THREAD_MANAGER_EMP)
    printf("enabling emp threads\n");
    set_thread_manager(  std::make_shared<EMPThreadManager>() );
  #else
    printf("only enabling serial computations\n");
    set_thread_manager(  std::make_shared<SerialThreadManager>() );
  #endif

}

index_t
nb_cores()
{
  return index_t(sysconf(_SC_NPROCESSORS_ONLN));
}

index_t
maximum_concurrent_threads()
  { return thread_manager_->maxConcurrentThreads(); }

bool
is_running_threads()
  { return running_threads_invocations_>0; }

void
run_threads( ThreadGroup& threads )
{
  running_threads_invocations_++;
  thread_manager_->runThreads(threads);
  running_threads_invocations_--;
}

void
terminate()
{}

} // ProcessCPU

} // avro
