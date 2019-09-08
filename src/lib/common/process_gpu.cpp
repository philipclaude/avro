// ursa: Adaptive Voronoi Remesher
// Copyright 2017-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "common/process.h"
#include "common/thread.h"
#include "common/tools.h"

#ifdef URSA_GPU_THREAD_MANAGER_OPENCL
  #if __APPLE__
    #include <OpenCL/opencl.h>
  #else
    #include <CL/cl.h>
#endif

#include <math.h>

#ifdef URSA_GPU_THREAD_MANAGER_OPENCL

static std::string
clErrorCode(int err)
{
  // some common error codes i keep faulting on
  std::string s;
  if (err==-1) s = "CL_DEVICE_NOT_FOUND";
  else if (err==-2) s = "CL_DEVICE_NOT_AVAILABLE";
  else if (err==-3) s = "CL_COMPILER_NOT_AVAILABLE";
  else if (err==-30) s = "CL_INVALID_VALUE";
  else if (err==-44) s = "CL_INVALID_PROGRAM";
  else if (err==-45) s = "CL_INVALID_PROGRAM_EXECUTABLE";
  else if (err==-46) s = "CL_INVALID_KERNEL_NAME";
  else if (err==-47) s = "CL_INVALID_KERNEL_DEFINITION";
  else if (err==-48) s = "CL_INVALID_KERNEL";
  else if (err==-51) s = "CL_INVALID_ARG_SIZE";
  else if (err==-52) s = "CL_INVALID_KERNEL_ARGS";
  return s;
}

#define CL_CHECK(X) { \
   cl_int err = (X); \
   if (err!=CL_SUCCESS) \
   printf("OpenCL error %d (%s) at line %d of %s, running %s\n", err , \
          clErrorCode(err).c_str(), \
          __LINE__,__FILE__,(#X)); \
  }
#endif
#endif

namespace ursa
{

namespace ProcessGPU
{
  std::shared_ptr<ThreadManager> thread_manager_;

  #ifdef URSA_GPU_THREAD_MANAGER_CUDA
  class CUDAThreadManager : public ThreadManager
  {
  public:
    CUDAThreadManager() {}

    index_t maxConcurrentThreads()
      { ursa_implement; }
    void enterCriticalSection()
      { ursa_implement; }
    void leaveCriticalSection()
      { ursa_implement; }
    void runConcurrentThreads( ThreadGroup& threads , index_t max_threads )
      { ursa_implement; }
  };
  #endif

  #ifdef URSA_GPU_THREAD_MANAGER_OPENCL
  class OpenCLThreadManager : public ThreadManager
  {
  public:
    OpenCLThreadManager()
    {

      UNUSED( nb_devices_ );
      UNUSED( event_ );
      UNUSED( info_ );
      UNUSED( gpu_ );

      CL_CHECK( clGetPlatformIDs(1,&platform_,NULL) );

      // initialize the device
      CL_CHECK( clGetDeviceIDs(platform_,CL_DEVICE_TYPE_GPU,1,&device_,NULL) );
      if (device_==NULL)
      {
        printf("falling back on opencl-cpu\n");
        err_ = clGetDeviceIDs(platform_,CL_DEVICE_TYPE_CPU,1,&device_,NULL);
        ursa_assert(device_!=NULL);
      }
      printf ("gpu enabled!\n");
      clGetDeviceInfo(device_,CL_DEVICE_NAME,len_,name_,NULL);
      printf("initializing opencl using the %s with \n \
      \t %lu possible threads and %lu possible compute units.\n",
          name_,maxConcurrentThreads(),nb_units());

      // initialize the context
      properties_[0] = CL_CONTEXT_PLATFORM;
      properties_[1] = (cl_context_properties)platform_;
      properties_[2] = 0;
      context_ = clCreateContext(properties_,1,&device_,NULL,NULL,&err_);
      CL_CHECK(err_);

      // create the queue
      queue_ = clCreateCommandQueue( context_ , device_ , 0 , &err_ );
      CL_CHECK( err_ );

    }

    index_t maxConcurrentThreads()
      { return index_t(CL_DEVICE_MAX_WORK_GROUP_SIZE); }

    index_t nb_units()
      { return index_t(CL_DEVICE_MAX_COMPUTE_UNITS); }

    void enterCriticalSection()
      { ursa_implement; }
    void leaveCriticalSection()
      { ursa_implement; }

    void runConcurrentThreads( std::string& src , std::string& name , index_t n )
    {
      printf("running with opencl\n");

      // clear memory leftover from a previous computation
      for (index_t k=0;k<buffers_.size();k++)
        clReleaseMemObject( buffers_[k] );
      buffers_.clear();

      // get the source for this thread group
      const char* s = src.c_str();
      program_ = clCreateProgramWithSource(context_,1,(const char**)&s,NULL,&err_);
      CL_CHECK(err_);

      if (clBuildProgram(program_,0,NULL,NULL,NULL,NULL)!=CL_SUCCESS)
      {
        char buffer[2048];
        size_t length;
        clGetProgramBuildInfo(program_,device_,CL_PROGRAM_BUILD_LOG,sizeof(buffer),buffer,&length);
        printf("build log:\n%s",buffer);
      }

      kernel_ = clCreateKernel(program_,name.c_str(),&err_);
      CL_CHECK(err_);

      nb_threads_ = n;

    }

    void addArgument( void* arg , size_t sz , index_t n )
    {
      index_t k = buffers_.size();

      buffers_.push_back( clCreateBuffer(context_, CL_MEM_READ_WRITE, n*sz, NULL, &err_ ) );
      CL_CHECK( err_ );

      CL_CHECK( clEnqueueWriteBuffer( queue_, buffers_[k] , CL_TRUE, 0, n*sz , arg, 0, NULL, NULL ) );
      CL_CHECK( clSetKernelArg( kernel_ , k , sizeof(buffers_[k]) , &buffers_[k] ) );

      printf("added argument %lu with size %lu, nb = %lu.\n",k,sz,n);
    }

    void synchronize()
    {
      CL_CHECK( clEnqueueNDRangeKernel( queue_ , kernel_ , 1, NULL , &nb_threads_ , &nb_threads_ , 0 , NULL , NULL ) );
      clFinish(queue_);
    }

    void retrieveValue( index_t k , index_t N , index_t sz , void* value )
    {
      printf("read %lu values of size %lu in argument %lu\n",N,sz,k);
      CL_CHECK( clEnqueueReadBuffer( queue_ , buffers_[k] , CL_TRUE , 0 , N*sz , value ,0 , NULL , NULL ) );
      clFinish(queue_);
    }

    // this just needs to be defined
    void runConcurrentThreads( ThreadGroup& threads , index_t max_threads )
      { ursa_assert_not_reached; }

    virtual ~OpenCLThreadManager()
    {
      for (index_t k=0;k<buffers_.size();k++)
        clReleaseMemObject( buffers_[k] );
      buffers_.clear();
      CL_CHECK( clFinish(queue_) );
      CL_CHECK( clReleaseCommandQueue(queue_) );
      CL_CHECK( clReleaseContext(context_) );
    }

  private:
    cl_platform_id platform_;
    cl_device_id device_;
    cl_context context_;
    cl_context_properties properties_[3];
    cl_uint nb_devices_;
    cl_program program_;
    cl_command_queue queue_;
    std::vector<cl_device_id> devices_;
    cl_event event_;
    cl_kernel kernel_;
    index_t nb_threads_;

    char info_[128];
    char name_[128];
    const index_t len_ = 128;
    cl_int err_;
    bool gpu_;
    std::vector<cl_mem> buffers_;

  };
  #endif

void set_thread_manager( std::shared_ptr<ThreadManager> manager )
{
  thread_manager_ = manager;
}

void
initialize()
{
  #if defined(URSA_GPU_THREAD_MANAGER_CUDA)
    printf("enabling cuda for GPU threading\n");
    set_thread_manager( std::make_shared<CUDAThreadManager>() );
  #elif defined(URSA_GPU_THREAD_MANAGER_OPENCL)
    printf("enabling opencl for GPU threading\n");
    set_thread_manager( std::make_shared<OpenCLThreadManager>() );
  #else
    printf("only enabling serial computations\n");
    set_thread_manager( std::make_shared<SerialThreadManager>() );
  #endif
}

index_t maximum_concurrent_threads()
  { return thread_manager_->maxConcurrentThreads(); }

bool
is_running_threads()
  { return false; }

void
run_threads( std::string& src , std::string& name , index_t n )
{
  thread_manager_->runConcurrentThreads(src,name,n);
}

template<typename type>
void
add_kernel_arg( type* arg , index_t n )
{
  thread_manager_->addArgument( arg , sizeof(type) , n );
}

template<typename type>
void
retrieve_kernel_value( index_t k , type* value , index_t n )
{
  thread_manager_->retrieveValue(k,n,sizeof(type),value);
}

void
synchronize()
{
  thread_manager_->synchronize();
}

void
terminate()
{}

template void add_kernel_arg( float* , index_t );
template void add_kernel_arg( real* , index_t );
template void add_kernel_arg( int* , index_t );
template void retrieve_kernel_value(index_t,real*,index_t);
template void retrieve_kernel_value(index_t,float*,index_t);
template void retrieve_kernel_value(index_t,int*,index_t);

} // ProcessGPU

} // ursa
