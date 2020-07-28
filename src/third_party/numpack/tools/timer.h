// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef SANS_TIMER_H
#define SANS_TIMER_H


#if defined(_MPI)
#include <boost/mpi.hpp>
#elif defined(_OPENMP)
#include <omp.h>
#else
//#include <boost/date_time/posix_time/posix_time_types.hpp>
#include <chrono>
#endif

namespace numpack
{
#if defined(_MPI)

  class timer
  {
    boost::mpi::timer time;
  public:
    double elapsed() const
    {
      return time.elapsed();
    }
  };

#elif defined(_OPENMP)

  class timer
  {
    double starttime;
  public:
    timer()
    {
      starttime = omp_get_wtime();
    }

    double elapsed() const
    {
      return omp_get_wtime() - starttime;
    }
  };

#else

#if 0
  class timer
  {
    boost::posix_time::ptime starttime;
  public:
    timer()
    {
      starttime = boost::posix_time::microsec_clock::local_time();
    }

    double elapsed() const
    {
      boost::posix_time::ptime endtime = boost::posix_time::microsec_clock::local_time();
      boost::posix_time::time_duration diff = endtime - starttime;
      return (double)diff.total_milliseconds()/((double)1000.);
    }
  };
#else
  class timer
  {
    typedef std::chrono::high_resolution_clock clock_type;
    typedef clock_type::time_point time_type;
    time_type starttime;
  public:
    timer()
    {
      starttime = clock_type::now();
    }

    double elapsed() const
    {
      time_type endtime = clock_type::now();
      return std::chrono::duration<double>(endtime-starttime).count();
    }
  };
#endif


#endif
}


#endif //SANS_TIMER_H
