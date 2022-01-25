//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
// include this file within the enclosed namespaces
// needed so that ::maximum_concurrent_threads, ::is_running_threads, ::run_threads
// are defined in the appropriate namespaces

// parallel for loop
template<class Func>
inline void parallel_for( const Func& func , index_t from, index_t to ,
                          index_t threads_per_core = 1, bool interleaved = false )
{
  index_t nb_threads = std::min( to-from , maximum_concurrent_threads()*threads_per_core );
  nb_threads = std::max(1ul,nb_threads);

  index_t batch_size = (to -from)/nb_threads;

  if (is_running_threads() || nb_threads==1)
  {
    // run in serial
    for (index_t i=from;i<to;i++)
      const_cast<Func&>(func)(i);
  }
  else
  {
    ThreadGroup threads;
    if (interleaved)
    {
      for (index_t i=0;i<nb_threads;i++)
        threads.push_back( std::make_shared<ParallelForThread<Func>>(func,from+i,to,nb_threads) );
    }
    else
    {
      index_t cur = from;
      for (index_t i=0;i<nb_threads;i++)
      {
        if (i==nb_threads-1)
          threads.push_back( std::make_shared<ParallelForThread<Func>>(func,cur,to) );
        else
          threads.push_back( std::make_shared<ParallelForThread<Func>>(func,cur,cur+batch_size) );
        cur += batch_size;
      }
    }
    run_threads(threads);
  }
}
