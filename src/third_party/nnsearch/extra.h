#ifndef GEOGRAM_BASIC_EXTRA
#define GEOGRAM_BASIC_EXTRA

#include "common/process.h"

#include "common.h"
#include "defs.h"
#include "numeric.h"

#include <cmath>
#include <float.h>
#include <limits.h>
#include <algorithm> // for std::min / std::max

namespace GEO
{

/**
 * \def geo_restrict
 * \brief Informs the compiler that a given pointer has no aliasing
 * \details
 *  No aliasing means that no other pointer points to the same area of
 *  memory.
 * \code
 * double* geo_restrict p = ...;
 * \endcode
 */
#if   defined(GEO_COMPILER_INTEL)
#define geo_restrict __restrict
#elif defined(GEO_COMPILER_GCC_FAMILY)
#define geo_restrict __restrict__
#elif defined(GEO_COMPILER_MSVC)
#define geo_restrict __restrict
#elif defined(GEO_COMPILER_EMSCRIPTEN)
#define geo_restrict
#endif

namespace Geom
{

inline double
distance2( const double* x , const double* y, unsigned short dim )
{
  double result = 0.0;
  for (unsigned short d=0;d<dim;d++)
    result += (x[d] -y[d])*(x[d] -y[d]);
  return result;
}

} // Geom

/**
 * \brief Used by the implementation of GEO::parallel()
 * \see GEO::parallel()
 */
class ParallelThread : public avro::Thread {
public:
/**
* \brief ParallelThread constructor.
* \param[in] func a void function with no parameter.
*/
ParallelThread(
  std::function<void(void)> func
) : func_(func) {
}

/**
* \copydoc Thread::run()
*/
    void run() override {
  func_();
    }
private:
std::function<void()> func_;
};


/**
 * \brief Used by the implementation of GEO::parallel_for()
 * \see GEO::parallel_for()
 */
class ParallelForThread : public avro::Thread {
public:

/**
* \param[in] func a void function that takes an index_t
* \param[in] from the first iteration index
* \param[in] to one position past the last interation index
* \param[in] step iteration step
*/
ParallelForThread(
  std::function<void(index_t)> func,
  index_t from, index_t to, index_t step=1
) : func_(func), from_(from), to_(to), step_(step) {
}

/**
* \copydoc Thread::run()
*/
    void run() override {
        for(index_t i = from_; i < to_; i += step_) {
            func_(i);
        }
    }
private:
std::function<void(index_t)> func_;
index_t from_;
index_t to_;
index_t step_;
};

/**
 * \brief Used by the implementation of GEO::parallel_for_slice()
 * \see GEO::parallel_for_slice()
 */
class ParallelForSliceThread : public avro::Thread {
public:

/**
* \param[in] func a void function that takes two index_t arguments
* \param[in] from the first iteration index
* \param[in] to one position past the last interation index
*/
ParallelForSliceThread(
  std::function<void(index_t,index_t)> func,
  index_t from, index_t to
) : func_(func), from_(from), to_(to) {
}

/**
* \copydoc Thread::run()
*/
    void run() override {
  func_(from_, to_);
    }
private:
std::function<void(index_t,index_t)> func_;
index_t from_;
index_t to_;
};

// parallel

inline void parallel_for(
    index_t from, index_t to, std::function<void(index_t)> func,
    index_t threads_per_core, bool interleaved
) {
#ifdef GEO_OS_WINDOWS
    // TODO: This is a limitation of WindowsThreadManager, to be fixed.
    threads_per_core = 1;
#endif

    index_t nb_threads = std::min(
        avro::index_t(to - from),
        avro::ProcessCPU::maximum_concurrent_threads() * threads_per_core
    );

nb_threads = std::max(index_t(1), nb_threads);

    index_t batch_size = (to - from) / nb_threads;
    if(avro::ProcessCPU::is_running_threads() || nb_threads == 1) {
        for(index_t i = from; i < to; i++) {
            func(i);
        }
    } else {
        avro::ThreadGroup threads;
        if(interleaved) {
            for(index_t i = 0; i < nb_threads; i++) {
                threads.push_back(
                    std::make_shared<ParallelForThread>(
                        func, from + i, to, nb_threads
                    )
                );
            }
        } else {
            index_t cur = from;
            for(index_t i = 0; i < nb_threads; i++) {
                if(i == nb_threads - 1) {
                    threads.push_back(
                        std::make_shared<ParallelForThread>(
                            func, cur, to
                        )
                    );
                } else {
                    threads.push_back(
                        std::make_shared<ParallelForThread>(
                            func, cur, cur + batch_size
                        )
                    );
                }
                cur += batch_size;
            }
        }
        avro::ProcessCPU::run_threads(threads);
    }
}


inline void parallel_for_slice(
index_t from, index_t to, std::function<void(index_t, index_t)> func,
    index_t threads_per_core
) {
#ifdef GEO_OS_WINDOWS
    // TODO: This is a limitation of WindowsThreadManager, to be fixed.
    threads_per_core = 1;
#endif

    index_t nb_threads = std::min(
        avro::index_t(to - from),
        avro::ProcessCPU::maximum_concurrent_threads() * threads_per_core
    );

nb_threads = std::max(index_t(1), nb_threads);

    index_t batch_size = (to - from) / nb_threads;
    if(avro::ProcessCPU::is_running_threads() || nb_threads == 1) {
  func(from, to);
    } else {
        avro::ThreadGroup threads;
  index_t cur = from;
  for(index_t i = 0; i < nb_threads; i++) {
if(i == nb_threads - 1) {
    threads.push_back(
      std::make_shared<ParallelForSliceThread>(
      func, cur, to
    )
  );
} else {
    threads.push_back(
  std::make_shared<ParallelForSliceThread>(
      func, cur, cur + batch_size
                       )
                    );
}
cur += batch_size;
  }
        avro::ProcessCPU::run_threads(threads);
    }
}

inline void parallel(
std::function<void()> f1,
std::function<void()> f2
) {
    if(avro::ProcessCPU::is_running_threads()) {
  f1();
  f2();
    } else {
        avro::ThreadGroup threads;
  threads.push_back(std::make_shared<ParallelThread>(f1));
  threads.push_back(std::make_shared<ParallelThread>(f2));
        avro::ProcessCPU::run_threads(threads);
    }
}


inline void parallel(
std::function<void()> f1,
std::function<void()> f2,
std::function<void()> f3,
std::function<void()> f4
) {
    if(avro::ProcessCPU::is_running_threads()) {
  f1();
  f2();
  f3();
  f4();
    } else {
        avro::ThreadGroup threads;
  threads.push_back(std::make_shared<ParallelThread>(f1));
  threads.push_back(std::make_shared<ParallelThread>(f2));
  threads.push_back(std::make_shared<ParallelThread>(f3));
  threads.push_back(std::make_shared<ParallelThread>(f4));
        avro::ProcessCPU::run_threads(threads);
    }
}


inline void parallel(
std::function<void()> f1,
std::function<void()> f2,
std::function<void()> f3,
std::function<void()> f4,
std::function<void()> f5,
std::function<void()> f6,
std::function<void()> f7,
std::function<void()> f8
) {
    if(avro::ProcessCPU::is_running_threads()) {
  f1();
  f2();
  f3();
  f4();
  f5();
  f6();
  f7();
  f8();
    } else {
        avro::ThreadGroup threads;
  threads.push_back(std::make_shared<ParallelThread>(f1));
  threads.push_back(std::make_shared<ParallelThread>(f2));
  threads.push_back(std::make_shared<ParallelThread>(f3));
  threads.push_back(std::make_shared<ParallelThread>(f4));
  threads.push_back(std::make_shared<ParallelThread>(f5));
  threads.push_back(std::make_shared<ParallelThread>(f6));
  threads.push_back(std::make_shared<ParallelThread>(f7));
  threads.push_back(std::make_shared<ParallelThread>(f8));
        avro::ProcessCPU::run_threads(threads);
    }
}

} // GEO

#endif
