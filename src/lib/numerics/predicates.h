//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
/*
 *  Copyright (c) 2012-2014, Bruno Levy
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *  * Neither the name of the ALICE Project-Team nor the names of its
 *  contributors may be used to endorse or promote products derived from this
 *  software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact: Bruno Levy
 *
 *     Bruno.Levy@inria.fr
 *     http://www.loria.fr/~levy
 *
 *     ALICE Project
 *     LORIA, INRIA Lorraine,
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX
 *     FRANCE
 *
 */


/*
 *  This file is a PSM (pluggable software module)
 *   generated from the distribution of Geogram.
 *
 *  See Geogram documentation on:
 *   http://alice.loria.fr/software/geogram/doc/html/index.html
 *
 *  See documentation of the functions bundled in this PSM on:
 *   http://alice.loria.fr/software/geogram/doc/html/namespaceGEO_1_1PCK.html
 */



/******* extracted from ../api/defs.h *******/

#ifndef GEOGRAM_API_DEFS
#define GEOGRAM_API_DEFS

#include "avro_types.h"
using namespace avro;

/*
 * Deactivate warnings about documentation
 * We do that, because CLANG's doxygen parser does not know
 * some doxygen commands that we use (retval, copydoc) and
 * generates many warnings for them...
 */
#if defined(__INTEL_COMPILER)
// skip this
#elif defined(__clang__)
#pragma clang diagnostic ignored "-Wunknown-pragmas"
#pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
#endif


#if defined(GEO_DYNAMIC_LIBS)
   #if defined(_MSC_VER)
      #define GEO_IMPORT __declspec(dllimport)
      #define GEO_EXPORT __declspec(dllexport)
   #elif defined(__GNUC__)
      #define GEO_IMPORT
      #define GEO_EXPORT __attribute__ ((visibility("default")))
   #else
      #define GEO_IMPORT
      #define GEO_EXPORT
   #endif
#else
   #define GEO_IMPORT
   #define GEO_EXPORT
#endif

#ifdef geogram_EXPORTS
#define GEOGRAM_API GEO_EXPORT
#else
#define GEOGRAM_API GEO_IMPORT
#endif


#define NO_GEOGRAM_API

typedef int GeoMesh;

/*
typedef unsigned char geo_coord_index_t;

typedef unsigned int geo_index_t;

typedef int geo_signed_index_t;

typedef double geo_coord_t;

typedef int geo_boolean;
*/

// pcaplan: defined in common/types.h
typedef index_t geo_index_t;
typedef coord_t geo_coord_t;
typedef coord_t geo_coord_index_t;

// pcaplan: not defined in avro
typedef int geo_signed_index_t;
typedef int geo_boolean;

enum {
    GEO_FALSE = 0,
    GEO_TRUE = 1
};

#endif


/******* extracted from ../basic/common.h *******/

#ifndef GEOGRAM_BASIC_COMMON
#define GEOGRAM_BASIC_COMMON


// iostream should be included before anything else,
// otherwise 'cin', 'cout' and 'cerr' will be uninitialized.
#include <iostream>



namespace GEO {

    void GEOGRAM_API initialize();

    void GEOGRAM_API terminate();
}


#ifdef NDEBUG
#undef GEO_DEBUG
#undef GEO_PARANOID
#else
#define GEO_DEBUG
#define GEO_PARANOID
#endif

// =============================== LINUX defines ===========================

#if defined(__ANDROID__)
#define GEO_OS_ANDROID
#endif

#if defined(__linux__)

#define GEO_OS_LINUX
#define GEO_OS_UNIX

#ifndef GEO_OS_ANDROID
#define GEO_OS_X11
#endif

#if defined(_OPENMP)
#  define GEO_OPENMP
#endif

#if defined(__INTEL_COMPILER)
#  define GEO_COMPILER_INTEL
#elif defined(__clang__)
#  define GEO_COMPILER_CLANG
#elif defined(__GNUC__)
#  define GEO_COMPILER_GCC
#else
#  error "Unsupported compiler"
#endif

// The following works on GCC and ICC
#if defined(__x86_64)
#  define GEO_ARCH_64
#else
#  define GEO_ARCH_32
#endif

// =============================== WINDOWS defines =========================

#elif defined(WIN32) || defined(_WIN64)

#define GEO_OS_WINDOWS

#if defined(_OPENMP)
#  define GEO_OPENMP
#endif

#if defined(_MSC_VER)
#  define GEO_COMPILER_MSVC
#else
#  error "Unsupported compiler"
#endif

#if defined(_WIN64)
#  define GEO_ARCH_64
#else
#  define GEO_ARCH_32
#endif

// =============================== APPLE defines ===========================

#elif defined(__APPLE__)

#define GEO_OS_APPLE
#define GEO_OS_UNIX

#if defined(_OPENMP)
#  define GEO_OPENMP
#endif

#if defined(__clang__)
#  define GEO_COMPILER_CLANG
#elif defined(__GNUCC__)
#define GEO_COMPILER_GCC
#else
//#  error "Unsupported compiler"
#endif

#if defined(__x86_64) || defined(__ppc64__)
#  define GEO_ARCH_64
#else
#  define GEO_ARCH_32
#endif

// =============================== Emscripten defines  ======================

#elif defined(__EMSCRIPTEN__)

#define GEO_OS_UNIX
#define GEO_OS_LINUX
#define GEO_OS_EMSCRIPTEN
#define GEO_ARCH_64
#define GEO_COMPILER_EMSCRIPTEN

// =============================== Unsupported =============================
#else

#error "Unsupported operating system"

#endif

#ifdef DOXYGEN_ONLY
// Keep doxygen happy
#define GEO_OS_WINDOWS
#define GEO_OS_APPLE
#define GEO_OS_ANDROID
#define GEO_ARCH_32
#define GEO_COMPILER_INTEL
#define GEO_COMPILER_MSVC
#endif

#define CPP_CONCAT_(A, B) A ## B

#define CPP_CONCAT(A, B) CPP_CONCAT_(A, B)

#if defined(GOMGEN)
#define GEO_NORETURN
#elif defined(GEO_COMPILER_CLANG) || \
    defined(GEO_COMPILER_GCC)   || \
    defined(GEO_COMPILER_EMSCRIPTEN)
#define GEO_NORETURN __attribute__((noreturn))
#else
#define GEO_NORETURN
#endif

#if defined(GOMGEN)
#define GEO_NORETURN_DECL
#elif defined(GEO_COMPILER_MSVC)
#define GEO_NORETURN_DECL __declspec(noreturn)
#else
#define GEO_NORETURN_DECL
#endif

#if defined(GEO_COMPILER_CLANG)
#if __has_feature(cxx_noexcept)
#define GEO_NOEXCEPT noexcept
#endif
#endif

#ifndef GEO_NOEXCEPT
#define GEO_NOEXCEPT throw()
#endif

#endif


/******* extracted from ../basic/argused.h *******/

#ifndef GEOGRAM_BASIC_ARGUSED
#define GEOGRAM_BASIC_ARGUSED



namespace GEO {

    template <class T>
    inline void geo_argused(const T&) {
    }
}

#endif


/******* extracted from ../basic/numeric.h *******/

#ifndef GEOGRAM_BASIC_NUMERIC
#define GEOGRAM_BASIC_NUMERIC

#include <math.h>
#include <float.h>
#include <limits.h>

// Visual C++ ver. < 2010 does not have C99 stdint.h,
// using a fallback portable one.
#if defined(GEO_OS_WINDOWS) && (_MSC_VER < 1600)
#else
#include <stdint.h>
#endif

#include <limits>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


namespace GEO {

    namespace Numeric {


        typedef void* pointer;


        typedef int8_t int8;


        typedef int16_t int16;


        typedef int32_t int32;


        typedef int64_t int64;


        typedef uint8_t uint8;


        typedef uint16_t uint16;


        typedef uint32_t uint32;


        typedef uint64_t uint64;


        typedef float float32;


        typedef double float64;

        inline float32 max_float32() {
            return std::numeric_limits<float32>::max();
        }

        inline float32 min_float32() {
            // Note: numeric_limits<>::min() is not
            // what we want (it returns the smallest
            // positive non-denormal).
            return -max_float32();
        }

        inline float64 max_float64() {
            return std::numeric_limits<float64>::max();
        }

        inline float64 min_float64() {
            // Note: numeric_limits<>::min() is not
            // what we want (it returns the smallest
            // positive non-denormal).
            return -max_float64();
        }

        bool GEOGRAM_API is_nan(float32 x);

        bool GEOGRAM_API is_nan(float64 x);

        void GEOGRAM_API random_reset();

        int32 GEOGRAM_API random_int32();

        float32 GEOGRAM_API random_float32();

        float64 GEOGRAM_API random_float64();

        template <class T, bool is_numeric>
        struct LimitsHelper : std::numeric_limits<T> {
        };

        template <class T>
        struct LimitsHelper<T, true> : std::numeric_limits<T> {

            static const size_t size = sizeof(T);

            static const size_t numbits = 8 * sizeof(T);
        };

        template <class T>
        struct Limits :
            LimitsHelper<T, std::numeric_limits<T>::is_specialized> {
        };
    }



    template <class T>
    inline T geo_max(T x1, T x2) {
        return x1 < x2 ? x2 : x1;
    }

    template <class T>
    inline T geo_min(T x1, T x2) {
        return x2 < x1 ? x2 : x1;
    }

    enum Sign {

        NEGATIVE = -1,

        ZERO = 0,

        POSITIVE = 1
    };

    template <class T>
    inline Sign geo_sgn(const T& x) {
        return (x > 0) ? POSITIVE : (
            (x < 0) ? NEGATIVE : ZERO
        );
    }

    template <class T>
    inline T geo_abs(T x) {
        return (x < 0) ? -x : x;
    }

    template <class T>
    inline T geo_sqr(T x) {
        return x * x;
    }

    template <class T>
    inline void geo_clamp(T& x, T min, T max) {
        if(x < min) {
            x = min;
        } else if(x > max) {
            x = max;
        }
    }

    template <class T>
    inline void geo_swap(T& x, T& y) {
        T z = x;
        x = y;
        y = z;
    }

    //typedef geo_index_t index_t;

    inline index_t max_index_t() {
        return std::numeric_limits<index_t>::max();
    }

    typedef geo_signed_index_t signed_index_t;

    inline signed_index_t max_signed_index_t() {
        return std::numeric_limits<signed_index_t>::max();
    }

    inline signed_index_t min_signed_index_t() {
        return std::numeric_limits<signed_index_t>::min();
    }

    typedef geo_coord_index_t coord_index_t;
}

#endif


/******* extracted from ../basic/psm.h *******/

#ifndef GEOGRAM_BASIC_PSM
#define GEOGRAM_BASIC_PSM


#include <assert.h>
#include <iostream>
#include <string>

#include "common/error.h"

#define GEOGRAM_PSM

//#define geo_assert(x) if(!(x)) {}
#define geo_assert(x) avro_assert(x)
//#define geo_assert(x) if (!(x)) throw("error");
#define geo_range_assert(x, min_val, max_val) \
    assert((x) >= (min_val) && (x) <= (max_val))
#define geo_assert_not_reached assert(0)

#ifdef GEO_DEBUG
#define geo_debug_assert(x) assert(x)
#define geo_debug_range_assert(x, min_val, max_val) \
    assert((x) >= (min_val) && (x) <= (max_val))
#else
#define geo_debug_assert(x)
#define geo_debug_range_assert(x, min_val, max_val)
#endif

#ifdef GEO_PARANOID
#define geo_parano_assert(x) geo_assert(x)
#define geo_parano_range_assert(x, min_val, max_val) \
    geo_range_assert(x, min_val, max_val)
#else
#define geo_parano_assert(x)
#define geo_parano_range_assert(x, min_val, max_val)
#endif


namespace GEO {
    namespace Process {

        typedef int spinlock;

        inline void acquire_spinlock(spinlock& x) {
            // Not implemented yet for PSMs
            geo_argused(x);
            //geo_assert_not_reached; // pcaplan
        }

        inline void release_spinlock(spinlock& x) {
            // Not implemented yet for PSMs
            geo_argused(x);
            //geo_assert_not_reached; // pcaplan
        }
    }


    namespace Logger {
        inline std::ostream& out(const std::string& name) {
            return std::cout << " [" << name << "]";
        }

        inline std::ostream& err(const std::string& name) {
            return std::cerr << "E[" << name << "]";
        }

        inline std::ostream& warn(const std::string& name) {
            return std::cerr << "W[" << name << "]";
        }
    }

}

#ifndef FPG_UNCERTAIN_VALUE
#define FPG_UNCERTAIN_VALUE 0
#endif

#endif

namespace GEO {

  class GEOGRAM_API expansion;

  #define expansion_sum(a, b)            \
      new_expansion_on_stack(           \
          expansion::sum_capacity(a, b)   \
      )->assign_sum(a, b)

  #define expansion_sum3(a, b, c)          \
      new_expansion_on_stack(            \
          expansion::sum_capacity(a, b, c) \
      )->assign_sum(a, b, c)


  #define expansion_sum4(a, b, c, d)          \
      new_expansion_on_stack(              \
          expansion::sum_capacity(a, b, c, d) \
      )->assign_sum(a, b, c, d)

  #define expansion_diff(a, b)             \
      new_expansion_on_stack(             \
          expansion::diff_capacity(a, b)   \
      )->assign_diff(a, b)

  #define expansion_product(a, b)            \
      new_expansion_on_stack(               \
          expansion::product_capacity(a, b)  \
      )->assign_product(a, b)

  #define expansion_product3(a, b, c)           \
      new_expansion_on_stack(                 \
          expansion::product_capacity(a, b, c)  \
      )->assign_product(a, b, c)

  #define expansion_square(a)             \
      new_expansion_on_stack(             \
          expansion::square_capacity(a)   \
      )->assign_square(a)

      // =============== determinants =====================================

  #define expansion_det2x2(a11, a12, a21, a22)          \
      new_expansion_on_stack(                        \
          expansion::det2x2_capacity(a11, a12, a21, a22) \
      )->assign_det2x2(a11, a12, a21, a22)

  #define expansion_det3x3(a11, a12, a13, a21, a22, a23, a31, a32, a33)   \
      new_expansion_on_stack(                                             \
          expansion::det3x3_capacity(a11,a12,a13,a21,a22,a23,a31,a32,a33) \
      )->assign_det3x3(a11, a12, a13, a21, a22, a23, a31, a32, a33)

  #define expansion_det_111_2x3(a21, a22, a23, a31, a32, a33)           \
      new_expansion_on_stack(                                      \
          expansion::det_111_2x3_capacity(a21, a22, a23, a31, a32, a33) \
      )->assign_det_111_2x3(a21, a22, a23, a31, a32, a33)

      // =============== geometric functions ==============================

  #define expansion_det4x4(a11,a12,a13,a14,a21,a22,a23,a24,a31,a32,a33,a34,a41,a42,a43,a44)  \
    new_expansion_on_stack(                                             \
        expansion::det4x4_capacity(a11,a12,a13,a14,a21,a22,a23,a24,a31,a32,a33,a34,a41,a42,a43,a44) \
    )->assign_det4x4(a11, a12, a13, a14, a21, a22, a23, a24, a31, a32, a33, a34, a41, a42, a43, a44 )

  #define expansion_det5x5(a11,a12,a13,a14,a15,a21,a22,a23,a24,a25,a31,a32,a33,a34,a35,a41,a42,a43,a44,a45,a51,a52,a53,a54,a55) \
  expansion_sum3( \
    expansion_diff( expansion_product(a11,expansion_det4x4(a22,a23,a24,a25,a32,a33,a34,a35,a42,a43,a44,a45,a52,a53,a54,a55) ) , \
                    expansion_product(a12,expansion_det4x4(a21,a23,a24,a25,a31,a33,a34,a35,a41,a43,a44,a45,a51,a53,a54,a55) ) ) , \
    expansion_diff( expansion_product(a13,expansion_det4x4(a21,a22,a24,a25,a31,a32,a34,a35,a41,a42,a44,a45,a51,a52,a54,a55) ) , \
                    expansion_product(a14,expansion_det4x4(a21,a22,a23,a25,a31,a32,a33,a35,a41,a42,a43,a45,a51,a52,a53,a55) ) ) , \
                    expansion_product(a15,expansion_det4x4(a21,a22,a23,a24,a31,a32,a33,a34,a41,a42,a43,a44,a51,a52,a53,a54) ) \
  )

  #define expansion_det_1111_3x4(a21,a22,a23,a24,a31,a32,a33,a34,a41,a42,a43,a44) \
    new_expansion_on_stack(                                      \
        expansion::det_1111_3x4_capacity(a21, a22, a23, a24, a31, a32, a33, a34, a41, a42, a43, a44) \
    )->assign_det_1111_3x4(a21, a22, a23, a24, a31, a32, a33, a34, a41, a42, a43, a44) 

  #define expansion_sq_dist(a, b, dim)           \
      new_expansion_on_stack(                  \
          expansion::sq_dist_capacity(dim)     \
      )->assign_sq_dist(a, b, dim)

  #define expansion_dot_at(a, b, c, dim)           \
      new_expansion_on_stack(                   \
          expansion::dot_at_capacity(dim)       \
      )->assign_dot_at(a, b, c, dim)


  #define expansion_length2(x,y,z)              \
      new_expansion_on_stack(                   \
         expansion::length2_capacity(x,y,z)     \
      )->assign_length2(x,y,z)



      Sign GEOGRAM_API sign_of_expansion_determinant(
          const expansion& a00,const expansion& a01,
          const expansion& a10,const expansion& a11
      );

      Sign GEOGRAM_API sign_of_expansion_determinant(
          const expansion& a00,const expansion& a01,const expansion& a02,
          const expansion& a10,const expansion& a11,const expansion& a12,
          const expansion& a20,const expansion& a21,const expansion& a22
      );

      Sign GEOGRAM_API sign_of_expansion_determinant(
          const expansion& a00,const expansion& a01,
          const expansion& a02,const expansion& a03,
          const expansion& a10,const expansion& a11,
          const expansion& a12,const expansion& a13,
          const expansion& a20,const expansion& a21,
          const expansion& a22,const expansion& a23,
          const expansion& a30,const expansion& a31,
          const expansion& a32,const expansion& a33
      );
}

/******* extracted from predicates.h *******/

#ifndef GEOGRAM_NUMERICS_PREDICATES
#define GEOGRAM_NUMERICS_PREDICATES



namespace GEO {

    namespace PCK {

       // predicate signatures
       Sign GEOGRAM_API side1_SOS(const double* p0,const double* p1,const double* q0, coord_index_t DIM);
       Sign GEOGRAM_API side2_SOS(const double* p0,const double* p1,const double* p2,const double* q0,const double* q1, coord_index_t DIM);
       Sign GEOGRAM_API side3_SOS(const double* p0,const double* p1,const double* p2,const double* p3,const double* q0,const double* q1,const double* q2, coord_index_t DIM);
       Sign GEOGRAM_API side4_SOS(const double* p0,const double* p1,const double* p2,const double* p3,const double* p4,const double* q0,const double* q1,const double* q2,const double* q3, coord_index_t DIM);
       Sign GEOGRAM_API side5_SOS(const double* p0,const double* p1,const double* p2,const double* p3,const double* p4,const double* p5,const double* q0,const double* q1,const double* q2,const double* q3,const double* q4, coord_index_t DIM);

       void GEOGRAM_API initialize();

    }
}

#endif
