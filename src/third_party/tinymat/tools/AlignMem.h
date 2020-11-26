// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef ALIGN_MEM_H_INCLUDED
#define ALIGN_MEM_H_INCLUDED

#include "tools/SANSException.h"
#include "CacheLineSize.h"

#if (defined __INTEL_COMPILER)
#define ALIGN_MEM __declspec(align(16))

#elif (defined __GNUC__)
#define ALIGN_MEM __attribute__((aligned))

#elif (defined _MSC_VER)
#define ALIGN_MEM

#else
#define ALIGN_MEM
#warning "Please define the equivalent of the gnu __attribute__((aligned)) for this compiler in (src/tools/AlignMem.h)"
#endif


//Memory alligned allocation

#if (defined __INTEL_COMPILER)
#define ALIGNED_MALLOC(T, value, size) value = (T*)_mm_malloc(size*sizeof(T), 16); SANS_ASSERT(value)

#elif (defined __GNUC__)
#define ALIGNED_MALLOC(T, value, size) SANS_ASSERT(!posix_memalign((void**)&value, 16, size*sizeof(T) ))

#elif (defined _MSC_VER)
#define ALIGNED_MALLOC(T, value, size) value = (T*)_aligned_malloc(size*sizeof(T), 16); SANS_ASSERT(value)

#else
#define ALIGNED_MALLOC(T, value, size) value = new T[size]
#warning "Please define the equivalent of posix_memalign for this compiler in (tools/AlignMem.h)"
#endif

//Memory alligned allocation in constructor

//#if (defined __INTEL_COMPILER)
//#define ALIGNED_INIT_MALLOC(T, value, size) value( (T*)_mm_malloc(size*sizeof(T), CACHE_LINE_SIZE) )
//
//#elif (defined __GNUC__)
//#define ALIGNED_INIT_MALLOC(T, value, size) SANS_ASSERT(!posix_memalign((void**)&value, CACHE_LINE_SIZE, size*sizeof(T) ))
//
//#elif (defined _MSC_VER)
//#define ALIGNED_INIT_MALLOC(T, value, size) value( (T*)_aligned_malloc(size*sizeof(T), CACHE_LINE_SIZE) )
//
//#else
//#define ALIGNED_INIT_MALLOC(T, value, size) value(new T[size])
//#warning "Please define the equivalent of posix_memalign for this compiler in (src/tools/AlignMem.h)"
//#endif


//Memory alligned deallocation

#if (defined __INTEL_COMPILER)
#define ALIGNED_FREE(value) _mm_free(value)

#elif (defined __GNUC__)
#define ALIGNED_FREE(value) free(value)

#elif (defined _MSC_VER)
#define ALIGNED_FREE(value) _aligned_free(value)

#else
#define ALIGNED_FREE(value) delete [] value
#warning "Please define the equivalent of the gnu posix_memalign for this compiler in (src/tools/AlignMem.h)"
#endif


#endif //ALIGN_MEM_H_INCLUDED
