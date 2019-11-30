#include "numerics/predicates.h"
#include "numerics/expansion.h"
#include <algorithm>

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



/******* extracted from ../basic/memory.h *******/

#ifndef GEOGRAM_BASIC_MEMORY
#define GEOGRAM_BASIC_MEMORY

#include <vector>
#include <string.h>
#include <stdlib.h>

#ifdef GEO_OS_WINDOWS

#include <windows.h>
#ifdef min
#undef min
#endif
#ifdef max
#undef max
#endif

#else

#include <unistd.h>

#endif


namespace GEO {

    namespace Memory {

        typedef unsigned char byte;


        typedef unsigned char word8;


        typedef unsigned short word16;


        typedef unsigned int word32;


        typedef byte* pointer;

//#define nil 0 // pcaplan: defined in common/types.h

        inline void clear(void* addr, size_t size) {
            ::memset(addr, 0, size);
        }

        inline void copy(void* to, const void* from, size_t size) {
            ::memcpy(to, from, size);
        }

#define GEO_MEMORY_ALIGNMENT 64

        template <int DIM>
        struct PointAlignment {
            static const size_t value = 1;
        };

        template <>
        struct PointAlignment<2> {
            static const size_t value = 16;
        };

        template <>
        struct PointAlignment<3> {
            static const size_t value = 8;
        };

        template <>
        struct PointAlignment<4> {
            static const size_t value = 32;
        };

        template <>
        struct PointAlignment<6> {
            static const size_t value = 16;
        };

        template <>
        struct PointAlignment<8> {
            static const size_t value = 64;
        };

#define geo_dim_alignment(dim) GEO::Memory::PointAlignment<dim>::value

        inline void* aligned_malloc(
            size_t size, size_t alignment = GEO_MEMORY_ALIGNMENT
        ) {
#if   defined(GEO_OS_ANDROID)
            // Alignment not supported under Android.
            geo_argused(alignment);
            return malloc(size);
#elif defined(GEO_COMPILER_INTEL)
            return _mm_malloc(size, alignment);
#elif defined(GEO_COMPILER_GCC) || defined(GEO_COMPILER_CLANG)
            void* result;
            return posix_memalign(&result, alignment, size) == 0
                   ? result : 0;
#elif defined(GEO_COMPILER_MSVC)
            return _aligned_malloc(size, alignment);
#else
            geo_argused(alignment);
            return malloc(size);
#endif
        }

        inline void aligned_free(void* p) {
#if   defined(GEO_OS_ANDROID)
            // Alignment not supported under Android.
            free(p);
#elif defined(GEO_COMPILER_INTEL)
            _mm_free(p);
#elif defined(GEO_COMPILER_GCC) || defined(GEO_COMPILER_CLANG)
            free(p);
#elif defined(GEO_COMPILER_MSVC)
            _aligned_free(p);
#else
            free(p);
#endif
        }

#if   defined(GEO_OS_ANDROID)
#define geo_decl_aligned(var) var
#elif defined(GEO_COMPILER_INTEL)
#define geo_decl_aligned(var) __declspec(aligned(GEO_MEMORY_ALIGNMENT)) var
#elif defined(GEO_COMPILER_GCC) || defined(GEO_COMPILER_CLANG)
#define geo_decl_aligned(var) var __attribute__((aligned(GEO_MEMORY_ALIGNMENT)))
#elif defined(GEO_COMPILER_MSVC)
#define geo_decl_aligned(var) __declspec(align(GEO_MEMORY_ALIGNMENT)) var
#elif defined(GEO_COMPILER_EMSCRIPTEN)
#define geo_decl_aligned(var) var
#endif

#if   defined(GEO_OS_ANDROID)
#define geo_assume_aligned(var, alignment)
#elif defined(GEO_COMPILER_INTEL)
#define geo_assume_aligned(var, alignment) \
    __assume_aligned(var, alignment)
#elif defined(GEO_COMPILER_CLANG)
#define geo_assume_aligned(var, alignment)
        // GCC __builtin_assume_aligned is not yet supported by clang-3.3
#elif defined(GEO_COMPILER_GCC)
#if __GNUC__ >= 4 && __GNUC_MINOR__ >= 7
#define geo_assume_aligned(var, alignment) \
        *(void**) (&var) = __builtin_assume_aligned(var, alignment)
        // the GCC way of specifiying that a pointer is aligned returns
        // the aligned pointer (I can't figure out why). It needs to be
        // affected otherwise it is not taken into account (verified by
        // looking at the output of gcc -S)
#else
#define geo_assume_aligned(var, alignment)
#endif
#elif defined(GEO_COMPILER_MSVC)
#define geo_assume_aligned(var, alignment)
        // TODO: I do not know how to do that with MSVC
#elif defined(GEO_COMPILER_EMSCRIPTEN)
#define geo_assume_aligned(var, alignment)
#endif

#if   defined(GEO_COMPILER_INTEL)
#define geo_restrict restrict
#elif defined(GEO_COMPILER_GCC) || defined(GEO_COMPILER_CLANG)
#define geo_restrict __restrict__
#elif defined(GEO_COMPILER_MSVC)
#define geo_restrict __restrict
#elif defined(GEO_COMPILER_EMSCRIPTEN)
#define geo_restrict
#endif

        inline bool is_aligned(
            void* p, size_t alignment = GEO_MEMORY_ALIGNMENT
        ) {
            return (reinterpret_cast<size_t>(p) & (alignment - 1)) == 0;
        }

        inline void* align(void* p) {
            size_t offset = (
                GEO_MEMORY_ALIGNMENT -
                (reinterpret_cast<size_t>(p) & (GEO_MEMORY_ALIGNMENT - 1))
            ) & (GEO_MEMORY_ALIGNMENT - 1);
            return reinterpret_cast<char*>(p) + offset;
        }

#define geo_aligned_alloca(size) \
    GEO::Memory::align(alloca(size + GEO_MEMORY_ALIGNMENT - 1))

        template <class T, int ALIGN = GEO_MEMORY_ALIGNMENT>
        class aligned_allocator {
        public:

            typedef T value_type;


            typedef T* pointer;


            typedef T& reference;


            typedef const T* const_pointer;


            typedef const T& const_reference;


            typedef ::std::size_t size_type;


            typedef ::std::ptrdiff_t difference_type;

            template <class U>
            struct rebind {

                typedef aligned_allocator<U> other;
            };

            pointer address(reference x) {
                return &x;
            }

            const_pointer address(const_reference x) {
                return &x;
            }

            pointer allocate(
                size_type nb_elt, ::std::allocator<void>::const_pointer hint = 0
            ) {
                geo_argused(hint);
                pointer result = static_cast<pointer>(
                    aligned_malloc(sizeof(T) * nb_elt, ALIGN)
                );
                return result;
            }

            void deallocate(pointer p, size_type nb_elt) {
                geo_argused(nb_elt);
                aligned_free(p);
            }

            size_type max_size() const {
                ::std::allocator<char> a;
                return a.max_size() / sizeof(T);
            }

            void construct(pointer p, const_reference val) {
                new (static_cast<void*>(p))value_type(val);
            }

            void destroy(pointer p) {
                p->~value_type();
#ifdef GEO_COMPILER_MSVC
                (void) p; // to avoid a "unreferenced variable" warning
#endif
            }

            template <class T2, int A2> operator aligned_allocator<T2, A2>() {
                return aligned_allocator<T2,A2>();
            }
        };

        template <typename T1, int A1, typename T2, int A2>
        inline bool operator== (
            const aligned_allocator<T1, A1>&, const aligned_allocator<T2, A2>&
        ) {
            return true;
        }

        template <typename T1, int A1, typename T2, int A2>
        inline bool operator!= (
            const aligned_allocator<T1, A1>&, const aligned_allocator<T2, A2>&
        ) {
            return false;
        }
    }



    template <class T>
    class vector : public ::std::vector<T, Memory::aligned_allocator<T> > {
        typedef ::std::vector<T, Memory::aligned_allocator<T> > baseclass;

    public:
        vector() :
            baseclass() {
        }

        explicit vector(index_t size) :
            baseclass(size) {
        }

        explicit vector(index_t size, const T& val) :
            baseclass(size, val) {
        }

        index_t size() const {
            //   casts baseclass::size() from size_t (64 bits)
            //   to index_t (32 bits), because all
            //   indices in Vorpaline are supposed to fit in 32 bits (index_t).
            // TODO: geo_debug_assert(baseclass::size() < max index_t)
            return index_t(baseclass::size());
        }

        T& operator[] (index_t i) {
            geo_debug_assert(i < size());
            return baseclass::operator[] (i);
        }

        const T& operator[] (index_t i) const {
            geo_debug_assert(i < size());
            return baseclass::operator[] (i);
        }

        T& operator[] (signed_index_t i) {
            geo_debug_assert(i >= 0 && index_t(i) < size());
            return baseclass::operator[] (index_t(i));
        }

        const T& operator[] (signed_index_t i) const {
            geo_debug_assert(i >= 0 && index_t(i) < size());
            return baseclass::operator[] (index_t(i));
        }

        T* data() {
            return size() == 0 ? nil : &(*this)[0];
        }

        const T* data() const {
            return size() == 0 ? nil : &(*this)[0];
        }

    };

    template <>
    class vector<bool> : public ::std::vector<bool> {
        typedef ::std::vector<bool> baseclass;

    public:

        vector() :
            baseclass() {
        }


        explicit vector(index_t size) :
            baseclass(size) {
        }


        explicit vector(index_t size, bool val) :
            baseclass(size, val) {
        }


        index_t size() const {
            //   casts baseclass::size() from size_t (64 bits)
            //   to index_t (32 bits), because all
            //   indices in Vorpaline are supposed to fit in 32 bits (index_t).
            // TODO: geo_debug_assert(baseclass::size() < max index_t)
            return index_t(baseclass::size());
        }

        // TODO: operator[] with bounds checking (more complicated
        // than just returning bool&, check implementation in STL).
    };
}

#endif


/******* extracted from ../basic/matrix.h *******/

#ifndef GEOGRAM_BASIC_MATRIX
#define GEOGRAM_BASIC_MATRIX



namespace GEO {



    inline double det2x2(
        double a11, double a12,
        double a21, double a22
    ) {
        return a11*a22-a12*a21 ;
    }

    inline double det3x3(
        double a11, double a12, double a13,
        double a21, double a22, double a23,
        double a31, double a32, double a33
    ) {
    return
         a11*det2x2(a22,a23,a32,a33)
        -a21*det2x2(a12,a13,a32,a33)
        +a31*det2x2(a12,a13,a22,a23);
    }


    inline double det4x4(
        double a11, double a12, double a13, double a14,
        double a21, double a22, double a23, double a24,
        double a31, double a32, double a33, double a34,
        double a41, double a42, double a43, double a44
    ) {
        double m12 = a21*a12 - a11*a22;
        double m13 = a31*a12 - a11*a32;
        double m14 = a41*a12 - a11*a42;
        double m23 = a31*a22 - a21*a32;
        double m24 = a41*a22 - a21*a42;
        double m34 = a41*a32 - a31*a42;

        double m123 = m23*a13 - m13*a23 + m12*a33;
        double m124 = m24*a13 - m14*a23 + m12*a43;
        double m134 = m34*a13 - m14*a33 + m13*a43;
        double m234 = m34*a23 - m24*a33 + m23*a43;

        return (m234*a14 - m134*a24 + m124*a34 - m123*a44);
    }


    template <index_t DIM, class FT>
    class Matrix {
    public:

        typedef Matrix<DIM, FT> matrix_type;


        typedef FT value_type;


        static const index_t dim = DIM;

        inline Matrix() {
            load_identity();
        }

        explicit Matrix(const FT* vals) {
            for(index_t i = 0; i < DIM; i++) {
                for(index_t j = 0; j < DIM; j++) {
                    coeff_[i][j] = *vals;
                    ++vals;
                }
            }
        }

        inline index_t dimension() const {
            return DIM;
        }

        inline void load_zero() {
            for(index_t i = 0; i < DIM; i++) {
                for(index_t j = 0; j < DIM; j++) {
                    coeff_[i][j] = FT(0);
                }
            }
        }

        inline void load_identity() {
            for(index_t i = 0; i < DIM; i++) {
                for(index_t j = 0; j < DIM; j++) {
                    coeff_[i][j] = (i == j) ? FT(1) : FT(0);
                }
            }
        }

        inline bool is_identity() const {
            for(index_t i = 0; i < DIM; i++) {
                for(index_t j = 0; j < DIM; j++) {
                    FT rhs = ((i == j) ? FT(1) : FT(0));
                    if(coeff_[i][j] != rhs) {
                        return false;
                    }
                }
            }
            return true;
        }

        inline FT& operator() (index_t i, index_t j) {
            geo_debug_assert(i < DIM);
            geo_debug_assert(j < DIM);
            return coeff_[i][j];
        }

        inline const FT& operator() (index_t i, index_t j) const {
            geo_debug_assert(i < DIM);
            geo_debug_assert(j < DIM);
            return coeff_[i][j];
        }

        inline matrix_type& operator+= (const matrix_type& m) {
            for(index_t i = 0; i < DIM; i++) {
                for(index_t j = 0; j < DIM; j++) {
                    coeff_[i][j] += m.coeff_[i][j];
                }
            }
            return *this;
        }

        inline matrix_type& operator-= (const matrix_type& m) {
            for(index_t i = 0; i < DIM; i++) {
                for(index_t j = 0; j < DIM; j++) {
                    coeff_[i][j] -= m.coeff_[i][j];
                }
            }
            return *this;
        }

        inline matrix_type& operator*= (FT val) {
            for(index_t i = 0; i < DIM; i++) {
                for(index_t j = 0; j < DIM; j++) {
                    coeff_[i][j] *= val;
                }
            }
            return *this;
        }

        inline matrix_type& operator/= (FT val) {
            for(index_t i = 0; i < DIM; i++) {
                for(index_t j = 0; j < DIM; j++) {
                    coeff_[i][j] /= val;
                }
            }
            return *this;
        }

        inline matrix_type operator+ (const matrix_type& m) const {
            matrix_type result = *this;
            result += m;
            return result;
        }

        inline matrix_type operator- (const matrix_type& m) const {
            matrix_type result = *this;
            result -= m;
            return result;
        }

        inline matrix_type operator* (FT val) const {
            matrix_type result = *this;
            result *= val;
            return result;
        }

        inline matrix_type operator/ (FT val) const {
            matrix_type result = *this;
            result /= val;
            return result;
        }

        matrix_type operator* (const matrix_type& m) const {
            matrix_type result;
            for(index_t i = 0; i < DIM; i++) {
                for(index_t j = 0; j < DIM; j++) {
                    result.coeff_[i][j] = FT(0);
                    for(index_t k = 0; k < DIM; k++) {
                        result.coeff_[i][j] += coeff_[i][k] * m.coeff_[k][j];
                    }
                }
            }
            return result;
        }

        matrix_type inverse() const {
            matrix_type result;
            bool invertible = compute_inverse(result);
            geo_assert(invertible);
            return result;
        }


        bool compute_inverse(matrix_type& result) const {
            FT val=FT(0.0), val2=FT(0.0);
            matrix_type tmp = (*this);

            result.load_identity();

            for(index_t i = 0; i != DIM; i++) {
                val = tmp(i, i);                     /* find pivot */
                index_t ind = i;
                for(index_t j = i + 1; j != DIM; j++) {
                    if(fabs(tmp(j, i)) > fabs(val)) {
                        ind = j;
                        val = tmp(j, i);
                    }
                }

                if(ind != i) {
                    for(index_t j = 0; j != DIM; j++) {
                        val2 = result(i, j);
                        result(i, j) = result(ind, j);
                        result(ind, j) = val2;           /* swap columns */
                        val2 = tmp(i, j);
                        tmp(i, j) = tmp(ind, j);
                        tmp(ind, j) = val2;
                    }
                }

                if(val == 0.0) {
                    return false;
                }

                for(index_t j = 0; j != DIM; j++) {
                    tmp(i, j) /= val;
                    result(i, j) /= val;
                }

                for(index_t j = 0; j != DIM; j++) {
                    if(j == i) {
                        continue;                       /* eliminate column */
                    }
                    val = tmp(j, i);
                    for(index_t k = 0; k != DIM; k++) {
                        tmp(j, k) -= tmp(i, k) * val;
                        result(j, k) -= result(i, k) * val;
                    }
                }
            }

            return true;
        }

        matrix_type transpose() const {
            matrix_type result;
            for(index_t i = 0; i < DIM; i++) {
                for(index_t j = 0; j < DIM; j++) {
                    result(i, j) = (* this)(j, i);
                }
            }
            return result;
        }



        inline const FT* data() const {
            return &(coeff_[0][0]);
        }



        inline FT* data() {
            return &(coeff_[0][0]);
        }

        void get_lower_triangle(FT* store) const {
            for(index_t i = 0; i < DIM; i++) {
                for(index_t j = 0; j <= i; j++) {
                    *store++ = coeff_[i][j];
                }
            }
        }

    private:
        FT coeff_[DIM][DIM];
    };



    template <index_t DIM, class FT>
    inline std::ostream& operator<< (
        std::ostream& output, const Matrix<DIM, FT>& m
    ) {
        const char* sep = "";
        for(index_t i = 0; i < DIM; i++) {
            for(index_t j = 0; j < DIM; j++) {
                output << sep << m(i, j);
                sep = " ";
            }
        }
        return output;
    }

    template <index_t DIM, class FT>
    inline std::istream& operator>> (
        std::istream& input, Matrix<DIM, FT>& m
    ) {
        for(index_t i = 0; i < DIM; i++) {
            for(index_t j = 0; j < DIM; j++) {
                input >> m(i, j);
            }
        }
        return input;
    }



    template <index_t DIM, class FT> inline
    void mult(const Matrix<DIM, FT>& M, const FT* x, FT* y) {
        for(index_t i = 0; i < DIM; i++) {
            y[i] = 0;
            for(index_t j = 0; j < DIM; j++) {
                y[i] += M(i, j) * x[j];
            }
        }
    }

}

#endif


// This makes sure the compiler will not optimize y = a*x+b
// with fused multiply-add, this would break the exact
// predicates.
#ifdef GEO_COMPILER_MSVC
#pragma fp_contract(off)
#endif


namespace {

    using namespace GEO;



    bool expansion_length_stat_ = false;
    std::vector<index_t> expansion_length_histo_;

    class Pools {
    public:

        Pools() : pools_(1024,static_cast<void*>(0)) {
            chunks_.reserve(1024);
        }

        ~Pools() {
            for(index_t i=0; i<chunks_.size(); ++i) {
                delete[] chunks_[i];
            }
        }

        void* malloc(size_t size) {
            if(size >= pools_.size()) {
                return ::malloc(size);
            }
            if(pools_[size] == nil) {
                new_chunk(size);
            }
            void* result = pools_[size];
            pools_[size] = *static_cast<void**>(pools_[size]);
            return result;
        }

        void free(void* ptr, size_t size) {
            if(size >= pools_.size()) {
                ::free(ptr);
                return;
            }
            *static_cast<void**>(ptr) = pools_[size];
            pools_[size] = ptr;
        }


    protected:
        static const index_t POOL_CHUNK_SIZE = 512;

        void new_chunk(size_t size_in) {
            size_t size = (size_in / 8 + 1)*8; // Align memory.
            Memory::pointer chunk = new Memory::byte[size * POOL_CHUNK_SIZE];
            for(index_t i=0; i<POOL_CHUNK_SIZE-1; ++i) {
                Memory::pointer cur = chunk + size * i;
                Memory::pointer next = cur + size;
                *reinterpret_cast<void**>(cur) = next;
            }
            *reinterpret_cast<void**>(chunk + (size-1)*POOL_CHUNK_SIZE) = pools_[size_in];
            pools_[size_in] = chunk;
            chunks_.push_back(chunk);
        }


    private:
        std::vector<void*> pools_;

        std::vector<Memory::pointer> chunks_;

    };

    static Pools pools_;



    inline void fast_two_sum(double a, double b, double& x, double& y) {
        x = a + b;
        double bvirt = x - a;
        y = b - bvirt;
    }

#ifdef REMOVE_ME
    inline void fast_two_diff(double a, double b, double& x, double& y) {
        x = a - b;
        double bvirt = a - x;
        y = bvirt - b;
    }
#endif

    inline void two_one_sum(
        double a1, double a0, double b, double& x2, double& x1, double& x0
    ) {
        double _i;
        two_sum(a0, b, _i, x0);
        two_sum(a1, _i, x2, x1);
    }

    inline void two_two_sum(
        double a1, double a0, double b1, double b0,
        double& x3, double& x2, double& x1, double& x0
    ) {
        double _j, _0;
        two_one_sum(a1, a0, b0, _j, _0, x0);
        two_one_sum(_j, _0, b1, x3, x2, x1);
    }

    inline void two_product_presplit(
        double a, double b, double bhi, double blo, double& x, double& y
    ) {
        x = a * b;
        double ahi;
        double alo;
        split(a, ahi, alo);
        double err1 = x - (ahi * bhi);
        double err2 = err1 - (alo * bhi);
        double err3 = err2 - (ahi * blo);
        y = (alo * blo) - err3;
    }

    inline void two_product_2presplit(
        double a, double ahi, double alo,
        double b, double bhi, double blo,
        double& x, double& y
    ) {
        x = a * b;
        double err1 = x - (ahi * bhi);
        double err2 = err1 - (alo * bhi);
        double err3 = err2 - (ahi * blo);
        y = (alo * blo) - err3;
    }

    inline void two_square(
        double a1, double a0,
        double* x
    ) {
        double _0, _1, _2;
        double _j, _k, _l;
        square(a0, _j, x[0]);
        _0 = a0 + a0;
        two_product(a1, _0, _k, _1);
        two_one_sum(_k, _1, _j, _l, _2, x[1]);
        square(a1, _j, _1);
        two_two_sum(_j, _1, _l, _2, x[5], x[4], x[3], x[2]);
    }

    void two_two_product(
        const double* a,
        const double* b,
        double* x
    ) {
        double _0, _1, _2;
        double _i, _j, _k, _l, _m, _n;

        // If the target processor supports the FMA (Fused Multiply Add)
        // instruction, then the product of two doubles into a length-2
        // expansion can be implemented as follows. Thanks to Marc Glisse
        // for the information.
        // Note: under gcc, automatic generations of fma() for a*b+c needs
        // to be deactivated, using -ffp-contract=off, else it may break
        // other functions such as fast_expansion_sum_zeroelim().
#ifdef FP_FAST_FMA
        two_product(a[0],b[0],_i,x[0]);
        two_product(a[1],b[0],_j,_0);
        two_sum(_i, _0, _k, _1);
        fast_two_sum(_j, _k, _l, _2);
        two_product(a[0], b[1], _i, _0);
        two_sum(_1, _0, _k, x[1]);
        two_sum(_2, _k, _j, _1);
        two_sum(_l, _j, _m, _2);
        two_product(a[1], b[1], _j, _0);
        two_sum(_i, _0, _n, _0);
        two_sum(_1, _0, _i, x[2]);
        two_sum(_2, _i, _k, _1);
        two_sum(_m, _k, _l, _2);
        two_sum(_j, _n, _k, _0);
        two_sum(_1, _0, _j, x[3]);
        two_sum(_2, _j, _i, _1);
        two_sum(_l, _i, _m, _2);
        two_sum(_1, _k, _i, x[4]);
        two_sum(_2, _i, _k, x[5]);
        two_sum(_m, _k, x[7], x[6]);
#else
        double a0hi, a0lo;
        split(a[0], a0hi, a0lo);
        double bhi, blo;
        split(b[0], bhi, blo);
        two_product_2presplit(
            a[0], a0hi, a0lo, b[0], bhi, blo, _i, x[0]
        );
        double a1hi, a1lo;
        split(a[1], a1hi, a1lo);
        two_product_2presplit(
            a[1], a1hi, a1lo, b[0], bhi, blo, _j, _0
        );
        two_sum(_i, _0, _k, _1);
        fast_two_sum(_j, _k, _l, _2);
        split(b[1], bhi, blo);
        two_product_2presplit(
            a[0], a0hi, a0lo, b[1], bhi, blo, _i, _0
        );
        two_sum(_1, _0, _k, x[1]);
        two_sum(_2, _k, _j, _1);
        two_sum(_l, _j, _m, _2);
        two_product_2presplit(
            a[1], a1hi, a1lo, b[1], bhi, blo, _j, _0
        );
        two_sum(_i, _0, _n, _0);
        two_sum(_1, _0, _i, x[2]);
        two_sum(_2, _i, _k, _1);
        two_sum(_m, _k, _l, _2);
        two_sum(_j, _n, _k, _0);
        two_sum(_1, _0, _j, x[3]);
        two_sum(_2, _j, _i, _1);
        two_sum(_l, _i, _m, _2);
        two_sum(_1, _k, _i, x[4]);
        two_sum(_2, _i, _k, x[5]);
        two_sum(_m, _k, x[7], x[6]);
#endif
    }

    void grow_expansion_zeroelim(
        const expansion& e, double b, expansion& h
    ) {
        double Q, hh;
        double Qnew;
        index_t eindex, hindex;
        index_t elen = e.length();

        hindex = 0;
        Q = b;
        for(eindex = 0; eindex < elen; eindex++) {
            double enow = e[eindex];
            two_sum(Q, enow, Qnew, hh);
            Q = Qnew;
            if(hh != 0.0) {
                h[hindex++] = hh;
            }
        }
        if((Q != 0.0) || (hindex == 0)) {
            h[hindex++] = Q;
        }
        h.set_length(hindex);
    }

    void scale_expansion_zeroelim(
        const expansion& e, double b, expansion& h
    ) {
        double Q, sum;
        double hh;
        double product1;
        double product0;
        index_t eindex, hindex;

        // If the target processor supports the FMA (Fused Multiply Add)
        // instruction, then the product of two doubles into a length-2
        // expansion can be implemented as follows. Thanks to Marc Glisse
        // for the information.
        // Note: under gcc, automatic generations of fma() for a*b+c needs
        // to be deactivated, using -ffp-contract=off, else it may break
        // other functions such as fast_expansion_sum_zeroelim().
#ifndef FP_FAST_FMA
        double bhi, blo;
#endif
        index_t elen = e.length();

        // Sanity check: e and h cannot be the same.
        geo_debug_assert(&e != &h);

#ifdef FP_FAST_FMA
        two_product(e[0], b, Q, hh);
#else
        split(b, bhi, blo);
        two_product_presplit(e[0], b, bhi, blo, Q, hh);
#endif

        hindex = 0;
        if(hh != 0) {
            h[hindex++] = hh;
        }
        for(eindex = 1; eindex < elen; eindex++) {
            double enow = e[eindex];
#ifdef FP_FAST_FMA
            two_product(enow, b,  product1, product0);
#else
            two_product_presplit(enow, b, bhi, blo, product1, product0);
#endif
            two_sum(Q, product0, sum, hh);
            if(hh != 0) {
                h[hindex++] = hh;
            }
            fast_two_sum(product1, sum, Q, hh);
            if(hh != 0) {
                h[hindex++] = hh;
            }
        }
        if((Q != 0.0) || (hindex == 0)) {
            h[hindex++] = Q;
        }
        h.set_length(hindex);
    }

    void fast_expansion_sum_zeroelim(
        const expansion& e, const expansion& f, expansion& h
    ) {
        double Q;
        double Qnew;
        double hh;
        index_t eindex, findex, hindex;
        double enow, fnow;
        index_t elen = e.length();
        index_t flen = f.length();

        // sanity check: h cannot be e or f
        geo_debug_assert(&h != &e);
        geo_debug_assert(&h != &f);

        enow = e[0];
        fnow = f[0];
        eindex = findex = 0;
        if((fnow > enow) == (fnow > -enow)) {
            Q = enow;
            enow = e[++eindex];
        } else {
            Q = fnow;
            fnow = f[++findex];
        }
        hindex = 0;
        if((eindex < elen) && (findex < flen)) {
            if((fnow > enow) == (fnow > -enow)) {
                fast_two_sum(enow, Q, Qnew, hh);
                enow = e[++eindex];
            } else {
                fast_two_sum(fnow, Q, Qnew, hh);
                fnow = f[++findex];
            }
            Q = Qnew;
            if(hh != 0.0) {
                h[hindex++] = hh;
            }
            while((eindex < elen) && (findex < flen)) {
                if((fnow > enow) == (fnow > -enow)) {
                    two_sum(Q, enow, Qnew, hh);
                    enow = e[++eindex];
                } else {
                    two_sum(Q, fnow, Qnew, hh);
                    fnow = f[++findex];
                }
                Q = Qnew;
                if(hh != 0.0) {
                    h[hindex++] = hh;
                }
            }
        }
        while(eindex < elen) {
            two_sum(Q, enow, Qnew, hh);
            enow = e[++eindex];
            Q = Qnew;
            if(hh != 0.0) {
                h[hindex++] = hh;
            }
        }
        while(findex < flen) {
            two_sum(Q, fnow, Qnew, hh);
            fnow = f[++findex];
            Q = Qnew;
            if(hh != 0.0) {
                h[hindex++] = hh;
            }
        }
        if((Q != 0.0) || (hindex == 0)) {
            h[hindex++] = Q;
        }
        h.set_length(hindex);
    }

    void fast_expansion_diff_zeroelim(
        const expansion& e, const expansion& f, expansion& h
    ) {
        double Q;
        double Qnew;
        double hh;
        index_t eindex, findex, hindex;
        double enow, fnow;
        index_t elen = e.length();
        index_t flen = f.length();

        // sanity check: h cannot be e or f
        geo_debug_assert(&h != &e);
        geo_debug_assert(&h != &f);

        enow = e[0];
        fnow = -f[0];
        eindex = findex = 0;
        if((fnow > enow) == (fnow > -enow)) {
            Q = enow;
            enow = e[++eindex];
        } else {
            Q = fnow;
            fnow = -f[++findex];
        }
        hindex = 0;
        if((eindex < elen) && (findex < flen)) {
            if((fnow > enow) == (fnow > -enow)) {
                fast_two_sum(enow, Q, Qnew, hh);
                enow = e[++eindex];
            } else {
                fast_two_sum(fnow, Q, Qnew, hh);
                fnow = -f[++findex];
            }
            Q = Qnew;
            if(hh != 0.0) {
                h[hindex++] = hh;
            }
            while((eindex < elen) && (findex < flen)) {
                if((fnow > enow) == (fnow > -enow)) {
                    two_sum(Q, enow, Qnew, hh);
                    enow = e[++eindex];
                } else {
                    two_sum(Q, fnow, Qnew, hh);
                    fnow = -f[++findex];
                }
                Q = Qnew;
                if(hh != 0.0) {
                    h[hindex++] = hh;
                }
            }
        }
        while(eindex < elen) {
            two_sum(Q, enow, Qnew, hh);
            enow = e[++eindex];
            Q = Qnew;
            if(hh != 0.0) {
                h[hindex++] = hh;
            }
        }
        while(findex < flen) {
            two_sum(Q, fnow, Qnew, hh);
            fnow = -f[++findex];
            Q = Qnew;
            if(hh != 0.0) {
                h[hindex++] = hh;
            }
        }
        if((Q != 0.0) || (hindex == 0)) {
            h[hindex++] = Q;
        }
        h.set_length(hindex);
    }
}



namespace GEO {

    double expansion_splitter_;
    double expansion_epsilon_;

    void expansion::initialize() {
        // Taken from Jonathan Shewchuk's exactinit.
        double half;
        double check, lastcheck;
        int every_other;

        every_other = 1;
        half = 0.5;
        expansion_epsilon_ = 1.0;
        expansion_splitter_ = 1.0;
        check = 1.0;
        // Repeatedly divide `epsilon' by two until it is too small to add to
        // one without causing roundoff.  (Also check if the sum is equal to
        // the previous sum, for machines that round up instead of using exact
        // rounding.  Not that this library will work on such machines anyway.
        do {
            lastcheck = check;
            expansion_epsilon_ *= half;
            if(every_other) {
                expansion_splitter_ *= 2.0;
            }
            every_other = !every_other;
            check = 1.0 + expansion_epsilon_;
        } while((check != 1.0) && (check != lastcheck));
        expansion_splitter_ += 1.0;
    }

    expansion* expansion::new_expansion_on_heap(index_t capa) {
        static Process::spinlock lock = 0;
        if(expansion_length_stat_) {
            Process::acquire_spinlock(lock);
            if(capa >= expansion_length_histo_.size()) {
                expansion_length_histo_.resize(capa + 1);
            }
            expansion_length_histo_[capa]++;
            Process::release_spinlock(lock);
        }
        Memory::pointer addr = Memory::pointer(
            pools_.malloc(expansion::bytes(capa))
        );
        expansion* result = new(addr)expansion(capa);
        return result;
    }

    void expansion::delete_expansion_on_heap(expansion* e) {
        pools_.free(e, expansion::bytes(e->capacity()));
    }

    // ====== Initialization from expansion and double ===============

    expansion& expansion::assign_sum(const expansion& a, double b) {
        geo_debug_assert(capacity() >= sum_capacity(a, b));
        grow_expansion_zeroelim(a, b, *this);
        return *this;
    }

    expansion& expansion::assign_diff(const expansion& a, double b) {
        geo_debug_assert(capacity() >= diff_capacity(a, b));
        grow_expansion_zeroelim(a, -b, *this);
        return *this;
    }

    expansion& expansion::assign_product(const expansion& a, double b) {
        // TODO: implement special case where the double argument
        // is a power of two.
        geo_debug_assert(capacity() >= product_capacity(a, b));
        scale_expansion_zeroelim(a, b, *this);
        return *this;
    }

    // =============  expansion sum and difference =========================

    expansion& expansion::assign_sum(
        const expansion& a, const expansion& b
    ) {
        geo_debug_assert(capacity() >= sum_capacity(a, b));
        fast_expansion_sum_zeroelim(a, b, *this);
        return *this;
    }

    expansion& expansion::assign_sum(
        const expansion& a, const expansion& b, const expansion& c
    ) {
        geo_debug_assert(capacity() >= sum_capacity(a, b, c));
        expansion& ab = expansion_sum(a, b);
        this->assign_sum(ab, c);
        return *this;
    }

    expansion& expansion::assign_sum(
        const expansion& a, const expansion& b,
        const expansion& c, const expansion& d
    ) {
        geo_debug_assert(capacity() >= sum_capacity(a, b, c));
        expansion& ab = expansion_sum(a, b);
        expansion& cd = expansion_sum(c, d);
        this->assign_sum(ab, cd);
        return *this;
    }

    expansion& expansion::assign_diff(const expansion& a, const expansion& b) {
        geo_debug_assert(capacity() >= diff_capacity(a, b));
        fast_expansion_diff_zeroelim(a, b, *this);
        return *this;
    }

    // =============  expansion product ==================================

    // Recursive helper function for product implementation
    expansion& expansion::assign_sub_product(
        const double* a, index_t a_length, const expansion& b
    ) {
        geo_debug_assert(
            capacity() >= sub_product_capacity(a_length, b.length())
        );
        if(a_length == 1) {
            scale_expansion_zeroelim(b, a[0], *this);
        } else {
            // "Distillation" (see Shewchuk's paper) is computed recursively,
            // by splitting the list of expansions to sum into two halves.
            const double* a1 = a;
            index_t a1_length = a_length / 2;
            const double* a2 = a1 + a1_length;
            index_t a2_length = a_length - a1_length;
            expansion& a1b = expansion_sub_product(a1, a1_length, b);
            expansion& a2b = expansion_sub_product(a2, a2_length, b);
            this->assign_sum(a1b, a2b);
        }
        return *this;
    }

    expansion& expansion::assign_product(
        const expansion& a, const expansion& b
    ) {
        geo_debug_assert(capacity() >= product_capacity(a, b));
        if(a.length() == 0 || b.length() == 0) {
            x_[0] = 0.0;
            set_length(0);
        } else if(a.length() == 1 && b.length() == 1) {
            two_product(a[0], b[0], x_[1], x_[0]);
            set_length(2);
        } else if(a.length() == 1) {
            scale_expansion_zeroelim(b, a[0], *this);
        } else if(b.length() == 1) {
            scale_expansion_zeroelim(a, b[0], *this);
        } else if(a.length() == 2 && b.length() == 2) {
            two_two_product(a.data(), b.data(), x_);
            set_length(8);
        } else {
            // Recursive distillation: the shortest expansion
            // is split into two parts.
            if(a.length() < b.length()) {
                const double* a1 = a.data();
                index_t a1_length = a.length() / 2;
                const double* a2 = a1 + a1_length;
                index_t a2_length = a.length() - a1_length;
                expansion& a1b = expansion_sub_product(a1, a1_length, b);
                expansion& a2b = expansion_sub_product(a2, a2_length, b);
                this->assign_sum(a1b, a2b);
            } else {
                const double* b1 = b.data();
                index_t b1_length = b.length() / 2;
                const double* b2 = b1 + b1_length;
                index_t b2_length = b.length() - b1_length;
                expansion& ab1 = expansion_sub_product(b1, b1_length, a);
                expansion& ab2 = expansion_sub_product(b2, b2_length, a);
                this->assign_sum(ab1, ab2);
            }
        }
        return *this;
    }

    expansion& expansion::assign_product(
        const expansion& a, const expansion& b, const expansion& c
    ) {
        const expansion& bc = expansion_product(b, c);
        this->assign_product(a, bc);
        return *this;
    }

    expansion& expansion::assign_square(const expansion& a) {
        geo_debug_assert(capacity() >= square_capacity(a));
        if(a.length() == 1) {
            square(a[0], x_[1], x_[0]);
            set_length(2);
        } else if(a.length() == 2) {
            two_square(a[1], a[0], x_);
            set_length(6);
        } else {
            this->assign_product(a, a);
        }
        return *this;
    }

    // =============  determinants ==========================================

    expansion& expansion::assign_det2x2(
        const expansion& a11, const expansion& a12,
        const expansion& a21, const expansion& a22
    ) {
        const expansion& a11a22 = expansion_product(a11, a22);
        const expansion& a12a21 = expansion_product(a12, a21);
        return this->assign_diff(a11a22, a12a21);
    }

    expansion& expansion::assign_det3x3(
        const expansion& a11, const expansion& a12, const expansion& a13,
        const expansion& a21, const expansion& a22, const expansion& a23,
        const expansion& a31, const expansion& a32, const expansion& a33
    ) {
        // Development w.r.t. first row
        const expansion& c11 = expansion_det2x2(a22, a23, a32, a33);
        const expansion& c12 = expansion_det2x2(a23, a21, a33, a31);
        const expansion& c13 = expansion_det2x2(a21, a22, a31, a32);
        const expansion& a11c11 = expansion_product(a11, c11);
        const expansion& a12c12 = expansion_product(a12, c12);
        const expansion& a13c13 = expansion_product(a13, c13);
        return this->assign_sum(a11c11, a12c12, a13c13);
    }

    expansion& expansion::assign_det_111_2x3(
        const expansion& a21, const expansion& a22, const expansion& a23,
        const expansion& a31, const expansion& a32, const expansion& a33
    ) {
        const expansion& c11 = expansion_det2x2(a22, a23, a32, a33);
        const expansion& c12 = expansion_det2x2(a23, a21, a33, a31);
        const expansion& c13 = expansion_det2x2(a21, a22, a31, a32);
        return this->assign_sum(c11, c12, c13);
    }

    // =============  geometric operations ==================================

    expansion& expansion::assign_sq_dist(
        const double* p1, const double* p2, coord_index_t dim
    ) {
        geo_debug_assert(capacity() >= sq_dist_capacity(dim));
        if(dim == 1) {
            double d0, d1;
            two_diff(p1[0], p2[0], d1, d0);
            two_square(d1, d0, x_);
            set_length(6);
        } else {
            // "Distillation" (see Shewchuk's paper) is computed recursively,
            // by splitting the list of expansions to sum into two halves.
            coord_index_t dim1 = dim / 2;
            coord_index_t dim2 = coord_index_t(dim - dim1);
            const double* p1_2 = p1 + dim1;
            const double* p2_2 = p2 + dim1;
            expansion& d1 = expansion_sq_dist(p1, p2, dim1);
            expansion& d2 = expansion_sq_dist(p1_2, p2_2, dim2);
            this->assign_sum(d1, d2);
        }
        return *this;
    }

    expansion& expansion::assign_dot_at(
        const double* p1, const double* p2, const double* p0,
        coord_index_t dim
    ) {
        geo_debug_assert(capacity() >= dot_at_capacity(dim));
        if(dim == 1) {

            double v[2];
            two_diff(p1[0], p0[0], v[1], v[0]);
            double w[2];
            two_diff(p2[0], p0[0], w[1], w[0]);
            two_two_product(v, w, x_);
            set_length(8);
        } else {
            // "Distillation" (see Shewchuk's paper) is computed recursively,
            // by splitting the list of expansions to sum into two halves.
            coord_index_t dim1 = dim / 2;
            coord_index_t dim2 = coord_index_t(dim - dim1);
            const double* p1_2 = p1 + dim1;
            const double* p2_2 = p2 + dim1;
            const double* p0_2 = p0 + dim1;
            expansion& d1 = expansion_dot_at(p1, p2, p0, dim1);
            expansion& d2 = expansion_dot_at(p1_2, p2_2, p0_2, dim2);
            this->assign_sum(d1, d2);
        }
        return *this;
    }

    expansion& expansion::assign_length2(
        const expansion& x, const expansion& y, const expansion& z
    ) {
        const expansion& x2 = expansion_square(x);
        const expansion& y2 = expansion_square(y);
        const expansion& z2 = expansion_square(z);
        this->assign_sum(x2,y2,z2);
        return *this;
    }



    Sign sign_of_expansion_determinant(
        const expansion& a00,const expansion& a01,
        const expansion& a10,const expansion& a11
    ) {
        const expansion& result = expansion_det2x2(a00, a01, a10, a11);
        return result.sign();
    }

    Sign sign_of_expansion_determinant(
        const expansion& a00,const expansion& a01,const expansion& a02,
        const expansion& a10,const expansion& a11,const expansion& a12,
        const expansion& a20,const expansion& a21,const expansion& a22
    ) {
        // First compute the det2x2
        const expansion& m01 =
            expansion_det2x2(a00, a10, a01, a11);
        const expansion& m02 =
            expansion_det2x2(a00, a20, a01, a21);
        const expansion& m12 =
            expansion_det2x2(a10, a20, a11, a21);

        // Now compute the minors of rank 3
        const expansion& z1 = expansion_product(m01,a22);
        const expansion& z2 = expansion_product(m02,a12).negate();
        const expansion& z3 = expansion_product(m12,a02);

        const expansion& result = expansion_sum3(z1,z2,z3);
        return result.sign();
    }

    Sign sign_of_expansion_determinant(
        const expansion& a00,const expansion& a01,
        const expansion& a02,const expansion& a03,
        const expansion& a10,const expansion& a11,
        const expansion& a12,const expansion& a13,
        const expansion& a20,const expansion& a21,
        const expansion& a22,const expansion& a23,
        const expansion& a30,const expansion& a31,
        const expansion& a32,const expansion& a33
    ) {

        // First compute the det2x2
        const expansion& m01 =
            expansion_det2x2(a10,a00,a11,a01);
        const expansion& m02 =
            expansion_det2x2(a20,a00,a21,a01);
        const expansion& m03 =
            expansion_det2x2(a30,a00,a31,a01);
        const expansion& m12 =
            expansion_det2x2(a20,a10,a21,a11);
        const expansion& m13 =
            expansion_det2x2(a30,a10,a31,a11);
        const expansion& m23 =
            expansion_det2x2(a30,a20,a31,a21);

        // Now compute the minors of rank 3
        const expansion& m012_1 = expansion_product(m12,a02);
        expansion& m012_2 = expansion_product(m02,a12); m012_2.negate();
        const expansion& m012_3 = expansion_product(m01,a22);
        const expansion& m012 = expansion_sum3(m012_1, m012_2, m012_3);

        const expansion& m013_1 = expansion_product(m13,a02);
        expansion& m013_2 = expansion_product(m03,a12); m013_2.negate();

        const expansion& m013_3 = expansion_product(m01,a32);
        const expansion& m013 = expansion_sum3(m013_1, m013_2, m013_3);

        const expansion& m023_1 = expansion_product(m23,a02);
        expansion& m023_2 = expansion_product(m03,a22); m023_2.negate();
        const expansion& m023_3 = expansion_product(m02,a32);
        const expansion& m023 = expansion_sum3(m023_1, m023_2, m023_3);

        const expansion& m123_1 = expansion_product(m23,a12);
        expansion& m123_2 = expansion_product(m13,a22); m123_2.negate();
        const expansion& m123_3 = expansion_product(m12,a32);
        const expansion& m123 = expansion_sum3(m123_1, m123_2, m123_3);

        // Now compute the minors of rank 4
        const expansion& m0123_1 = expansion_product(m123,a03);
        const expansion& m0123_2 = expansion_product(m023,a13);
        const expansion& m0123_3 = expansion_product(m013,a23);
        const expansion& m0123_4 = expansion_product(m012,a33);

        const expansion& z1 = expansion_sum(m0123_1, m0123_3);
        const expansion& z2 = expansion_sum(m0123_2, m0123_4);

        const expansion& result = expansion_diff(z1,z2);
        return result.sign();
    }



}


namespace GEO {

    namespace PCK {

      // filter predicates
      #include "numerics/predicates/side_filters.h"

      // exact predicates
      #include "numerics/predicates/side_exact.h"

      #include "numerics/predicates/side_impl.h"

      void initialize() {
        expansion::initialize();
      }

    }
}
