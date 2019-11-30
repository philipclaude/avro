#ifndef LUNA_NUMERICS_EXPANSION_H_
#define LUNA_NUMERICS_EXPANSION_H_

#include "numerics/predicates.h"

#include <algorithm>
#include <new>

namespace GEO {

  extern double expansion_splitter_;
  extern double expansion_epsilon_;

  inline void two_sum(double a, double b, double& x, double& y) {
      x = a + b;
      double bvirt = x - a;
      double avirt = x - bvirt;
      double bround = b - bvirt;
      double around = a - avirt;
      y = around + bround;
  }

  inline void two_diff(double a, double b, double& x, double& y) {
      x = a - b;
      double bvirt = a - x;
      double avirt = x + bvirt;
      double bround = bvirt - b;
      double around = a - avirt;
      y = around + bround;
  }

  inline void split(double a, double& ahi, double& alo) {
      double c = expansion_splitter_ * a;
      double abig = c - a;
      ahi = c - abig;
      alo = a - ahi;
  }

  inline void two_product(double a, double b, double& x, double& y) {
#ifdef FP_FAST_FMA
      // If the target processor supports the FMA (Fused Multiply Add)
      // instruction, then the product of two doubles into a length-2
      // expansion can be implemented as follows. Thanks to Marc Glisse
      // for the information.
      // Note: under gcc, automatic generations of fma() for a*b+c needs
      // to be deactivated, using -ffp-contract=off, else it may break
      // other functions such as fast_expansion_sum_zeroelim().
      x = a*b;
      y = fma(a,b,-x);
#else
      x = a * b;
      double ahi, alo;
      split(a, ahi, alo);
      double bhi, blo;
      split(b, bhi, blo);
      double err1 = x - (ahi * bhi);
      double err2 = err1 - (alo * bhi);
      double err3 = err2 - (ahi * blo);
      y = (alo * blo) - err3;
#endif
  }

  inline void square(double a, double& x, double& y) {
#ifdef FP_FAST_FMA
      // If the target processor supports the FMA (Fused Multiply Add)
      // instruction, then the product of two doubles into a length-2
      // expansion can be implemented as follows. Thanks to Marc Glisse
      // for the information.
      // Note: under gcc, automatic generations of fma() for a*b+c needs
      // to be deactivated, using -ffp-contract=off, else it may break
      // other functions such as fast_expansion_sum_zeroelim().
      x = a*a;
      y = fma(a,a,-x);
#else
      x = a * a;
      double ahi, alo;
      split(a, ahi, alo);
      double err1 = x - (ahi * ahi);
      double err3 = err1 - ((ahi + ahi) * alo);
      y = (alo * alo) - err3;
#endif
  }


class GEOGRAM_API expansion {

public:
    index_t length() const {
        return length_;
    }

    index_t capacity() const {
        return capacity_;
    }

    void set_length(index_t new_length) {
        geo_debug_assert(new_length <= capacity());
        length_ = new_length;
    }

    const double& operator[] (index_t i) const {
        // Note: we allocate capacity+1 storage
        // systematically, since basic functions
        // may access one additional value (without
        // using it)
        geo_debug_assert(i <= capacity_);
        return x_[i];
    }

    double& operator[] (index_t i) {
        // Note: we allocate capacity+1 storage
        // systematically, since basic functions
        // may access one additional value (without
        // using it)
        geo_debug_assert(i <= capacity_);
        return x_[i];
    }

    double* data() {
        return x_;
    }

    const double* data() const {
        return x_;
    }

    static size_t bytes(index_t capa) {
        // --> 2*sizeof(double) because x_ is declared of size [2]
        // to avoid compiler's warning.
        // --> capa+1 to have an additional 'sentry' at the end
        // because fast_expansion_sum_zeroelim() may access
        // an entry past the end (without using it).
        return
            sizeof(expansion) - 2 * sizeof(double) +
            (capa + 1) * sizeof(double);
    }

    expansion(index_t capa) :
        length_(0),
        capacity_(capa) {
    }

#ifdef CPPCHECK
    // cppcheck does not understand that the result
    // of alloca() is passed to the placement syntax
    // of operator new.
expansion& new_expansion_on_stack(index_t capa);
#else
#define new_expansion_on_stack(capa)                           \
(new (alloca(expansion::bytes(capa)))expansion(capa))
//(new expansion(capa))
#endif

    static expansion* new_expansion_on_heap(index_t capa);

    static void delete_expansion_on_heap(expansion* e);

    // ========================== Initialization from doubles

    static index_t sum_capacity(double a, double b) {
        geo_argused(a);
        geo_argused(b);
        return 2;
    }

    expansion& assign_sum(double a, double b) {
        set_length(2);
        two_sum(a, b, x_[1], x_[0]);
        return *this;
    }

    static index_t diff_capacity(double a, double b) {
        geo_argused(a);
        geo_argused(b);
        return 2;
    }

    expansion& assign_diff(double a, double b) {
        set_length(2);
        two_diff(a, b, x_[1], x_[0]);
        return *this;
    }

    static index_t product_capacity(double a, double b) {
        geo_argused(a);
        geo_argused(b);
        return 2;
    }

    expansion& assign_product(double a, double b) {
        set_length(2);
        two_product(a, b, x_[1], x_[0]);
        return *this;
    }

    static index_t square_capacity(double a) {
        geo_argused(a);
        return 2;
    }

    expansion& assign_square(double a) {
        set_length(2);
        square(a, x_[1], x_[0]);
        return *this;
    }

    // ====== Initialization from expansion and double

    static index_t sum_capacity(const expansion& a, double b) {
        geo_argused(b);
        return a.length() + 1;
    }

    expansion& assign_sum(const expansion& a, double b);

    static index_t diff_capacity(const expansion& a, double b) {
        geo_argused(b);
        return a.length() + 1;
    }

    expansion& assign_diff(const expansion& a, double b);

    static index_t product_capacity(const expansion& a, double b) {
        geo_argused(b);
        // TODO: implement special case where the double argument
        // is a power of two.
        return a.length() * 2;
    }

    expansion& assign_product(const expansion& a, double b);

    // ========================== Initialization from expansions

    static index_t sum_capacity(const expansion& a, const expansion& b) {
        return a.length() + b.length();
    }

    expansion& assign_sum(const expansion& a, const expansion& b);

    static index_t sum_capacity(
        const expansion& a, const expansion& b, const expansion& c
    ) {
        return a.length() + b.length() + c.length();
    }

    expansion& assign_sum(
        const expansion& a, const expansion& b, const expansion& c
    );

    static index_t sum_capacity(
        const expansion& a, const expansion& b,
        const expansion& c, const expansion& d
    ) {
        return a.length() + b.length() + c.length() + d.length();
    }

    expansion& assign_sum(
        const expansion& a, const expansion& b,
        const expansion& c, const expansion& d
    );

    static index_t diff_capacity(const expansion& a, const expansion& b) {
        return a.length() + b.length();
    }

    expansion& assign_diff(const expansion& a, const expansion& b);

    static index_t product_capacity(
        const expansion& a, const expansion& b
    ) {
        return a.length() * b.length() * 2;
    }

    expansion& assign_product(const expansion& a, const expansion& b);

    static index_t product_capacity(
        const expansion& a, const expansion& b, const expansion& c
    ) {
        return a.length() * b.length() * c.length() * 4;
    }

    expansion& assign_product(
        const expansion& a, const expansion& b, const expansion& c
    );

    static index_t square_capacity(const expansion& a) {
        if(a.length() == 2) {
            return 6;
        }                                  // see two_square()
        return a.length() * a.length() * 2;
    }

    expansion& assign_square(const expansion& a);

    // ====== Determinants =============================

    static index_t det2x2_capacity(
        const expansion& a11, const expansion& a12,
        const expansion& a21, const expansion& a22
    ) {
        return
            product_capacity(a11, a22) +
            product_capacity(a21, a12);
    }

    expansion& assign_det2x2(
        const expansion& a11, const expansion& a12,
        const expansion& a21, const expansion& a22
    );

    static index_t det3x3_capacity(
        const expansion& a11, const expansion& a12, const expansion& a13,
        const expansion& a21, const expansion& a22, const expansion& a23,
        const expansion& a31, const expansion& a32, const expansion& a33
    ) {
        // Development w.r.t. first row
        index_t c11_capa = det2x2_capacity(a22, a23, a32, a33);
        index_t c12_capa = det2x2_capacity(a21, a23, a31, a33);
        index_t c13_capa = det2x2_capacity(a21, a22, a31, a32);
        return 2 * (
            a11.length() * c11_capa +
            a12.length() * c12_capa +
            a13.length() * c13_capa
        );
    }

    expansion& assign_det3x3(
        const expansion& a11, const expansion& a12, const expansion& a13,
        const expansion& a21, const expansion& a22, const expansion& a23,
        const expansion& a31, const expansion& a32, const expansion& a33
    );

    static index_t det_111_2x3_capacity(
        const expansion& a21, const expansion& a22, const expansion& a23,
        const expansion& a31, const expansion& a32, const expansion& a33
    ) {
        return
            det2x2_capacity(a22, a23, a32, a33) +
            det2x2_capacity(a23, a21, a33, a31) +
            det2x2_capacity(a21, a22, a31, a32);
    }

    expansion& assign_det_111_2x3(
        const expansion& a21, const expansion& a22, const expansion& a23,
        const expansion& a31, const expansion& a32, const expansion& a33
    );

    // ======= Geometry-specific initializations =======

    static index_t sq_dist_capacity(coord_index_t dim) {
        return index_t(dim) * 6;
    }

    expansion& assign_sq_dist(
        const double* p1, const double* p2, coord_index_t dim
    );

    static index_t dot_at_capacity(coord_index_t dim) {
        return index_t(dim) * 8;
    }

    expansion& assign_dot_at(
        const double* p1, const double* p2, const double* p0,
        coord_index_t dim
    );


    static index_t length2_capacity(
        const expansion& x, const expansion& y, const expansion& z
    ) {
        return square_capacity(x) + square_capacity(y) + square_capacity(z);
    }

    expansion& assign_length2(
        const expansion& x, const expansion& y, const expansion& z
    );

    // =============== some general purpose functions =========

    static void initialize();

    expansion& negate() {
        for(index_t i = 0; i < length_; ++i) {
            x_[i] = -x_[i];
        }
        return *this;
    }

    expansion& scale_fast(double s) {
        // TODO: debug assert is_power_of_two(s)
        for(index_t i = 0; i < length_; ++i) {
            x_[i] *= s;
        }
        return *this;
    }

    double estimate() const {
        double result = 0.0;
        for(index_t i = 0; i < length(); ++i) {
            result += x_[i];
        }
        return result;
    }

    Sign sign() const {
        if(length() == 0) {
            return ZERO;
        }
        return geo_sgn(x_[length() - 1]);
    }

    double value() const {
      if (length()==0) return 0.0;
      return x_[length() - 1];
    }

    std::ostream& show(std::ostream& os) const {
        for(index_t i = 0; i < length(); ++i) {
            os << i << ':' << x_[i] << ' ';
        }
        return os << std::endl;
    }

protected:
    static index_t sub_product_capacity(
        index_t a_length, index_t b_length
    ) {
        return a_length * b_length * 2;
    }

    expansion& assign_sub_product(
        const double* a, index_t a_length, const expansion& b
    );

#define expansion_sub_product(a, a_length, b)           \
new_expansion_on_stack(                       \
    sub_product_capacity(a_length, b.length()) \
)->assign_sub_product(a, a_length, b)

private:
    expansion(const expansion& rhs);

    expansion& operator= (const expansion& rhs);

private:
    index_t length_;
    index_t capacity_;
    double x_[2];  // x_ is in fact of size [capacity_]

    friend class expansion_nt;
};

} // luna

#endif
