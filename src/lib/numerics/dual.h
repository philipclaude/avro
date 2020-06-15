//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef AVRO_NUMERICS_DUAL_H_
#define AVRO_NUMERICS_DUAL_H_

#include "common/tools.h"
#include "common/error.h"

#include <stdio.h>
#include <cmath>
#include <iostream>

class dual {

	public:

		/* Constructors and destructors */
		dual() {u=0.;v=0.;}
		dual(double u0) {u=u0;v=0.;}
		dual(double u0,double v0) {u=u0;v=v0;}
		~dual() {}

		/* Variables */
		double u,v;

		/* Basic functions */
		void   print();
		double re();
		double du();

		/* Assignment operator */
		dual& operator= (double r);

		/* Addition operators */
		dual operator+ () const;
		dual operator+ (const dual d) const;
    dual operator+= (const dual d);
		friend dual operator+ (const double r, const dual d);
		friend dual operator+ (const dual d, const double r);

		/* Subtraction operators */
		dual operator- () const;
		dual operator- (const dual d) const;
    dual operator-= (const dual d);
		friend dual operator- (const double r, const dual d);
		friend dual operator- (const dual d, const double r);

		/* Multiplication operators */
		dual operator* (const dual d) const;
    dual operator*= (const dual d);
		friend dual operator* (const double r, const dual d);
		friend dual operator* (const dual d, const double r);

		/* Division operators */
		friend dual operator/ (const dual d1, const dual d2);
		friend dual operator/ (const double r, const dual d);
		friend dual operator/ (const dual d, const double r);

		/* Math functions */
		friend dual pow(const dual &x, const double &p);
		friend dual pow(const dual x, const dual y);
		friend dual exp(const dual x);
		friend dual log(const dual x);
		friend dual sqrt(const dual &x);
		friend dual sin(const dual x);
		friend dual cos(const dual x);
		friend dual tan(const dual x);
		friend double fabs(const dual x);

		/* Greater than */
		friend bool operator> (dual d1,dual d2);
		friend bool operator> (double r,dual d);
		friend bool operator> (dual d,double r);

		/* Less than */
		friend bool operator< (dual d1,dual d2);
		friend bool operator< (double r,dual d);
		friend bool operator< (dual d,double r);

    friend bool operator==( dual d1, dual d2) { return (d1.u==d2.u && d1.v==d2.v); }
    friend bool operator==( dual d, double x ) { return d.u==x; }

		/* I/O functions */
		friend std::ostream& operator<< (std::ostream &output, const dual &d);
};

/* Basic functions */
inline void
dual::print() {
	printf("(%f,%f)\n",u,v);
}

inline double
dual::re() {
	return u;
}

inline double
dual::du() {
	return v;
}

/* Assignment operator */
inline dual&
dual::operator= (double r) {
	u = r;
	v = 0.;
	return *this;
}

/* Addition operators */
inline dual
dual::operator+ () const {
	return *this;
}

inline dual
dual::operator+ (const dual d) const {
	return dual(u+d.u,v+d.v);
}

inline dual
dual::operator+= ( const dual d) {
  u += d.u;
  v += d.v;
  return *this;
}

inline dual
operator+ (const double r, const dual d) {
	return dual(r+d.u,d.v);
}

inline dual
operator+ (const dual d, const double r) {
	return dual(r+d.u,d.v);
}

/* Subtraction operators */
inline dual
dual::operator- () const {
	return dual(-u,-v);
}

inline dual
dual::operator- (const dual d) const {
	return dual(u-d.u,v-d.v);
}

inline dual
dual::operator-= ( const dual d) {
  u -= d.u;
  v -= d.v;
  return *this;
}

inline dual
operator- (const double r, const dual d) {
	return dual(r-d.u,-d.v);
}

inline dual
operator- (const dual d, const double r) {
	return dual(d.u-r,d.v);
}

/* Multiplication operators */
inline dual
dual::operator* (const dual d) const {
	return dual(u*d.u,v*d.u+u*d.v);
}

inline dual
dual::operator*= ( const dual d) {
  u *= d.u;
  v *= d.v;
  return *this;
}

inline dual
operator* (const double r, const dual d) {
	return dual(r*d.u,r*d.v);
}

inline dual
operator* (const dual d, const double r) {
	return dual(r*d.u,r*d.v);
}

/* Division operators */
inline dual
operator/ (const dual d1, const dual d2) {
	if (fabs(d2.u)<1e-14) {
		std::cout << d1 << " / " << d2 << " yields division by zero.\n";
		avro_assert(false);
	}
	return dual(d1.u/d2.u,d1.v/d2.u-d1.u*d2.v/(d2.u*d2.u));
}

inline dual
operator/ (const double r, const dual d) {
	return dual(r/d.u,-r*d.v/(d.u*d.u));
}

inline dual
operator/ (const dual d, const double r) {
	return dual(d.u/r,d.v/r);
}

/* Math functions */
inline dual
pow(const dual &x, const double &p) {
	return dual(std::pow(x.u,p),x.v*p*std::pow(x.u,p-1));
}

inline dual
pow(const dual x, const dual y) {
	return dual(std::pow(x.u,y.u),x.v*y.u*std::pow(x.u,y.u-1)+y.v*std::pow(x.u,y.u)*log(x.u));
}

inline dual
exp(const dual x) {
	return dual(exp(x.u),exp(x.u)*x.v);
}

inline dual
log(const dual x) {
	return dual(log(x.u),x.v/x.u);
}

inline dual
sin(const dual x) {
	return dual(sin(x.u),cos(x.u)*x.v);
}

inline dual
cos(const dual x) {
	return dual(cos(x.u),-sin(x.u)*x.v);
}

inline dual
tan(const dual x) {
	return sin(x)/cos(x);
}

inline double
fabs(const dual x) {
	return fabs(x.u);
}

inline dual
sqrt(const dual &x) {
	return dual(std::sqrt(x.u),.5*x.v/std::sqrt(x.u));
}

/* Greater than */
inline bool
operator> (dual d1,dual d2) {
	return d1.u>d2.u;
}

inline bool
operator> (double r,dual d) {
	return r>d.u;
}

inline bool
operator> (dual d,double r) {
	return d.u>r;
}

/* Less than */
inline bool
operator< (dual d1,dual d2) {
	return d1.u<d2.u;
}

inline bool
operator< (double r,dual d) {
	return r<d.u;
}

inline bool
operator< (dual d,double r) {
	return d.u<r;
}

/* I/O functions */
inline std::ostream&
operator<< (std::ostream &output, const dual &d) {
	output << "(" << d.u << ","<< d.v << ")";
	return output;
}

#endif
