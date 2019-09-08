// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define MATRIXS_INVERSELU_INSTANTIATE
#include "MatrixS_InverseLU_impl.h"

//#include <boost/preprocessor/repetition/repeat_from_to.hpp>

namespace SANS
{
namespace DLA
{

// Explicitly instantiate all datatypes used here
template struct MatrixSLUSolver< 2, 2, Real, MatrixS<2,1,Real> >;
template struct MatrixSLUSolver< 2, 2, Real, MatrixS<2,3,Real> >;

template struct MatrixSLUSolver< 3, 3, Real, MatrixS<3,1,Real> >;
template struct MatrixSLUSolver< 4, 4, Real, MatrixS<4,1,Real> >;
template struct MatrixSLUSolver< 5, 5, Real, MatrixS<5,1,Real> >;
template struct MatrixSLUSolver< 6, 6, Real, MatrixS<6,1,Real> >;

template struct MatrixSLUSolver< 2, 2, MatrixS<2,2,Real>, MatrixS<2,2,MatrixS<2,2,Real>> >;
template struct MatrixSLUSolver< 2, 2, MatrixS<3,3,Real>, MatrixS<2,2,MatrixS<3,3,Real>> >;
template struct MatrixSLUSolver< 2, 2, MatrixS<4,4,Real>, MatrixS<2,2,MatrixS<4,4,Real>> >;
template struct MatrixSLUSolver< 2, 2, MatrixS<5,5,Real>, MatrixS<2,2,MatrixS<5,5,Real>> >;

template struct MatrixSLUSolver< 3, 3, MatrixS<2,2,Real>, MatrixS<3,3,MatrixS<2,2,Real>> >;
template struct MatrixSLUSolver< 3, 3, MatrixS<3,3,Real>, MatrixS<3,3,MatrixS<3,3,Real>> >;
template struct MatrixSLUSolver< 3, 3, MatrixS<4,4,Real>, MatrixS<3,3,MatrixS<4,4,Real>> >;
template struct MatrixSLUSolver< 3, 3, MatrixS<5,5,Real>, MatrixS<3,3,MatrixS<5,5,Real>> >;
template struct MatrixSLUSolver< 3, 3, MatrixS<6,6,Real>, MatrixS<3,3,MatrixS<6,6,Real>> >;

#define DECL(z, n, text) template struct MatrixSLUSolver< n, n, Real, MatrixS<n,n,Real> >;
//BOOST_PP_REPEAT_FROM_TO(1, 16, DECL, )

DECL( ,1 , )
DECL( ,2 , )
DECL( ,3 , )
DECL( ,4 , )
DECL( ,5 , )
DECL( ,6 , )
DECL( ,7 , )
DECL( ,8 , )
DECL( ,9 , )
DECL( ,10 , )
DECL( ,11 , )
DECL( ,12 , )
DECL( ,13 , )
DECL( ,14 , )
DECL( ,15 , )
//DECL( ,16 , )

}
}
