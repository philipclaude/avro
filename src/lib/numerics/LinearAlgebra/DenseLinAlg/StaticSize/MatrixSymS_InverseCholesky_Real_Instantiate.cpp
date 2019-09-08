// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define MATRIXSYMS_INVERSECHOLESKY_INSTANTIATE
#include "MatrixSymS_InverseCholesky_impl.h"

//#include <boost/preprocessor/repetition/repeat_from_to.hpp>

namespace SANS
{
namespace DLA
{

// Explicitly instantiate all datatypes used here
template struct MatrixSymSCholeskySolver< 2, 2, Real, MatrixS<2,1,Real> >;
template struct MatrixSymSCholeskySolver< 2, 2, Real, MatrixS<2,3,Real> >;

template struct MatrixSymSCholeskySolver< 3, 3, Real, MatrixS<3,1,Real> >;

#define DECL(z, n, text) template struct MatrixSymSCholeskySolver< n, n, Real, MatrixS<n,n,Real> >;
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
#undef DECL

#define DECL(z, n, text) template struct MatrixSymSCholeskySolver< n, n, Real, MatrixSymS<n,Real> >;
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

#undef DECL

}
}
