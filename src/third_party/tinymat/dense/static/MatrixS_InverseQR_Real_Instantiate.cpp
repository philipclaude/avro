// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define MATRIXS_INVERSEQR_INSTANTIATE
#include "MatrixS_InverseQR_impl.h"

//#include <boost/preprocessor/repetition/repeat_from_to.hpp>

namespace tinymat 
{
namespace DLA
{

INSTANTIATE(2, 2, 1, Real, MATRIXS(2, 1, Real) )
INSTANTIATE(2, 2, 3, Real, MATRIXS(2, 3, Real) )

INSTANTIATE(3, 2, 1, Real, MATRIXS(2, 1, Real) )
INSTANTIATE(3, 3, 1, Real, MATRIXS(3, 1, Real) )

#define DECL(z, n, text) INSTANTIATE(n, n, n, Real, MATRIXS(n, n, Real) )
//BOOST_PP_REPEAT_FROM_TO(1, 8, DECL, )

DECL( ,1, )
DECL( ,2, )
DECL( ,3, )
DECL( ,4, )
DECL( ,5, )
DECL( ,6, )
DECL( ,7, )
//DECL( ,8, )

}
}
