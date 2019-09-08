// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define MATRIXSYMS_DECOMPOSE_CHOLESKY_INSTANTIATE
#include "MatrixSymS_Decompose_Cholesky_impl.h"

#include "VectorS.h"

//#include <boost/preprocessor/repetition/repeat_from_to.hpp>

namespace SANS
{
namespace DLA
{

//Explicit instantiation
#define DECL(z, n, text) template struct MatrixSymSCholesky<n,Real>;
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
