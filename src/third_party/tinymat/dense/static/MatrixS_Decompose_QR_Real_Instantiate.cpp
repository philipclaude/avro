// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define MATRIXS_DECOMPOSE_QR_INSTANTIATE
#include "MatrixS_Decompose_QR_impl.h"

//#include <boost/preprocessor/repetition/repeat_from_to.hpp>

namespace tinymat
{
namespace DLA
{

template struct MatrixSQR<3,2,Real>;

//Explicit instantiation
#define DECL(z, n, text) template struct MatrixSQR<n,n,Real>;
//BOOST_PP_REPEAT_FROM_TO(1, 8, DECL, )

DECL( ,1 , )
DECL( ,2 , )
DECL( ,3 , )
DECL( ,4 , )
DECL( ,5 , )
DECL( ,6 , )
DECL( ,7 , )
//DECL( ,8 , )

}
}
