// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef block_TYPE_H
#define block_TYPE_H

namespace numpack 
{
namespace BLA
{

// Base class of a generic block linear algebraic expression
template< class Derived>
class blockType
{
public:
  //A convenient method for casting to the derived type
  inline       Derived& cast()       { return static_cast<      Derived&>(*this); }
  inline const Derived& cast() const { return static_cast<const Derived&>(*this); }
};

// Base class of a block matrix
template< class Derived>
class BlockMatrixType
{
public:
  //A convenient method for casting to the derived type
  inline       Derived& cast()       { return static_cast<      Derived&>(*this); }
  inline const Derived& cast() const { return static_cast<const Derived&>(*this); }

  int m() const { return cast().m(); }
  int n() const { return cast().n(); }
};

// Base class of a block vector
template< class Derived>
class BlockVectorType
{
public:
  //A convenient method for casting to the derived type
  inline       Derived& cast()       { return static_cast<      Derived&>(*this); }
  inline const Derived& cast() const { return static_cast<const Derived&>(*this); }

  int m() const { return cast().m(); }
};

template<class T>
struct isBlockMatrix { static const bool value = false; };

template<class T>
struct isBlockVector { static const bool value = false; };

// Forward declarations
//---------------------------------------------------------------------------//
template<class Matrix00, class Matrix01,
         class Matrix10, class Matrix11>
class MatrixBlock_2x2;

template<class Vector0,
         class Vector1>
class VectorBlock_2;

template<class Matrix00, class Matrix01,
         class Matrix10, class Matrix11>
struct isBlockMatrix<MatrixBlock_2x2<Matrix00, Matrix01,
                                     Matrix10, Matrix11>> { static const bool value = true; };

template<class Vector0,
         class Vector1>
struct isBlockVector<VectorBlock_2<Vector0,Vector1>> { static const bool value = true; };

//---------------------------------------------------------------------------//
template<class Matrix00, class Matrix01, class Matrix02,
         class Matrix10, class Matrix11, class Matrix12,
         class Matrix20, class Matrix21, class Matrix22>
class MatrixBlock_3x3;

template<class Vector0,
         class Vector1,
         class Vector2>
class VectorBlock_3;

template<class Matrix00, class Matrix01, class Matrix02,
         class Matrix10, class Matrix11, class Matrix12,
         class Matrix20, class Matrix21, class Matrix22>
struct isBlockMatrix<MatrixBlock_3x3<Matrix00, Matrix01, Matrix02,
                                     Matrix10, Matrix11, Matrix12,
                                     Matrix20, Matrix21, Matrix22>> { static const bool value = true; };

template<class Vector0,
         class Vector1,
         class Vector2>
struct isBlockVector<VectorBlock_3<Vector0,Vector1,Vector2>> { static const bool value = true; };

//---------------------------------------------------------------------------//
template<class Vector0,
         class Vector1,
         class Vector2,
         class Vector3>
class VectorBlock_4;

template<class Matrix00, class Matrix01, class Matrix02, class Matrix03,
         class Matrix10, class Matrix11, class Matrix12, class Matrix13,
         class Matrix20, class Matrix21, class Matrix22, class Matrix23,
         class Matrix30, class Matrix31, class Matrix32, class Matrix33>
class MatrixBlock_4x4;

template<class Matrix00, class Matrix01, class Matrix02, class Matrix03,
         class Matrix10, class Matrix11, class Matrix12, class Matrix13,
         class Matrix20, class Matrix21, class Matrix22, class Matrix23,
         class Matrix30, class Matrix31, class Matrix32, class Matrix33>
struct isBlockMatrix<MatrixBlock_4x4<Matrix00, Matrix01, Matrix02, Matrix03,
                                     Matrix10, Matrix11, Matrix12, Matrix13,
                                     Matrix20, Matrix21, Matrix22, Matrix23,
                                     Matrix30, Matrix31, Matrix32, Matrix33>> { static const bool value = true; };

template<class Vector0,
         class Vector1,
         class Vector2,
         class Vector3>
struct isBlockVector<VectorBlock_4<Vector0,Vector1,Vector2,Vector3>> { static const bool value = true; };
} // namespace BLA
} // namespace numpack 

#endif //block_TYPE_H
