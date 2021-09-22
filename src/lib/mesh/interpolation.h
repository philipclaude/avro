//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_MESH_FIELD_INTERPOLATION_H_
#define avro_LIB_MESH_FIELD_INTERPOLATION_H_

#include "avro_types.h"

#include <memory>
#include <vector>

namespace avro
{

class Points;
template<typename type,typename T> class Field;
template<typename type> class ElementSearch;

template<typename type,typename T>
class FieldInterpolation
{

public:
  FieldInterpolation( const Field<type,T>* field=nullptr );
  virtual ~FieldInterpolation() {}

  virtual int eval( const Points& points , index_t p , const std::vector<index_t>& guesses , T& tp );

  bool analytic() const { return analytic_; }
  const ElementSearch<type>& searcher() const {
    return *searcher_.get();
  }

  const Field<type,T>& field() const { return *pfield_; }

protected:
  bool analytic_;

private:
  const Field<type,T>* pfield_;
  std::shared_ptr<ElementSearch<type>> searcher_;
};

} // avro

#endif
