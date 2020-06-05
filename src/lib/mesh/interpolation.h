#ifndef avro_LIB_MESH_FIELD_INTERPOLATION_H_
#define avro_LIB_MESH_FIELD_INTERPOLATION_H_

#include "common/types.h"

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

protected:
  bool analytic_;

private:
  const Field<type,T>* pfield_;
  std::shared_ptr<ElementSearch<type>> searcher_;
};

} // avro

#endif
