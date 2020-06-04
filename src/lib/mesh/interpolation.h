#ifndef avro_LIB_MESH_FIELD_INTERPOLATION_H_
#define avro_LIB_MESH_FIELD_INTERPOLATION_H_

#include "mesh/search.h"

namespace avro
{

class Points;
template<typename type,typename T> class Field;

template<typename type,typename T>
class FieldInterpolation
{

public:
  FieldInterpolation( const Field<type,T>& field );
  virtual ~FieldInterpolation() {}

  virtual int eval( const Points& points , index_t p , const std::vector<index_t>& guesses , T& tp );

private:
  const Field<type,T>& field_;
  ElementSearch<type> searcher_;
};

template<typename type>
class GeometryMetric : public FieldInterpolation<type,real_t>
{

public:
  GeometryMetric( const Field<type,real_t>& fld );

  int eval( const Points& points , index_t p , const std::vector<index_t>& guesses , real_t& mp ) override;

private:

};

} // avro

#endif
