#ifndef LUNA_LIB_ADAPTATION_METRIC_H_
#define LUNA_LIB_ADAPTATION_METRIC_H_

#include "comon/types.h"

#include "mesh/field.h"

namespace luna
{

template<typename type>
class MetricField : public Field< type , MatrixSymD<real_t> >
{

public:
  MetricField( coord_t dim );

private:
  coord_t dim_;
};

} // luna

#endif
