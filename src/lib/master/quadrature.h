#ifndef URSA_LIB_MASTER_QUADRATURE_H_
#define URSA_LIB_MASTER_QUADRATURE_H_

#include "common/types.h"

#include <vector>

namespace ursa
{

class Quadrature
{
public:
  Quadrature( coord_t dim , const int order=-1 ); // default to max quad
  virtual void retrieve( std::vector<real_t>& x , std::vector<real_t>& w ) = 0;
  virtual void define( const int order ) = 0;

protected:
  virtual ~Quadrature() {}

  coord_t dim_;
  index_t rule_;
  index_t order_;
  index_t nb_quad_;

  bool defined_;
};

class ConicalProductQuadrature : public Quadrature
{
public:
  using Quadrature::Quadrature;

  void define( const int order );
  void retrieve( std::vector<real_t>& x , std::vector<real_t>& w );
};

class GrundmannMollerQuadrature : public Quadrature
{
public:
  using Quadrature::Quadrature;

  void define( const int order );
  void retrieve( std::vector<real_t>& x , std::vector<real_t>& w );
};

} // ursa

#endif
