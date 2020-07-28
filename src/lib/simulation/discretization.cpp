#include "element/simplex.h"

#include "mesh/field.h"

#include "simulation/discretization.h"
#include "simulation/jacobian.h"

namespace avro
{

class NavierStokes
{
public:
  typedef std::vector<real_t> field_t;

};

template<typename discretization_t,typename pde_t>
void
Assembler<discretization_t,pde_t>::assemble()
{

  const std::string& name = field_.element().name();

  if (name == "simplex")
  {
    const Field<Simplex,field_t>& fld = field_cast<Simplex>();
    discretization_.assemble(fld);
  }
  else if (name == "cube")
  {
    avro_implement;
  }
  else if (name == "polytope")
  {
    avro_implement;
  }
  else
    avro_assert_not_reached;
}

template class Assembler<CG,NavierStokes>;

template class Discretization<CG>;
template class Discretization<FV>;
template class Discretization<DG>;
template class Discretization<FD>;

} // avro
