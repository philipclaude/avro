#ifndef AVRO_SRC_LIB_SIMULATION_H_
#define AVRO_SRC_LIB_SIMULATION_H_

template<typename type>
class PDE
{

public:
  bool has_advective_flux() const;
  bool has_diffusive_flux() const;
};

#endif
