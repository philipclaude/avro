#ifndef AVRO_LIB_VORONOI_PARTICLES_H_
#define AVRO_LIB_VORONOI_PARTICLES_H_

#include "avro_types.h"

#include "element/simplex.h"

#include "voronoi/new/diagram.h"

#include <memory>
#include <string>
#include <vector>

namespace avro
{

namespace voronoi
{

class PowerDiagram;

class ParticleSimulator {

public:
  ParticleSimulator( const std::string& domain_name , index_t nb_particles );


  void simulate( index_t nb_time_step );
  void sample( const std::string& sprinkle_method , const std::string& smoothing_method );
  void set_density();

  void euler_step();
  void add_gravity_force();
  void add_pressure_force();

  const PowerDiagram& diagram() const { return *diagram_.get(); }
  PowerDiagram& diagram() { return *diagram_.get(); }

  const std::vector<real_t>& density() const { return density_; }
  const std::vector<real_t>& mass() const { return mass_; }

  void save_every( int x ) { save_every_ = x; }
  void save_frames( const std::string& filename ) const;

private:

  // a helper field class to visualize the colour of each voronoi cell
  class FluidField : public Field<Polytope,real_t> {
  public:
    FluidField( ParticleSimulator& particles ) :
      Field(particles.diagram(),0,DISCONTINUOUS)
    {
      build();

      const PowerDiagram& diagram = particles.diagram();

      for (index_t k = 0; k < diagram.nb(); k++) {
        index_t p = diagram.polytope2site()[k];
        this->value(k) = particles.density()[k];
      }
    }

    index_t nb_rank() const { return 1; }
  };


  void initialize( const std::string& domain_name );

  coord_t number_;
  coord_t dimension_;
  index_t nb_particles_;

  std::shared_ptr<Topology<Simplex>> domain_;
  std::shared_ptr<PowerDiagram> diagram_;

  std::vector<real_t> mass_;
  std::vector<real_t> density_;

  std::vector<real_t> position_;
  std::vector<real_t> velocity_;
  std::vector<real_t> force_;

  std::vector<real_t> gravity_;
  real_t tau_;
  real_t eps_;

  std::shared_ptr<FluidField> fluid_field_;

  real_t time_;
  index_t iteration_;

  int save_every_;
  std::vector<float> frame_coordinates_;
  std::vector<float> frame_density_;
};

}

}

#endif
