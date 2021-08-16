
#include "library/factory.h"
#include "mesh/points.h"
#include "mesh/topology.h"

#include "voronoi/diagram.h"
#include "voronoi/particles.h"

#include <json/json.hpp>

#include <fstream>

namespace avro
{

namespace voronoi
{

ParticleSimulator::ParticleSimulator( const std::string& domain_name , index_t nb_particles ) :
  nb_particles_(nb_particles),
  domain_(nullptr),
  diagram_(nullptr),
  save_every_(-1)
{
  initialize(domain_name);
}

void
ParticleSimulator::initialize( const std::string& domain_name ) {

  // read in the mesh
  std::shared_ptr<TopologyBase> ptopology = nullptr;
  std::shared_ptr<Mesh> pmesh = library::get_mesh(domain_name,ptopology);
  domain_ = std::static_pointer_cast<Topology<Simplex>>(ptopology);

  number_ = domain_->number();
  dimension_ = domain_->points().dim() + 1; // + 1 since we will optimize weights which requires codimension 1

  diagram_ = std::make_shared<PowerDiagram>(*domain_.get(),dimension_);
  diagram_->set_ambient_dimension(domain_->points().dim());

  // TODO dimension or dimension+1 ?
  position_.resize( nb_particles_ * (dimension_-1) , 0.0 );
  velocity_.resize( nb_particles_ * (dimension_-1),  0.0 ); // initial velocity is zero
  mass_.resize( nb_particles_ , 0.0 );
  density_.resize( nb_particles_ , 1.0 );
}

void
ParticleSimulator::sample( const std::string& sprinkle_method , const std::string& smoothing_method ) {

  if (sprinkle_method == "random") {

    // randomly sprinkle points within the domain
    // assume CKF for now...
    for (index_t k = 0; k < nb_particles_; k++) {
      for (coord_t d = 0; d < dimension_-1; d++) {
        position_[k*(dimension_-1)+d] = random_within(0.0,1.0);
      }
    }
  }
  else
    avro_implement;

  diagram_->allocate_sites(nb_particles_);
  diagram_->set_sites(position_.data());
  diagram_->initialize();

  if (smoothing_method == "lloyd")
    diagram_->optimize_points_lloyd(150);
  else
    avro_implement;

  // save the point coordinates
  for (index_t k = 0; k < nb_particles_; k++)
  for (coord_t d = 0; d < dimension_-1; d++)
    position_[k*(dimension_-1)+d] = diagram_->sites()[k][d];
}

void
ParticleSimulator::set_density() {

  // retrieve the initial volume of each cell
  const std::vector<real_t>& volumes = diagram_->cell_volume();

  // compute the initial mass of each cell
  for (index_t k = 0; k < nb_particles_; k++) {

    real_t density = 1.0;

    real_t x,y,z = 0.0;
    x = diagram_->sites()[k][0];
    y = diagram_->sites()[k][1];
    //if (dimension_ > 2) z = position_[k*(dimension_)+2];

    real_t f = 0.1*sin(3*M_PI*x);
    if (y-0.5 > f)
      density = 10.0;

    density_[k] = density;
    mass_[k]    = density * volumes[k];
  }
}

void
print_header() {
  printf("--------------------------------------------------------------------\n");
  printf("%6s|%9s|%9s|%9s|%9s|%4s|%9s\n","step  "," time ","dt" , " dx " , "resid." , "iter" , "min_vol" );
  printf("--------------------------------------------------------------------\n");
}

void
ParticleSimulator::euler_step() {

  // option to save the position at this time step
  if ( save_every_ > 0 && iteration_ % save_every_ == 0) {
    index_t n = frame_coordinates_.size();
    frame_coordinates_.resize( n + position_.size() );
    for (index_t i = 0; i < position_.size(); i++)
      frame_coordinates_[n+i] = position_[i];
  }

  // solve for the weights that achieve the target mass
  std::vector<real_t> volumes( nb_particles_ );
  for (index_t k = 0; k < nb_particles_; k++)
    volumes[k] = mass_[k] / density_[k];

  bool converged = diagram_->optimize_weights_kmt( 5 , volumes );

  real_t vmin = * std::min_element( diagram_->cell_volume().begin() , diagram_->cell_volume().end() );
  if (vmin < 1e-16) converged = false;

  if (!converged) {
    printf("[error] optimal transport solver did not converge (try decreasing time step)\n");
    printf("        minimum volume = %g\n",vmin);
    avro_implement;
  }

  // TODO start with a time step but make it smaller if the weight solver does not converge
  tau_ = 2e-4;
  eps_ = 4e-2;

  time_ += tau_;

  // retrieve the particle centroids
  const std::vector<real_t>& centroids = diagram_->centroids();

  coord_t dim = dimension_ -1;

  // compute the forces on particles
  real_t dx = 0.0;
  for (index_t k = 0; k < nb_particles_; k++) {

    real_t fx, fy;

    fx = ( centroids[k*dim  ] - position_[k*dim  ] ) / (eps_*eps_);
    fy = ( centroids[k*dim+1] - position_[k*dim+1] ) / (eps_*eps_);

    fy -= mass_[k] * 9.81;

    velocity_[k*dim  ] += tau_ * fx / mass_[k];
    velocity_[k*dim+1] += tau_ * fy / mass_[k];

    position_[k*dim  ] += tau_ * velocity_[k*dim  ];
    position_[k*dim+1] += tau_ * velocity_[k*dim+1];

    for (coord_t d = 0; d < dim; d++) {
      if (position_[k*dim+d] < 0.0) {
        position_[k*dim+d] = 0.0;
        velocity_[k*dim+d] *= -1.0;
      }
      if (position_[k*dim+d] > 1.0) {
        position_[k*dim+d] = 1.0;
        velocity_[k*dim+d] *= -1.0;
      }
    }

    dx += velocity_[k*dim ]*velocity_[k*dim] + velocity_[k*dim+1]*velocity_[k*dim+1];

  }
  dx = std::sqrt(tau_*tau_*dx);

  if (iteration_ % 100 == 0) print_header();
  printf("%6lu|%6.3e|%6.3e|%6.3e|%6.3e|%4lu|%6.3e\n",iteration_,time_,tau_,dx,diagram_->residual(),diagram_->iteration(),vmin);

  for (index_t k = 0; k < nb_particles_; k++)
  for (coord_t d = 0; d < dim; d++)
    diagram_->sites()[k][d] = position_[k*dim+d];
}


void
ParticleSimulator::simulate( index_t nb_time_step ) {

  diagram_->set_verbose(false);
  time_ = 0.0;
  for (index_t k = 0; k < nb_time_step; k++) {
    iteration_ = k;
    euler_step();
  }

  diagram_->accumulate();

  // overwrite the site field with the density data
  fluid_field_ = std::make_shared<FluidField>(*this);
  diagram_->fields().make("fluid",fluid_field_);

  //diagram_->create_field();
}

void
ParticleSimulator::save_frames( const std::string& filename ) const {

  // use float instead of double to represent coordinates
  using json_flt = nlohmann::basic_json<std::map, std::vector, std::string, bool,
                                       std::int64_t, std::uint64_t, float>;

  json_flt jp;

  jp["nb_frames"] = frame_coordinates_.size() / position_.size();
  jp["nb_particles"] = nb_particles_;
  jp["nb_properties"] = 0;

  std::vector<float> data(frame_coordinates_.begin(),frame_coordinates_.end());
  jp["coordinates"] = data;
  jp["density"]    = density_;
  jp["properties"] = frame_density_;
  jp["dimension"] = dimension_-1;

  std::ofstream jfile(filename);
  jfile << jp;
  jfile.close();
}

} // voronoi

} // avro
