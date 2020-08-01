#ifndef AVRO_API_AVRO_H_
#define AVRO_API_AVRO_H_

#include <adaptation/parameters.h>

#include <common/types.h>

#include <memory>
#include <vector>

namespace avro
{

static const real_t unset_value = 1e20;

class Entity;
class Model;
class Simplex;
template <typename type> class Topology;
class Points;

class Context
{
public:
  Context( coord_t number , coord_t dim , coord_t udim );

  // this should be called once
  void define_geometry( const std::string& geometry );
  void attach_geometry(); // only at the beginning

  // these should be called anytime an adaptation problem is defined
  void define_mesh( const std::string& mesh );
  void load_mesh( const std::vector<real_t>& x , const std::vector<index_t>& s );
  void load_geometry( const std::vector<int>& g , const std::vector<real_t>& u );

  // the actual adaptation function
  int adapt( const std::vector<real_t>& metric );

  // retrieval of the adapted mesh
  void retrieve_mesh( std::vector<real_t>& x , std::vector<index_t>& s ) const;
  void retrieve_geometry( std::vector<int>& g , std::vector<real_t>& u ) const;
  void retrieve_boundary( std::vector<std::vector<index_t>>& faces ,
                          std::vector<int>& geometry , bool interior=false ) const;

  AdaptationParameters& parameters() { return parameters_; }

private:
  Entity* id2geometry( int id ) const;

  coord_t number_;
  coord_t dim_;
  coord_t udim_;

  std::shared_ptr<Model> model_;
  std::vector<Entity*> entities_;
  AdaptationParameters parameters_;

  std::shared_ptr<Points> points_;
  std::shared_ptr<Topology<Simplex>> topology_;
  std::vector<real_t> metric_;

};

} // avro

#endif
