#include "mesh/builder.h"

#include <type_traits>

namespace ursa
{

typedef struct
{
  std::vector<index_t> indices;
  std::vector<index_t> parents;
  coord_t dim;
} Facet;

template<typename Master_t>
class FacetDecomposition
{
public:
  FacetDecomposition( const Data<index_t>& topology , const Master_t& master ) :
    topology_(topology),
    master_(master)
  {}

  coord_t nb_dim() const { return facets_.size(); }

  const std::vector<Facet>& operator[](index_t d)
    { return facets_[d]; }

  void build()
  {
    ursa_implement;
  }

private:
  std::vector< std::vector<Facet> > facets_;
  const Data<index_t>& topology_;
  const Master_t& master_;
};

template<typename MasterFrom_t,typename MasterTo_t,typename dof_t>
Builder<MasterFrom_t,MasterTo_t,dof_t>::Builder( const Topology<MasterFrom_t>& topology , const MasterTo_t& masterTo ) :
  topology_(topology),
  masterFrom_(topology.master()),
  masterTo_(masterTo)
{
  // make sure the type is std::vector<real_t>
  ursa_assert( typeid(dof_t)==typeid(std::vector<real_t>) );

  // save the dof we convert from
  for (index_t k=0;k<topology.vertices().nb();k++)
  {
    std::vector<real_t> x( topology.vertices()[k] , topology.vertices()[k]+ topology.vertices().dim() );
    dofFrom_.push_back(x);
  }

  build();
}


template<typename MasterFrom_t,typename MasterTo_t,typename dof_t>
Builder<MasterFrom_t,MasterTo_t,dof_t>::Builder( const Field<MasterFrom_t,dof_t>& field , const MasterTo_t& masterTo ) :
  topology_(field),
  masterFrom_(field.master()),
  masterTo_(masterTo)
{
  const std::vector<dof_t>& data = field.data();

  // save the dof we convert from
  for (index_t k=0;k<data.size();k++)
  {
    dofFrom_.push_back( data[k] );
  }

  build();
}

template<typename MasterFrom_t,typename MasterTo_t,typename dof_t>
void
Builder<MasterFrom_t,MasterTo_t,dof_t>::build()
{

  FacetDecomposition<MasterFrom_t> facets( topology_ , masterFrom_ );
  facets.build();

  // loop through the dimensional hierarchy
  for (coord_t dim=0;dim<facets.nb_dim();dim++)
  {
    const std::vector<Facet> facets_d = facets[dim];

    // loop through all the facets
    for (index_t k=0;k<facets_d.size();k++)
    {
      const Facet& f = facets_d[k];
      ursa_assert( f.dim == dim );

      const std::vector<index_t>& idx = f.indices;

      // get the dof of this facet
      std::vector<dof_t> f_dof( idx.size() );
      for (index_t j=0;j<idx.size();j++)
        f_dof[j] = dofFrom_[ idx[j] ];

      // compute the interior dof of the new facet
      std::vector<dof_t> g_dof;
      masterTo_.template convert<MasterFrom_t,dof_t>( masterFrom_ , f_dof , g_dof );

      // add the dof to the new list
      // TODO check which ones actually need to be added
      for (index_t j=0;j<g_dof.size();j++)
        dofTo_.push_back( g_dof[j] );
    }

  }

}

// builder for high-order meshes
template class Builder< Simplex<Lagrange> , Simplex<Lagrange> , std::vector<real_t> >;


} // ursa
