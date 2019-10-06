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

template<typename Shape_t>
class FacetDecomposition
{
public:
  FacetDecomposition( const Topology<Shape_t>& topology ) :
    topology_(topology)
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
};

template<typename Shape_t,typename Master_t>
Builder<Shape_t,Master_t>::Builder( const Topology<Shape_t>& topology , const Master_t& master ) :
  topology_(topology),
  master_(master)
{
  build();
}

template<typename Shape_t,typename Master_t>
void
Builder<Shape_t,Master_t>::transfer( Topology<Master_t>& f ) const
{
  ursa_assert( topology_.nb() == this->nb() );
  ursa_assert( f.vertices().nb()==0 );
  ursa_assert( f.vertices().dim()==topology_.vertices().dim() );
  ursa_assert( f.nb()==0 );

  // create all the vertices for the outgoing topology
  const std::vector<index_t>& elems = this->elements();
  index_t nb_vertices = *std::max_element( elems.begin() , elems.end() );
  std::vector<real_t> x0( topology_.vertices().dim() , 0. );
  for (index_t k=0;k<nb_vertices;k++)
    f.vertices().create( x0.data() );

  // map all the vertices from the topology to f's vertices
  std::vector<bool> visited( nb_vertices , false );
  std::vector<const real_t*> dof0,dof1;
  std::vector<index_t> idx;
  for (index_t k=0;k<topology_.nb();k++)
  {
    // get the vertices of the current element
    dof0.resize( topology_.nv(k) , NULL );
    for (index_t j=0;j<topology_.nv(k);j++)
      dof0[j] = topology_.vertices()[ topology_(k,j) ];

    // size the vertices to be added
    dof1.resize( this->nv(k) , NULL );
    //master_.transfer( topology_.master() , dof0 , dof1 );

    idx.resize( this->nv(k) , 0 );
    for (index_t j=0;j<dof1.size();j++)
    {
      idx[j] = (*this)(k,j);

      // skip visited vertices
      if (visited[ idx[j] ]) continue;
      visited[ idx[j] ] = true;

      for (index_t d=0;d<f.vertices().dim();d++)
        f.vertices()[idx[j]][d] = dof1[j][d];
    }
    // create the element in the topology
    f.add( idx.data() , idx.size() );
  }
}

/*
template<typename Shape_t,typename Master_t>
template<typename MasterFrom_t,typename T>
void
Builder<Shape_t,Master_t>::transfer( const Field<Shape_t,MasterFrom_t,T>& fx , Field<Shape_t,Master_t>& fy ) const
{
  ursa_implement;
}
*/


/*
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
*/

template<typename Shape_t,typename Master_t>
void
Builder<Shape_t,Master_t>::build()
{

  FacetDecomposition<Shape_t> facets( topology_ );
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

      // sprinkle the new dof into place

    }

  }
}

// builder for high-order meshes
template class Builder< Simplex<Lagrange> , Simplex<Lagrange> >;
template class Builder< Simplex<Lagrange> , Simplex<Bezier> >;


} // ursa
