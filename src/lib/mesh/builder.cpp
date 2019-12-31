#include "common/tools.h"

#include "master/transfer.hpp"

#include "mesh/builder.h"

#include <set>
#include <type_traits>
#include <vector>

namespace avro
{

typedef struct
{
  std::vector<index_t> parents;
  std::vector<index_t> local;
  std::vector< std::vector<index_t> > canonical;
} FacetParent;

template<typename Shape_t>
class FacetDecomposition
{
public:
  FacetDecomposition( const Topology<Shape_t>& topology ) :
    topology_(topology)
  {}

  coord_t nb_dim() const { return facets_.size(); }

  const std::map<Element,FacetParent>& operator[](index_t d) const
    { return facets_[d]; }

  std::map<Element,FacetParent>& operator[](index_t d)
    { return facets_[d]; }

  void build();

private:
  std::vector< std::map<Element,FacetParent> > facets_;
  const Topology<Shape_t>& topology_;
};

template<typename Shape_t>
void
FacetDecomposition<Shape_t>::build()
{
  const coord_t nb_facet_dim = topology_.master().number()+1; // always!

  // create the maps for each dimension
  for (index_t j=0;j<nb_facet_dim;j++)
    facets_.push_back( std::map<Element,FacetParent>() );

  // loop through the elements in the topology
  for (index_t k=0;k<topology_.nb();k++)
  {
    // loop through the facets of this element
    for (index_t j=0;j<nb_facet_dim;j++)
    {
      std::map<Element,FacetParent>& facets_j = facets_[j];

      for (index_t i=0;i<topology_.master().nb_facets(j);i++)
      {
        Element f;
        f.dim = j;
        topology_.master().get_facet_vertices( topology_(k) , topology_.nv(k) , i , f );

        std::vector<index_t> canonical(f.indices.size());
        topology_.master().get_canonical_indices( topology_(k) , topology_.nv(k) , f , canonical );

        // check if this facet exists
        std::map<Element,FacetParent>::iterator it = facets_j.find(f);
        if ( it==facets_j.end() )
        {
          std::vector<index_t> parents = {k};
          FacetParent p;
          p.parents.push_back(k);
          p.local.push_back(i);
          p.canonical.push_back( canonical );
          facets_j.insert( {f,p} );
        }
        else
        {
          it->second.parents.push_back(k);
          it->second.local.push_back(i);
          it->second.canonical.push_back(canonical);
        }
      }
    }
  }
}

template<>
Builder<Simplex>::Builder( const Topology<Simplex>& topology , coord_t order , BasisFunctionCategory category ) :
  Table<index_t>(TableLayout_Rectangular,0),
  topology_(topology),
  master_(topology.number(),order)
{
  Table<index_t>::set_rank( master_.nb_basis() );
  master_.set_basis( category );
  build();
}

template<typename type>
void
Builder<type>::transfer( Topology<type>& f ) const
{
  avro_assert( topology_.nb() == this->nb() );
  avro_assert_msg( f.points().nb()==0 , "nb_vertices = %lu" , f.points().nb() );
  avro_assert( f.points().dim()==topology_.points().dim() );
  avro_assert( f.nb()==0 );

  // create all the vertices for the outgoing topology
  const std::vector<index_t>& elems = this->data();
  index_t nb_vertices = *std::max_element( elems.begin() , elems.end() ) +1;
  std::vector<real_t> x0( topology_.points().dim() , 0. );
  for (index_t k=0;k<nb_vertices;k++)
    f.points().create( x0.data() );

  // copy the topology
  for (index_t k=0;k<this->nb();k++)
    f.add( this->operator()(k) , this->nv(k) );

  // map all the vertices from the topology to f's vertices
  std::vector<const real_t*> dof0;
  std::vector<real_t*> dof1;
  for (index_t k=0;k<topology_.nb();k++)
  {
    // get the vertices of the current element
    dof0.resize( topology_.nv(k) , NULL );
    for (index_t j=0;j<topology_.nv(k);j++)
      dof0[j] = topology_.points()[ topology_(k,j) ];

    // size the vertices to be added
    dof1.resize( this->nv(k) , NULL );
    for (index_t j=0;j<dof1.size();j++)
    {
      index_t idx = (*this)(k,j);
      dof1[j] = f.points()[ idx ];
    }
    master_.transfer( topology_.master() , dof0 , dof1 , topology_.points().dim() );
  }
}

template<typename type>
void
Builder<type>::build()
{
  FacetDecomposition<type> facets( topology_ );
  facets.build();

  // allocate enough space for all the elements
  std::vector<index_t> elem;
  for (index_t k=0;k<topology_.nb();k++)
  {
    elem.resize( master_.nb_basis() , 0 );
    add( elem.data() , elem.size() );
  }

  // loop through the dimensional hierarchy
  index_t n = 0;
  for (coord_t dim=0;dim<facets.nb_dim();dim++)
  {
    const std::map<Element,FacetParent>& facets_d = facets[dim];
    std::map<Element,FacetParent>::const_iterator it;

    ReferenceElement<type> reference(dim,master_.order());

    // loop through all the facets
    for (it=facets_d.begin();it!=facets_d.end();++it)
    {
      const Element& f = it->first;
      avro_assert( f.dim == dim );

      const std::vector<index_t>& parents = it->second.parents;
      const std::vector<index_t>& local = it->second.local;
      avro_assert( parents.size() == local.size() );
      avro_assert( parents.size() == it->second.canonical.size() );

      // sprinkle the new dof into place
      std::vector<index_t> dof( master_.nb_interior(dim) );
      for (index_t j=0;j<dof.size();j++)
        dof[j] = n++;

      // assign the dof to all parents of this facet
      for (index_t i=0;i<parents.size();i++)
      {
        index_t k = parents[i];
        const std::vector<index_t>& canonical = it->second.canonical[i];

        for (index_t j=0;j<dof.size();j++)
        {

          // retrieve the barycentric coordinates of this interior point in the facet
          const index_t* lf = reference.get_lattice_coordinate( reference.interior(j) );

          // compute the barycentric coordinates of this point in the parent simplex
          // using the canonical vertices of the facet
          std::vector<index_t> ls( master_.number() , 0 );
          for (index_t ii=0;ii<ls.size();ii++)
          {
            for (index_t jj=0;jj<canonical.size();jj++)
              ls[ii] += lf[jj]*topology_.master().reference().get_lattice_coordinate( canonical[jj] )[ii];
          }

          // determine which index in the element this corresponds to
          int idx = master_.reference().find_index( ls.data() );
          avro_assert( idx>=0 );

          (*this)(k,idx) = dof[j];
        }
      }
    }
  }
}

// builder for high-order meshes and transfer from lagrange to bezier
template class Builder<Simplex>;

} // avro
