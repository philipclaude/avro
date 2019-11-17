#include "common/tools.h"

#include "master/transfer.hpp"

#include "mesh/builder.h"

#include <set>
#include <type_traits>
#include <vector>

namespace luna
{

typedef struct
{
  std::vector<index_t> parents;
  std::vector<index_t> local;
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

        printInline( f.indices );

        // check if this facet exists
        std::map<Element,FacetParent>::iterator it = facets_j.find(f);
        if ( it==facets_j.end() )
        {
          std::vector<index_t> parents = {k};
          FacetParent p;
          p.parents.push_back(k);
          p.local.push_back(i);
          facets_j.insert( {f,p} );
        }
        else
        {
          printf("    --> exists! adding parent (%lu,%lu)\n",k,i);
          it->second.parents.push_back(k);
          it->second.local.push_back(i);
        }
      }
    }
  }
}

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
  luna_assert( topology_.nb() == this->nb() );
  luna_assert_msg( f.points().nb()==0 , "nb_vertices = %lu" , f.points().nb() );
  luna_assert( f.points().dim()==topology_.points().dim() );
  luna_assert( f.nb()==0 );

  // create all the vertices for the outgoing topology
  const std::vector<index_t>& elems = this->elements();
  index_t nb_vertices = *std::max_element( elems.begin() , elems.end() ) +1;
  std::vector<real_t> x0( topology_.points().dim() , 0. );
  for (index_t k=0;k<nb_vertices;k++)
    f.points().create( x0.data() );

  // copy the topology
  for (index_t k=0;k<this->nb();k++)
    f.add( this->operator()(k) , this->nv(k) );

  printf("nb vertices = %lu\n",nb_vertices);

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

template<typename Shape_t,typename Master_t>
void
Builder<Shape_t,Master_t>::build()
{
  FacetDecomposition<Shape_t> facets( topology_ );
  facets.build();

  // allocate enough space for all the elements
  std::vector<index_t> elem;
  for (index_t k=0;k<topology_.nb();k++)
  {
    elem.resize( master_.nb_basis() , 0 );
    add( elem.data() , elem.size() );
  }

  printf("nb elements = %lu, nb_basis = %lu, order = %u\n",this->nb(),master_.nb_basis(),master_.order());

  // loop through the dimensional hierarchy
  index_t n = 0;
  for (coord_t dim=0;dim<facets.nb_dim();dim++)
  {
    const std::map<Element,FacetParent>& facets_d = facets[dim];
    std::map<Element,FacetParent>::const_iterator it;

    // loop through all the facets
    for (it=facets_d.begin();it!=facets_d.end();++it)
    {
      const Element& f = it->first;
      luna_assert( f.dim == dim );

      const std::vector<index_t>& parents = it->second.parents;
      const std::vector<index_t>& local = it->second.local;

      // sprinkle the new dof into place
      std::vector<index_t> dof( master_.nb_interior(dim) );

      for (index_t j=0;j<dof.size();j++)
        dof[j] = n++;

      // assign the dof to all parents of this facet
      for (index_t i=0;i<parents.size();i++)
      {
        index_t k = parents[i];

        for (index_t j=0;j<dof.size();j++)
        {
          index_t idx = master_.get_index( f.dim , local[i] , j );
          (*this)(k,idx) = dof[j];
        }
      }
    }
  }
}

// builder for high-order meshes and transfer from lagrange to bezier
template class Builder< Simplex<Lagrange> , Simplex<Lagrange> >;
template class Builder< Simplex<Lagrange> , Simplex<Bezier> >;


} // luna
