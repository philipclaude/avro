#include "common/tools.h"

#include "geometry/entity.h"
#include "geometry/egads/object.h"

#include "mesh/boundary.h"
#include "mesh/facets.h"
#include "mesh/points.h"

#include <set>

namespace avro
{

class FacetTable
{
public:
  FacetTable() {}

  bool contains( std::vector<index_t>& f )
  {
    avro_assert( f.size()>0 );
    std::sort( f.begin() , f.end() );
    std::string s = stringify(f[0]);
    for (index_t j=1;j<f.size();j++)
      s += "|"+stringify(f[j]);
    if (contains(s)) return true;
    facet_.insert(s);
    return false;
  }

  bool contains( std::string& s )
  {
    if (std::find(facet_.begin(),facet_.end(),s)==facet_.end())
      return false;
    return true;
  }

private:
  std::set<std::string> facet_;

};

bool
interior( const Points& points , index_t* v , index_t nv )
{
  avro_assert( nv>0 );

  for (index_t j=0;j<nv;j++)
  {
    if (points.entity(v[j])!=NULL)
      return false;
  }
  return true;
}

namespace BoundaryUtils
{

Entity*
geometryFacet( const Points& points , index_t* v , index_t nv , bool elem )
{
  avro_assert( nv>0 );

  Entity* dummy = (Entity*)v;

  // first check if the facet is interior
  if (interior(points,v,nv)) return NULL;
  else
  {
    if (elem)
      return dummy; // some non-null pointer, not important what it is
  }

  // some part of the facet is on a geometry entity, find the intersection
  Entity* e;

  // get all the entities these points are on
  std::vector<Entity*> entities( nv );
  for (index_t k=0;k<nv;k++)
    entities[k] = points.entity( v[k] );

  if (nv==1)
  {
    if (entities[0]==NULL) return NULL;
    return points.entity( v[0] );
  }
  else if (nv==2)
  {
    if (entities[0]==NULL || entities[1]==NULL) return NULL;
    e = entities[0]->intersect( entities[1] );
    return e;
  }
  else if (nv==3)
  {
    if (entities[0]==NULL || entities[1]==NULL || entities[2]==NULL) return NULL;
    e = entities[0]->intersect( entities[1] , entities[2] );
    return e;
  }
  else if (nv==4)
  {
    if (entities[0]==NULL || entities[1]==NULL ||
        entities[2]==NULL || entities[3]==NULL ) return NULL;
    e = entities[0]->intersect( entities[1] , entities[2] , entities[3] );
    return e;
  }
  avro_implement;
  return NULL;
}

} // BoundaryUtils

Entity*
lookForParent( Entity* e , const Points& points , index_t* v , index_t nv )
{
  avro_assert(e!=NULL);
  for (index_t k=0;k<e->nb_parents();k++)
  {
    if (!e->parents()[k]->tessellatable()) continue;
    index_t count = 0;
    for (index_t j=0;j<nv;j++)
    {
      Entity* ej = points.entity(v[j]);
      if (e->parents()[k]->above(ej))
        count++;
    }
    if (count==nv) return e->parents()[k];
  }
  avro_assert_not_reached;
  return NULL;
}

template<typename type>
void
Boundary<type>::extractall()
{
  typedef std::shared_ptr<Topology<type>> Topology_ptr;
  std::set<Entity*> entities;

  // first count how many geometries we are dealing with
  for (index_t k=0;k<topology_.points().nb();k++)
  {
    // retrieve the entity this vertex is on
    Entity* e = topology_.points().entity(k);

    // interior points are not associated with a boundary
    if (e==NULL) continue;

    // add the entity to the unique set
    entities.insert( e );
  }

  // add all the parents of the entities
  std::vector<Entity*> E;
  std::set<Entity*>::iterator eit;
  for (eit=entities.begin();eit!=entities.end();++eit)
  {
    Entity* e = *eit;
    E.push_back( e );
    for (index_t j=0;j<e->nb_parents();j++)
    {
      Entity* parent = e->parents(j);
      if (!parent->tessellatable())
        continue;
      E.push_back( e->parents(j) );
    }
  }
  uniquify(E);

  // create the sub-topologies to hold the facets
  entity2child_.clear();
  entity_.clear();
  for (index_t k=0;k<E.size();k++)
  {
    Entity* e = E[k];
    Topology_ptr t =
        std::make_shared<Topology<type>>(topology_.points(),e->number());
    entity2child_.insert( std::pair<Entity*,index_t>(e,this->nb_children() ) );
    entity_.push_back(e);
    this->add_child(t);
  }

  // loop through all the simplices
  FacetTable table;
  for (index_t k=0;k<topology_.nb();k++)
  {
    // check if this is entirely volume first
    if (BoundaryUtils::geometryFacet(topology_.points(),topology_(k),topology_.nv(k),true)==NULL)
      continue;

    index_t nf = topology_.master().nb_facets();
    std::vector<index_t> facet;

    // loop through all lower-dimensional facets of this simplex
    for (index_t j=0;j<nf;j++)
    {
      // get the indices of this facet
      // TODO implement this!!!
      avro_implement;
      //topology_.master().unwindFacet( j , topology_(k) , facet );

      // get the geometry this facet is on and lookup which topology this is
      Entity* e = BoundaryUtils::geometryFacet( topology_.points() ,facet.data() , facet.size() );
      if (e==NULL) continue;

      // the entity number must match the facet number
      if (facet.size()-1!=e->number())
        continue;

      // check if all facets share a ghost element
      if (topology_.closed())
      {
        // get all elements in common
        std::vector<index_t> common;
        topology_.all_with( facet , common );
        bool hasghost = false;
        for (index_t i=0;i<common.size();i++)
        {
          if (topology_.ghost(common[i]))
          {
            hasghost = true;
            break;
          }
        }
        if (!hasghost) continue;
      }
      else
      {
        // we can only ensure the facet is counted once in the topology
        // only true if we want dim-1 facets!
        std::vector<index_t> common;
        topology_.all_with( facet , common );
        if (common.size()!=1) continue;
      }

      // this is a geometry facet
      // check if the table contains it
      if (table.contains(facet)) continue;

      // add the facet to the topology
      index_t id = entity2child_[e];
      this->child(id).add( facet.data() , facet.size() );

    }
  }
}

template<typename type>
void
Boundary<type>::extract( bool interior )
{
  typedef std::shared_ptr<Topology<type>> Topology_ptr;

  // first list all the entities attached to the points
  std::set<Entity*> entities;
  const coord_t num = topology_.number();

  // first count how many geometries we are dealing with
  for (index_t k=0;k<topology_.points().nb();k++)
  {
    // retrieve the entity this vertex is on
    Entity* e = topology_.points().entity(k);

    // interior points are not associated with a boundary
    if (e==NULL) continue;

    // add the entity to the unique set
    entities.insert( e );
  }

  // add all the parents of the entities
  std::vector<Entity*> E0;
  std::set<Entity*>::iterator eit;
  for (eit=entities.begin();eit!=entities.end();++eit)
  {
    Entity* e = *eit;
    E0.push_back( e );
    for (index_t j=0;j<e->nb_parents();j++)
    {
      Entity* parent = e->parents(j);
      if (!parent->tessellatable())
        continue;
      E0.push_back( e->parents(j) );
    }
  }
  uniquify(E0);

  // only add the entities of the appropriate topological number
  entity_.clear();
  for (index_t j=0;j<E0.size();j++)
  {
    if (E0[j]->number()==num-1)
        entity_.push_back( E0[j] );
  }

  // create the sub-topologies to hold the facets
  entity2child_.clear();
  for (index_t k=0;k<entity_.size();k++)
  {
    // each child topology t will be associated with the entity e
    Entity* e = entity_[k];
    Topology_ptr t =
        std::make_shared<Topology<type>>(topology_.points(),e->number());
    entity2child_.insert( std::pair<Entity*,index_t>(e,k) );
    this->add_child(t);
  }

  // build up the list of facets on the boundary
  Topology<type> bnd( topology_.points() , num-1 );

  if (topology_.closed() || topology_.points().nb_ghost()>0)
  {
    std::vector<index_t> facet(num);
    for (index_t k=0;k<topology_.nb();k++)
    {

      if (!topology_.ghost(k)) continue;
      avro_assert( topology_(k,0)<topology_.points().nb_ghost() );

      for (index_t j=1;j<topology_.nv(k);j++)
      {
        avro_assert( topology_(k,j)>=topology_.points().nb_ghost() );
        facet[j-1] = topology_(k,j);
      }
      bnd.add( facet.data() , facet.size() );
    }
    if (bnd.nb()==0)
    {
      // the topology was probably not closed
      topology_.get_boundary(bnd);
    }
  }
  else
    topology_.get_boundary(bnd);

  for (index_t k=0;k<bnd.nb();k++)
  {
    // get the geometry this bnd facet is on and lookup which topology this is
    Entity* e = BoundaryUtils::geometryFacet( topology_.points() , bnd(k) , bnd.nv(k) );
    if (e==NULL)
    {
      print_inline( bnd.get(k) , "facet is not on geometry entity: " );
      for (index_t j=0;j<bnd.nv(k);j++)
      {
        Entity* ent = topology_.points().entity( bnd(k,j) );
        printf("vertex %lu:\n",bnd(k,j) );
        topology_.points().print( bnd(k,j) , true );
        if (ent==NULL) continue;
        ent->print();
        printf("parents:\n");
        for (index_t ii=0;ii<ent->nb_parents();ii++)
          ent->parents()[ii]->print();
      }

      //Gamma<type> gamma;
      //library::Plottable<type> plot(bnd);
      //gamma.writeMesh(plot,"bnd_debug.mesh");
    }
    avro_assert(e!=NULL);

    if (e->number()!=bnd.number())
    {
      e = lookForParent(e,topology_.points(),bnd(k),bnd.nv(k));
    }

    // add the facet to the topology
    if (entity2child_.find(e)==entity2child_.end())
    {
      printf("entity not allocated!\n");
      e->print();
      print_inline( bnd.get(k) , "facet is not on geometry entity: " );
      for (index_t j=0;j<bnd.nv(k);j++)
      {
        Entity* ent = topology_.points().entity( bnd(k,j) );
        printf("vertex %lu:\n",bnd(k,j) );
        topology_.points().print( bnd(k,j) , true );
        ent->print();
        printf("parents:\n");
        for (index_t ii=0;ii<ent->nb_parents();ii++)
          ent->parents()[ii]->print();
      }
    }
    index_t id = entity2child_[e];
    this->child(id).add( bnd(k) , bnd.nv(k) );
    if (e!=entity_[id])
      e->print();
    avro_assert( e==entity_[id] );
  }

  if (interior)
  {
    // identify interior facets which might be on the geometry
    Facets facets(topology_);
    facets.compute();
    const index_t nv = topology_.number();
    std::vector<index_t> F(nv);

    // loop through facets and check which ones are on geometry and interior
    for (index_t k=0;k<facets.nb();k++)
    {
      if (topology_.closed())
      {
        // skip actual boundary facets
        if (topology_.ghost( facets.side0(k) ) || topology_.ghost( facets.side1(k) ) )
          continue;
      }
      else
      {
        if (facets.boundary(k))
          continue;
      }

      // retrieve the facet indices
      facets.retrieve(k,F);

      // determine if this is a geometry facet

      Entity* e = NULL;
      try
      {
        e = BoundaryUtils::geometryFacet( topology_.points() , F.data() , F.size() );
      }
      catch(...)
      {
      }

      // the entity is null for interior facets
      // and should have the correct topological number
      if (e==NULL) continue;
      if (e->number()!=topology_.number()-1) continue;
      if (!e->interior()) continue;

      // should be a geometry facet!
      index_t id = entity2child_[e];
      this->child(id).add( F.data() , F.size() );
    }
  }

}

template<typename type>
bool
Boundary<type>::check() const
{
  // check all (n-1) facets are counted once

  Topology<type> boundary( this->points_ , topology_.number()-1 );
  this->get_elements( boundary );

  index_t offenders = 0;
  std::vector<index_t> cells;
  for (index_t k=0;k<boundary.nb();k++)
  {
    topology_.all_with( boundary.get(k) , cells );
    if (cells.size()!=1)
    {

      if (cells.size()==2)
      {
        if (topology_.ghost(cells[0]) || topology_.ghost(cells[1]))
          continue;
      }

      print_inline( boundary.get(k) , "counted " + stringify(cells.size()) +" times:" );
      print_inline( cells , "sharing cells: " );
      for (index_t i=0;i<cells.size();i++)
      {
        if (topology_.ghost(cells[i])) printf("ghost ->");
        else printf("real ->");
        print_inline( topology_.get(cells[i]));
      }
      offenders++;
    }
  }
  if (offenders!=0)
  {
    printf("there are %lu offenders\n",offenders);
    return false;
  }

  return true;
}

template<typename type>
index_t
Boundary<type>::indexof( Entity* e0 )
{
  EGADS::Object* e = (EGADS::Object*)e0;
  std::map<Entity*,index_t>::const_iterator it;
  e->print(false);
  for (it=entity2child_.begin();it!=entity2child_.end();++it)
  {
    //avro_implement; // TODO lookup object?
    EGADS::Object* obj = (EGADS::Object*)it->first;
    obj->print(false);
    if (obj->object() == e->object()) return it->second;
    //if (it->first->object()==e->object()) return it->second;
  }
  printf("could not find entity:\n");
  e->print();

  for (it=entity2child_.begin();it!=entity2child_.end();++it)
    it->first->print();
  avro_assert_not_reached;
  return 0;
}

template<typename type>
void
Boundary<type>::print() const
{
  std::map<Entity*,index_t>::const_iterator it;
  for (it=entity2child_.begin();it!=entity2child_.end();++it)
  {
    it->first->print();
    printf("\tat index %lu with %lu elements\n",it->second,this->child(it->second).nb());
  }
}

template<typename type>
Entity*
Boundary<type>::entity( index_t k )
{
  avro_assert( k < entity_.size() );
  return entity_[k];
}

BoundaryPoints::BoundaryPoints( const Points& points , bool meta ) :
  Points(points.dim())
{
  index_t count = 0;
  for (index_t k=0;k<points.nb();k++)
  {
    // skip interior points
    if (!points.entity(k)) continue;

    create( points[k] );

    if (meta)
    {
      set_entity(count,points.entity(k) );
    }
    local2global_.insert({count,k});
    count++;
  }

  if (count==0)
  {
    printf("[warning] no boundary points found! adding all points!\n");
    for (index_t k=0;k<points.nb();k++)
    {
      if (k<points.nb_ghost()) continue;
      create(points[k]);
      local2global_.insert({count,k});
      count++;
    }
  }
}

template class Boundary<Simplex>;

} // avro
