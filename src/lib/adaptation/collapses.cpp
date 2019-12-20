#include "adaptation/adapt.h"
#include "adaptation/metric.h"
#include "adaptation/primitive.h"

#include "geometry/entity.h"
#include "geometry/egads/object.h"

namespace luna
{

template<typename type>
Collapse<type>::Collapse( Topology<type>& _topology ) :
  Primitive<type>(_topology)
{
  this->setName("collapser");
  nb_accepted_.resize(_topology.number()+1);
  nb_rejected_.resize(_topology.number()+1);
  //this->curved_ = false;
}

template<typename type>
bool
Collapse<type>::visibleParameterSpace( index_t p , index_t q , Entity* g0 , bool edge )
{
  EGADS::Object* g = (EGADS::Object*) g0;
  luna_assert( g->number()==2 );
  if (!this->curved_) return true;
  if (g->interior()) return true;

  // extract the geometry cavity
  this->extractGeometry( g, {p} );
  luna_assert( this->G_.nb()>0 );

  // based on the type of face, we may need to flip the sign for the volume calculation
  int oclass,mtype;
  ego ref,prev,next;
  EGADS_ENSURE_SUCCESS( EG_getInfo(*g->object(), &oclass, &mtype,&ref, &prev, &next) );
  if (mtype==SREVERSE)
    this->gcavity_.sign() = -1;
  else
    this->gcavity_.sign() = 1;

  // compute the cavity about the original configuration
  GeometryOrientationChecker checker( this->points_ , this->u_ , this->u2v_ , g );
  bool accept = this->gcavity_.compute( this->v2u_[p] , this->u_[this->v2u_[p]] , this->S_ );
  luna_assert( accept );

  // ensure the normals are originally in the right direction
  int s = checker.signof( this->gcavity_ );
  luna_assert( s > 0 );
  luna_assert( checker.hasPositiveVolumes(this->gcavity_,mtype));

  // check visibility in the new configuration
  nb_parameter_tests_++;
  accept = this->gcavity_.compute( this->v2u_[q] , this->u_[this->v2u_[q]] , this->S_ );
  if (!accept)
  {
    if (edge) nb_rej_vis_Edge_++;
    nb_parameter_rejections_++;
    return false;
  }

  // check the normal orientations
  s = checker.signof( this->gcavity_ );
  if (s<0)
  {
    if (edge) nb_rej_sgn_Edge_++;
    return false;
  }

  // check the geometry topology can be extracted correctly
  if (checker.createsBadGeometry(this->gcavity_) || checker.createsBadGeometry(this->boundary()))
  {
    if (edge) nb_rej_geo_Edge_++;
    nb_invalid_geometry_++;
    return false;
  }

  // ensure we have all positive volumes
  luna_assert( checker.hasPositiveVolumes(this->gcavity_,mtype));

  return true;
}

template<typename type>
bool
Collapse<type>::valid( const index_t p , const index_t q )
{
  if (this->topology_.points().fixed(p)) return false;

  // check if a lower dimensional geometry entity is collapsed to a higher
  // dimensional one
  Entity* e0 = this->topology_.points().entity(p);
  Entity* e1 = this->topology_.points().entity(q);

  // check if the removed vertex is a volume one (always valid)
  if (e0==NULL) return true;

  if (e1==NULL)
  {
    // e0 is has a non-null geometry, so e1 cannot be a volume point
    return false;
  }

  // e0 must be higher in the geometry hierarchy than e1
  if (e0->above(e1))
  {
    if (this->geometry(p,q)==NULL) return false;
    return true;
  }

  // if the entities are not the same then we cannot collapse the edge
  // this also accounts for node-node collapses
  if (e0!=e1) return false;

  // check that the edge has a ghost
  if (this->geometry(p,q)==NULL) return false;

  if (e0!=NULL && e1!=NULL)
  {
    // check if the edge between the points contains a ghost
    bool contains = false;
    std::vector<index_t> shell;

    this->topology_.intersect( {p,q} , shell );
    for (index_t k=0;k<shell.size();k++)
    {
      if (this->topology_.ghost(shell[k]))
      {
        contains = true;
        break;
      }
    }
    Entity* g = this->geometry(p,q);
    if (!contains && !g->interior()) return false;
  }

  return true;
}

template<typename type>
bool
Collapse<type>::apply( const index_t p , const index_t q , bool delay )
{

  // attempt to collapse p onto q
  // assign the cavity
  this->topology_.intersect( {p} , this->C_ );

  this->info_ = "trying to collapse " + stringify<index_t>(p) +
                " onto " + stringify<index_t>(q);

  // count initial number of ghosts
  index_t nb_ghost0 = 0;
  for (index_t k=0;k<this->C_.size();k++)
  {
    if (this->topology_.ghost(this->C_[k]))
      nb_ghost0++;
  }

  // turn off enlarging
  this->enlarge_ = false;
  bool accept = this->compute( q , this->topology_.points()[q] , this->C_ );
  if (!accept)
  {
    // the point is not visible
    return false;
  }

  if (this->invalidatesTopology())
  {
    // the collapse invalidates the topology in the sense that
    // closed bodies (with number n) do not have n+1 points (i.e. disappear)
    return false;
  }

  if (!this->philipcondition())
  {
    // we don't want duplicate elements! can happen with ghosted topologies..
    return false;
  }

  // check if all produce elements have a positive determinant of implied metric
  if (!this->positiveImpliedMetrics())
    return false;

  // determine if the receiving point is visible in the parameter space
  Entity* g = this->geometry(p,q);
  if (g!=NULL && g->number()==2 && !visibleParameterSpace(p,q,g))
  {
    nb_rejected_[g->number()]++;
    return false;
  }

  // determine if the point is visible on every face parenting the Edge
  if (g!=NULL && this->topology_.number()==3 && g->number()==1 && this->curved_)
  {
    // we need to check for visibility on all Faces that parent this Edge
    for (index_t k=0;k<g->nb_parents();k++)
    {
      Entity* parent = g->parents(k);
      if (parent->number()!=2 || !parent->tessellatable())
      {
        continue;
      }
      if (!visibleParameterSpace(p,q,parent,true)) // flag this is an edge for the counters
      {
        nb_rejected_[g->number()]++;
        return false;
      }
    }
  }

  // check if only ghosts are created
  index_t only_ghost = true;
  for (index_t k=0;k<this->nb();k++)
  {
    if (!this->ghost(k))
    {
      only_ghost = false;
      break;
    }
  }
  if (only_ghost) return false;

  // apply the operator if no delay was requested
  if (!delay)
    this->topology_.apply(*this);

  if (g==NULL) nb_accepted_[this->topology_.number()]++;
  else nb_accepted_[g->number()]++;
  return true;
}

template class Collapse<Simplex>;

} // luna
