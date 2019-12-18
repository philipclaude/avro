#include "common/tools.h"

#include "geometry/entity.h"
#include "geometry/egads/object.h"

#include "library/metric.h"

#include "adaptation/metric.h"
#include "adaptation/primitive.h"

#include <egads.h>

#include <set>

namespace luna
{

#define CROSS(a,b,c)      a[0] = (b[1]*c[2]) - (b[2]*c[1]);\
                          a[1] = (b[2]*c[0]) - (b[0]*c[2]);\
                          a[2] = (b[0]*c[1]) - (b[1]*c[0])
#define DOT(a,b)         (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])

#define FLIP_MTYPE 1

class GeometryOrientationChecker
{
public:
  GeometryOrientationChecker( Points& v , Points& u,
                              const std::vector<index_t>& _u2v,
                              Entity* entity ) :
    v_(v),
    u_(u),
    u2v_(_u2v),
    entity_((EGADS::Object*)entity),
    normal_(TableLayout_Rectangular,3)
  {
    luna_assert( u_.dim()==v_.dim()-1 );
    luna_assert( v_.dim()==3 );
    luna_assert( u2v_.size()==u_.nb() );

    // compute the normal to every vertex
    std::vector<real_t> N(v_.dim());
    for (index_t k=0;k<u_.nb();k++)
    {
      if (k<u_.nb_ghost())
      {
        std::fill(N.begin(),N.end(),0);
        normal_.add(N.data(),3);
        continue;
      }
      normal( k , N );
      normal_.add( N.data() , 3 );
    }
  }

  template<typename type> int signof( Topology<type>& topology , bool verbose=false )
  {
    luna_assert( topology.number()==u_.dim() );

    std::vector<int> sign( topology.nb() );

    // loop through all the cells (facets) in the topology
    std::vector<real_t> N(v_.dim());
    int s = 1;
    index_t counted = 0;
    for (index_t k=0;k<topology.nb();k++)
    {
      if (topology.ghost(k)) continue;
      counted++;

      const real_t* v0 = v_[ u2v_.at(topology(k,0)) ];
      const real_t* v1 = v_[ u2v_.at(topology(k,1)) ];
      const real_t* v2 = v_[ u2v_.at(topology(k,2)) ];

      if (verbose)
      {
        printf("checking facet (%lu,%lu,%lu):\n",
                u2v_.at(topology(k,0)),u2v_.at(topology(k,1)),u2v_.at(topology(k,2)));
        printf("v0 = (%g,%g,%g)\n",v0[0],v0[1],v0[2]);
        printf("v1 = (%g,%g,%g)\n",v1[0],v1[1],v1[2]);
        printf("v2 = (%g,%g,%g)\n",v2[0],v2[1],v2[2]);
      }

      // compute the normal of this facet in physical space
      normal( v0 , v1 , v2 , N );

      sign[k] = 1;
      for (index_t j=0;j<3;j++)
      {
        luna_assert( topology(k,j)>=topology.points().nb_ghost() );
        const real_t* n = normal_(topology(k,j));
        real_t dp = DOT(N,n);
        if (verbose)
        {
          printf("vertex %lu dot product = %g\n",topology(k,j),dp);
          std::vector<real_t> n0(n,n+3);
          print_inline(N,"nt");
          print_inline(n0,"ng");
        }
        if (dp<0.1)
        {
          sign[k] = -1;
          s = -1;
        }
      }
    }
    luna_assert( counted>0 );
    return s;
  }

  template<typename type> bool hasPositiveVolumes( const Topology<type>& topology , int sign=1 ) const
  {
    std::vector<real_t> vg(topology.nb());
    topology.get_volumes(vg);
    for (index_t k=0;k<vg.size();k++)
    {
      if (topology.ghost(k)) continue;
      if (sign*vg[k]<=0.0)
      {
        print_inline(vg,"vols");
        return false;
      }
    }
    return true;
  }

  template<typename type> bool createsBadGeometry( const Topology<type>& topology ) const
  {
    for (index_t k=0;k<topology.nb();k++)
    {
      if (topology.ghost(k)) continue;

      // check for weird geometry configurations that cannot be recovered from
      Entity* e0 = topology.points().entity( topology(k,0) );
      Entity* e1 = topology.points().entity( topology(k,1) );
      Entity* e2 = topology.points().entity( topology(k,2) );

      if (e0==NULL || e1==NULL || e2==NULL) continue;

      Entity* face = e0->intersect(e1,e2,true); // only check, don't error!
      if (face==NULL) return true;
    }
    return false;
  }

  void normal( index_t k , std::vector<real_t>& N )
  {
    luna_assert( N.size()==3 );

    // evaluate the normal at the vertex using the stored parameter coordinates
    real_t result[18];
    EGADS_ENSURE_SUCCESS( EG_evaluate( *entity_->object() , u_[k] , result ) );
    real_t dx_du[3],dx_dv[3];
    dx_du[0] = result[3];
    dx_du[1] = result[4];
    dx_du[2] = result[5];
    dx_dv[0] = result[6];
    dx_dv[1] = result[7];
    dx_dv[2] = result[8];
    CROSS(N, dx_du, dx_dv);
    numerics::normalize( N.data() , N.size() );

    // flip the normal depending on the type of face (SFORWARD or SREVERSE)
    int oclass,mtype;
    ego ref,prev,next;
    EGADS_ENSURE_SUCCESS( EG_getInfo(*entity_->object(), &oclass, &mtype,&ref, &prev, &next) );
    N[0] *= mtype;
    N[1] *= mtype;
    N[2] *= mtype;
  }

  void normal( const real_t* v0 , const real_t* v1 , const real_t* v2 ,
               std::vector<real_t>& N )
  {
    // evaluate the normal of the oriented triangle
    luna_assert( N.size()==3 );
    real_t X01[3],X02[3];
    for (coord_t d=0;d<3;d++)
    {
      X01[d] = v1[d] -v0[d];
      X02[d] = v2[d] -v0[d];
    }
    CROSS(N,X01,X02);
    real_t ln = std::sqrt(N[0]*N[0]+N[1]*N[1]+N[2]*N[2]);
    luna_assert(ln!=0.0);
    numerics::normalize( N.data() , N.size() );
  }

private:
  Points& v_; // original points in physical space
  Points& u_; // points in parameter space
  const std::vector<index_t>& u2v_; // map from parameter to physical coordinates
  EGADS::Object* entity_; // the geometry entity in question

  Table<real_t> normal_; // vertex normals (using the geometry)
};

template<typename type>
Entity*
Primitive<type>::geometry( index_t p0 , index_t p1 )
{
  Entity* e0 = this->topology_.points().entity(p0);
  Entity* e1 = this->topology_.points().entity(p1);
  if (e0==NULL || e1==NULL) return NULL;

  Entity* g = e0->intersect(e1);
  if (g==NULL) return NULL;

  if (g->interior()) return g; // skip the ghost check

  // we need to make sure the edge is attached to some ghosts
  std::vector<index_t> shell;
  this->topology_.intersect( {p0,p1} , shell );

  for (index_t k=0;k<shell.size();k++)
  {
    if (this->topology_.ghost(shell[k]))
      return g;
  }
  // there are no ghosts, cannot be a geometry edge
  return NULL;
}

template<typename type>
void
Primitive<type>::extractGeometry( Entity* e , const std::vector<index_t>& f )
{
  luna_assert( e->number()==2 );
  u_.clear();
  G_.clear();
  v2u_.clear();
  u2v_.clear();
  gcavity_.clear();
  S_.clear();
  G_.neighbours().forceCompute(); // forces the neighbours to be recomputed
  G_.set_closed(false); // forces the re-closing of the mesh

  this->computeGeometry( e , G_ , v2u_ , u2v_ );
  if (G_.nb()==0) return;

  G_.inverse().build();
  if (f.size()==0)
  {
    // a size of zero is a code for extracting non-ghost elements
    for (index_t k=0;k<G_.nb();k++)
    {
      if (G_.ghost(k)) continue;
      S_.push_back(k);
    }
  }
  else if (f.size()==1)
  {
    // requested a vertex v, look up the u value
    index_t u = v2u_.at(f[0]);
    G_.inverse().ball(u,S_);
  }
  else if (f.size()==2)
  {
    // requested an edge, lookup non-ghost elements
    index_t u0 = v2u_.at(f[0]);
    index_t u1 = v2u_.at(f[1]);
    G_.inverse().shell(u0,u1,S_);
    luna_assert( S_.size()==2 );
  }
  else
  {
    print_inline(f,"unsupported facet: ");
    luna_assert_not_reached;
  }
}

template<typename type>
Collapse<type>::Collapse( Topology<type>& _topology ) :
  Primitive<type>(_topology)
{
  this->setName("collapser");
  nb_accepted_.resize(_topology.number()+1);
  nb_rejected_.resize(_topology.number()+1);
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
        continue;
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

template<typename type>
Insert<type>::Insert( Topology<type>& _topology ) :
  Primitive<type>(_topology)
{
  this->setName("inserter");
}

template<typename type>
bool
Insert<type>::visibleParameterSpace( real_t* x , real_t* params , Entity* ep0 )
{
  EGADS::Object* ep = (EGADS::Object*) ep0;

  luna_assert( ep->number()==2 );
  if (!this->curved_) return true;
  if (ep->interior()) return true;

  // based on the type of face, we may need to flip the sign for the volume calculation
  int oclass,mtype;
  ego ref,prev,next;
  EGADS_ENSURE_SUCCESS( EG_getInfo(*ep->object(), &oclass, &mtype,&ref, &prev, &next) );
  if (mtype==SREVERSE)
    this->gcavity_.sign() = -1;
  else
    this->gcavity_.sign() = 1;

  // extract the geometry cavity
  this->extractGeometry( ep ); // no facet information, will assign the cavity to non-ghosts
  luna_assert(this->G_.nb()>0);

  // add the parameter coordinates for the inserted point along with the mapped index
  this->u_.create( params );
  this->u2v_.push_back( this->points_.nb()-1 );

  // the insertion should be visible in parameter space
  bool accept = this->gcavity_.compute( this->u_.nb()-1 , params , this->S_ );
  luna_assert(accept);

  // check the orientation of the facets
  GeometryOrientationChecker checker(this->points_,this->u_,this->u2v_,ep);
  int s = checker.signof( this->gcavity_ );
  if (s<0) return false;

  luna_assert( checker.hasPositiveVolumes(this->gcavity_,mtype) );

  return true;
}


template<typename type>
bool
Insert<type>::apply( const index_t e0 , const index_t e1 , real_t* x , real_t* u , const std::vector<index_t>& shell )
{
  //std::vector<index_t> elems;
  elems_.clear();

  // insertions should not enlarge the initial set of cavity elements
  this->enlarge_ = false;
  this->checkVisibility_ = true;
  enlarged_ = false; // no enlargement (for geometry) yet

  // determine the cavity elements if none were provided
  if (shell.size()==0)
  {
    // likely an interior edge, cavity is the shell around the edge
    this->topology_.intersect( {e0,e1} , elems_ );
  }
  else
  {
    // likely a geometry insertion, use the provided shell
    elems_ = shell;
  }

  // index of the vertex to be inserted
  index_t ns = this->topology_.points().nb();

  // we need to add the vertex
  this->topology_.points().create(x);

  // this might be a boundary insertion,
  // compute intersection of geometry entities of each vertex
  Entity* entity0 = this->topology_.points().entity(e0);
  Entity* entity1 = this->topology_.points().entity(e1);
  Entity* entitys;
  int bodys;
  if (entity0==NULL || entity1==NULL)
  {
    // interior split
    entitys = NULL;
    bodys   = 0;
  }
  else
  {
    // set the geometry of the inserted vertex as the intersection
    // of the endpoint geometries
    entitys = this->geometry(e0,e1);

    // determine the body
    if (entitys==NULL)
    {
      // interior split
      bodys = 0;
    }
    else
    {
      int body0 = this->topology_.points().body(e0);
      int body1 = this->topology_.points().body(e1);

      if (body0==body1) bodys = body0;
      else bodys = 0; // not sure how this can even happen but just being safe
    }
  }

  if (entitys!=NULL)
  {
    // enlarge the cavity for boundary insertions
    // if the shell was given, then this should have been precomputed
    std::vector<index_t> elems0 = elems_;
    if (!this->findGeometry( x , elems_ ))
    {
      // we could not enlarge to find the inserted point
      this->topology_.remove_point(ns);
      return false;
    }
    if (elems_.size()!=elems0.size())
    {
      enlarged_ = true;
    }
  }

  // attempt the operator
  bool accept = this->compute( ns , x , elems_ );
  if (!accept)
  {
    // the cavity requested enlargment which is possible because of a
    // minimum volume constraint in the cavity operator
    // remove the vertex and inform the caller the insertion is not allowed
    this->topology_.remove_point(ns);
    return false;
  }

  // check if all produce elements have a positive determinant of implied metric
  #if 0
  // pcaplan REMOVE THIS in master branch (this ruins timing and should only be used for curved=true)
  if (!this->positiveImpliedMetrics())
  {
    this->topology_.remove_point(ns);
    return false;
  }
  #endif


  if (entitys!=NULL && entitys->number()==2)
  {
    // check if the insertion is visible to the boundary of the geometry
    // cavity in the parameter space of the geometry entity
    accept = visibleParameterSpace( x , u , entitys );
    if (!accept)
    {
      this->topology_.remove_point(ns);
      return false;
    }
  }

  if (entitys!=NULL && this->topology_.number()==3 && !entitys->interior())
  {
    index_t nb_ghost = 0;
    for (index_t k=0;k<this->nb();k++)
    {
      if (this->ghost(k))
        nb_ghost++;
    }
    PRIMITIVE_CHECK( nb_ghost==4 );
    if (nb_ghost!=4)
    {
      this->topology_.remove_point(ns);
      return false;
    }
  }

  if (this->nb_removedNodes()>0)
  {
    // when we are enlarging, do not allow points to be removed
    this->topology_.remove_point(ns);
    return false;
  }

  // set the geometry entity and body
  this->topology_.points().set_entity( ns , entitys );
  this->topology_.points().set_param( ns , u );
  this->topology_.points().body(ns) = bodys;

  // apply the operator to the topology if the caller did not request a delay
  if (!this->delay_)
    this->topology_.apply(*this);
  return true;
}

template<typename type>
EdgeSwap<type>::EdgeSwap( Topology<type>& _topology ) :
  Primitive<type>(_topology)
{
  this->setName("edge swapper");
  nb_geometry_rejections_.resize( _topology.number() );
}

template<typename type>
bool
EdgeSwap<type>::visibleParameterSpace( index_t p , index_t e0 , index_t e1 , Entity* face0 )
{
  EGADS::Object* face = (EGADS::Object*) face0;

  // check if the edge swap is valid in the parameter space
  // 2d swaps area already checks, 4d is limited to linear geometries for now
  if (this->topology_.number()!=3) return true;
  if (!this->curved_) return true;
  if (face->interior()) return true;

  // based on the type of face, we may need to flip the sign for the volume calculation
  int oclass,mtype;
  ego ref,prev,next;
  EGADS_ENSURE_SUCCESS( EG_getInfo(*face->object(), &oclass, &mtype,&ref, &prev, &next) );
  if (mtype==SREVERSE)
    this->gcavity_.sign() = -1;
  else
    this->gcavity_.sign() = 1;

  // extract the geometry cavity
  this->extractGeometry( face, {e0,e1} );
  luna_assert( this->G_.nb()>0 );

  if (!this->G_.closed())
    luna_assert_msg( this->G_.nb()==2 ,
                    "|G| = %lu, |swap ghosts| = %lu" , this->G_.nb() , this->nb_ghost() );

  if (this->G_.closed())
    luna_assert( this->G_.nb_real()==2 );

  // ensure the e0 is visible to the cavity boundary
  bool accept = this->gcavity_.compute( this->v2u_[e0] , this->u_[this->v2u_[e0]] , this->S_ );
  luna_assert( accept );

  GeometryOrientationChecker checker( this->points_ , this->u_ , this->u2v_ , face );
  int s = checker.signof( this->gcavity_ );
  luna_assert_msg( s > 0 , "negative orientation for edge (%lu,%lu) with vertex %lu" , e0,e1,e0 );
  luna_assert( checker.hasPositiveVolumes(this->gcavity_,mtype));

  // ensure e1 is visible to the cavity boundary
  accept = this->gcavity_.compute( this->v2u_[e1] , this->u_[this->v2u_[e1]] , this->S_ );
  luna_assert( accept );
  s = checker.signof( this->gcavity_ );
  luna_assert_msg( s > 0 , "negative orientation for edge (%lu,%lu) with vertex %lu" , e0,e1,e1 );
  luna_assert( checker.hasPositiveVolumes(this->gcavity_,mtype));

  // check for visibility of p
  nb_parameter_tests_++;
  accept = this->gcavity_.compute( this->v2u_[p] , this->u_[this->v2u_[p]] , this->S_ );
  luna_assert( this->gcavity_.nb_real()==2 );
  if (!accept)
  {
    //printf("swap not visible in parameter space!\n");
    nb_parameter_rejections_++;
    return false;
  }

  s = checker.signof( this->gcavity_ );
  if (s<0)
  {
    //printf("swap creates negative normals!\n");
    return false;
  }

  if (checker.createsBadGeometry(this->gcavity_))
  {
    //printf("swap creates bad geometry!\n");
    nb_invalid_geometry_++;
    return false;
  }

  // the geometry cavity should have positive volumes
  luna_assert( checker.hasPositiveVolumes(this->gcavity_,mtype));

  return true;
}

template<typename type>
bool
EdgeSwap<type>::valid( const index_t p , const index_t e0 , const index_t e1 )
{
  // topology checks
  if (this->topology_.points().fixed(e0) ||
      this->topology_.points().fixed(e1))
    return false;
  if (p<this->topology_.points().nb_ghost()) return false;
  if (e0<this->topology_.points().nb_ghost()) return false;
  if (e1<this->topology_.points().nb_ghost()) return false;

  std::vector<index_t> edge = {e0,e1};
  if (this->C_.empty())
    this->topology_.intersect(edge,this->C_);

  std::vector<index_t>& elems = this->C_;

  for (index_t k=0;k<elems.size();k++)
  for (index_t j=0;j<this->topology_.nv(elems[k]);j++)
  {
    if (this->topology_.points().fixed(this->topology_(elems[k],j)))
      return false;
  }

  // check if the two points lie on an Edge entity
  Entity* ge = this->geometry(e0,e1);

  // check if this edge is in the volume
  if (ge==NULL)
  {
    // edge is in the volume, always valid
    return true;
  }

  // check if this edge is on a geometry Edge (can't swap these)
  if (ge!=NULL && ge->number()<2)
  {
    nb_geometry_rejections_[ge->number()]++;
    return false;
  }

  // check the entity on the re-inserted vertex
  if (ge!=NULL)
  {
    Entity* ep = this->topology_.points().entity(p);
    if (ep==NULL) return false;
    if (ep!=ge)
    {
      if (ge->above(ep))
      {} // the swap is okay
      else
      {
        nb_geometry_rejections_[ge->number()]++;
        return false;
      }
    }
  }

  if (ge->interior())
  {
    nb_wake_++;
  }

  return true;
}

template<typename type>
bool
EdgeSwap<type>::apply( const index_t p , const index_t e0 , const index_t e1 )
{
  // check if the swap is valid in terms of geometry
  if (!valid(p,e0,e1)) return false; // computes cavity C_

  this->info_ = "trying to swap edge (" + stringify<index_t>(e0) + "/"
                + stringify<index_t>(e1) +")"+
                " with reinsertion " + stringify<index_t>(p);

  // compute the original number of ghost elements
  index_t nb_ghost0 = 0;
  for (index_t j=0;j<this->C_.size();j++)
    if (this->topology_.ghost(this->C_[j]))
      nb_ghost0++;

  // apply the operator, checking visibility
  this->enlarge_ = false;
  bool accept = this->compute( p , this->topology_.points()[p] , this->C_ );
  if (!accept) return false;

  // check if all produce elements have a positive determinant of implied metric
  if (!this->positiveImpliedMetrics())
    return false;

  // check visibility in the parametric space
  Entity* ge = this->geometry(e0,e1);
  if (ge!=NULL && ge->number()==2)
  {
    if (!visibleParameterSpace(p,e0,e1,ge))
    {
      //printf("swap not visible in parameter space!\n");
      return false;
    }
  }

  // count the resulting number of ghost elements
  index_t nb_ghost = 0;
  for (index_t j=0;j<this->nb();j++)
    if (this->ghost(j))
      nb_ghost++;

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

  // do not allow swaps which change the number of ghosts for 3-simplices
  if (this->topology_.number()==3 && nb_ghost0!=nb_ghost)
    return false;

  if (!accept) return false;

  return accept;
}

template<typename type>
RidgeSwap<type>::RidgeSwap( Topology<type>& _topology ) :
  Primitive<type>(_topology)
{
  this->setName("ridge swapper");
  nb_geometry_rejections_.resize( _topology.number() );
}

template<typename type>
bool
RidgeSwap<type>::valid( index_t p , const index_t t0 , const index_t t1 , const index_t t2 )
{
  Entity* ep = this->topology_.points().entity(p);
  Entity* e0 = this->topology_.points().entity(t0);
  Entity* e1 = this->topology_.points().entity(t1);
  Entity* e2 = this->topology_.points().entity(t2);

  // for now disregard any ridge swap touch a boundary
  if (ep!=NULL || e0!=NULL || e1!=NULL || e2!=NULL) return false;
  return true;
}

template<typename type>
bool
RidgeSwap<type>::apply( index_t p , const index_t t0 , const index_t t1 , const index_t t2 )
{
  luna_assert( this->topology_.number()==4 ); // should only be used in 4d though technically a face swap in 3d
  if (!valid(p,t0,t1,t2)) return false;

  // get all elements touching this triangle
  this->topology_.intersect( {t0,t1,t2} , this->C_ );

  this->enlarge_ = false;
  bool accept = this->compute( p , this->topology_.points()[p] , this->C_ );

  // check if any points are deleted
  if (this->nb_removedNodes()>0) return false;

  return accept;
}

template<typename type>
FacetSwap<type>::FacetSwap( Topology<type>& _topology ) :
  Primitive<type>(_topology)
{
  this->setName("facet swapper");
}

template<typename type>
bool
FacetSwap<type>::valid( const index_t p , const index_t k0 , const index_t k1 )
{
  if (p<=this->topology_.points().nb_ghost()) return false;
  if (this->topology_.ghost(k0) || this->topology_.ghost(k1))
    return false;
  return true;
}

template<typename type>
bool
FacetSwap<type>::apply( const index_t p , const index_t k1 , const index_t k2 )
{
  // check if the facet swap is valid
  if (!valid(p,k1,k2)) return false;
  this->enlarge_ = false;

  // check if the swap is visible
  std::vector<index_t> elems = {k1,k2};
  bool accept = this->compute( p , this->topology_.points()[p] , elems );

  // check if all produced elements have a positive determinant of implied metric
  if (!this->positiveImpliedMetrics())
    return false;
  return accept;
}

template<typename type>
Smooth<type>::Smooth( Topology<type>& _topology ) :
  Primitive<type >(_topology),
  delta_(0.0),
  delta_min_(1e20),
  delta_max_(-1),
  M0_(_topology.number())
{
  this->setName("smoother");
  resetRejections();
  nb_accepted_ = 0;
}

template<typename type>
bool
Smooth<type>::visibleParameterSpace( index_t p , real_t* x , real_t* params , Entity* ep0 )
{
  EGADS::Object* ep = (EGADS::Object*) ep0;

  luna_assert( ep->number()==2 );
  if (!this->curved_) return true;
  if (ep->interior()) return true;

  // based on the type of face, we may need to flip the sign for the volume calculation
  int oclass,mtype;
  ego ref,prev,next;
  EGADS_ENSURE_SUCCESS( EG_getInfo(*ep->object(), &oclass, &mtype,&ref, &prev, &next) );
  if (mtype==SREVERSE)
    this->gcavity_.sign() = -1;
  else
    this->gcavity_.sign() = 1;

  // extract the geometry cavity
  this->extractGeometry( ep , {p} );
  luna_assert( this->G_.nb()>0 );
  //luna_assert_msg( this->G_.nb_real()==this->nb_ghost() , "|g| = %lu, |c| = %lu\n" , this->G_.nb_real() , this->nb_ghost() );

  for (coord_t d=0;d<2;d++)
    this->u_[this->v2u_[p]][d] = params[d];

  // the following assertion won't hold when p is on an Edge
  if (this->topology_.points().entity(p)->number()==2 && !this->G_.closed())
  {
    luna_assert_msg( this->G_.nb()==this->nb_ghost() ,
                      "|G| = %lu, |smooth ghosts| = %lu" , this->G_.nb() , this->nb_ghost() );
  }

  GeometryOrientationChecker checker( this->points_ , this->u_ , this->u2v_ , ep );

  // check for visibility of the new coordinates
  nb_parameter_tests_++;
  bool accept = this->gcavity_.compute( this->v2u_[p] , params , this->S_ );
  if (!accept)
  {
    nb_parameter_rejections_++;
    return false;
  }

  // check the orientation of facets produced w.r.t. the geometry
  int s = checker.signof( this->gcavity_ );
  if (s<0) return false;

  // check the produced volumes are positive
  luna_assert( checker.hasPositiveVolumes(this->gcavity_,mtype));

  return true;
}

template<typename type>
bool
Smooth<type>::apply( const index_t p , MetricField<type>& metric , real_t Q0 )
{
  // compute the cavity around p
  this->C_.clear();
  this->topology_.intersect( {p} , this->C_ );

  const coord_t dim = this->topology_.points().dim();
  std::vector<real_t> F( dim,0. );
  std::vector<real_t> u( dim,0. );
  std::vector<real_t> x( dim,0. );
  std::vector<real_t> x0 (dim,0. );

  const coord_t udim = this->topology_.points().udim();
  std::vector<real_t> params0( udim,0. );
  std::vector<real_t> params( udim,0. );

  // get all the nodes in the cavity (excluding this one)
  std::set<index_t> Ns;
  index_t q;
  Entity* ep = this->topology_.points().entity(p);
  if (ep!=NULL && ep->number()==0) return false; // don't move points on Nodes!
  Entity* eq;
  for (index_t k=0;k<this->C_.size();k++)
  {
    for (index_t j=0;j<this->topology_.nv(this->C_[k]);j++)
    {
      q = this->topology_(this->C_[k],j);
      if (q==p) continue;
      if (q<this->topology_.points().nb_ghost()) continue;
      eq = this->topology_.points().entity(q);
      if (ep!=NULL)
      {
        // for geometry points, the smoothing can only be affected by
        // geometry points which are equal to or higher in the hierarchy
        // than the geometry of the smoothed vertex
        // todo: use "geometryEdge"
        if (eq==NULL) continue;
        if (ep!=eq)
        {
          // the entities are different
          if (!ep->above(eq))
          {
            // ep must at least be higher in the topological hierarchy
            continue; // skip vertex q
          }
        }
      }
      Ns.insert(q);
    }
  }

  if (Ns.size()==0)
  {
    // no points were found to do the smoothing
    // this really should not happen but just in case
    nb_zero_valency_++;
    return false;
  }

  // copy the set of points into the vector of points N
  std::vector<index_t> N;
  for (std::set<index_t>::iterator it=Ns.begin();it!=Ns.end();++it)
    N.push_back(*it);
  Ntot_ += N.size();

  // compute the tensile force
  std::vector<real_t> lens0;
  real_t len;
  real_t f;
  for (index_t k=0;k<N.size();k++)
  {
    // compute the metric length
    len = metric.length(this->topology_.points(),p,N[k]);
    lens0.push_back(len);

    // physical vector from the vertex to the neighbour
    numerics::vector( this->topology_.points()[p] ,
                        this->topology_.points()[N[k]] , dim , u.data() );

    // compute the force on the vertex
    f = std::pow(len,4);
    f = (1. -len)*std::exp(-len);

    for (coord_t d=0;d<dim;d++)
      F[d] += -f*u[d];
  }

  for (coord_t d=0;d<dim;d++)
    objective_ += std::pow( F[d] , 2. );

  // compute the new position of the vertex
  real_t omega = 0.2; // relaxation factor
  for (index_t d=0;d<dim;d++)
  {
    x0[d] = this->topology_.points()[p][d];
    x[d]  = this->topology_.points()[p][d] +omega*F[d];
  }
  real_t delta_p = numerics::distance( x0.data() , x.data() , dim );
  delta_ += delta_p;

  if (delta_p<delta_min_) delta_min_ = delta_p;
  if (delta_p>delta_max_) delta_max_ = delta_p;

  this->enlarge_ = true;

  // check if the cavity needs to be enlarged
  if (ep!=NULL)
  {
    for (index_t d=0;d<udim;d++)
      params0[d] = params[d] = this->topology_.points().u(p,d);

    // use the previous parameter value as the initial guess of the projection
    /*if (!this->curved_)
      ep->project(x);
    else
      ep->projectGuess(x,params);
    */
    if (!this->curved_)
      ep->inverse(x,params);
    else
    {
      ep->inverse_guess(x,params);

      // check the re-evaluated parameter coordinates are close to the new x coordinates
      #if 0
      real_t result[18];
      EGADS_ENSURE_SUCCESS( EG_evaluate(*ep->object(),params.data(),result));
      real_t dg = numerics::distance(result,x.data(),dim);
      real_t gtol = 0;
      EGADS_ENSURE_SUCCESS( EG_tolerance(*ep->object(),&gtol) );
      luna_assert_msg( dg < gtol , "x = (%g,%g,%g), x(u) = (%g,%g,%g), distance = %g\n",
                        x[0],x[1],x[2],result[0],result[1],result[2],dg);
      #endif
    }

    // make sure the cavity is not enlarged
    std::vector<index_t> C0 = this->C_;
    this->findGeometry( x.data() , this->C_ );
    if (this->C_.size()!=C0.size())
    {
      nb_enlarged_rejections_++;
      return false;
    }

    // there are more efficient ways of doing this but this is okay for now...
    std::sort( this->C_.begin() , this->C_.end() );
    std::sort( C0.begin() , C0.end() );
    for (index_t j=0;j<C0.size();j++)
    {
      if (this->C_[j]!=C0[j])
      {
        nb_enlarged_rejections_++;
        return false;
      }
    }

    nb_geometry_++;
  }

  this->enlarge_ = false;
  bool accept = this->compute( p , x.data() , this->C_ );
  if (this->nb_removedNodes()>0)
  {
    nb_removed_rejections_++;
    return false;
  }
  if (!accept)
  {
    // the point was not visible within the cavity
    nb_visibility_rejections_++;
    return false;
  }

  // apply the physical coordinates
  for (coord_t d=0;d<dim;d++)
    this->topology_.points()[p][d] = x[d];

  // apply the parameter space coordinates if necesary
  if (ep!=NULL)
  {
    for (coord_t d=0;d<udim;d++)
      this->topology_.points().u(p,d) = params[d];
  }

  // check if all produce elements have a positive determinant of implied metric
  if (!this->positiveImpliedMetrics())
  {
    nb_implied_metric_rejections_++;

    // revert the physical coordinates
    for (index_t d=0;d<dim;d++)
      this->topology_.points()[p][d] = x0[d];

    // revert the parameter space coordinates
    for (coord_t d=0;d<udim;d++)
      this->topology_.points().u(p,d) = params0[d];

    // signal the smoothing swas not applied
    return false;
  }

  // check if the point is visible on the cavity boundary
  if (ep!=NULL && ep->number()==2 && this->curved_)
  {
    if (!visibleParameterSpace( p , x.data(), params.data() , ep ))
    {
      // revert the physical coordinates
      for (index_t d=0;d<dim;d++)
        this->topology_.points()[p][d] = x0[d];

      // revert the parameter space coordinates
      for (coord_t d=0;d<udim;d++)
        this->topology_.points().u(p,d) = params0[d];

      // signal the smoothing was not applied
      return false;
    }
  }

  if (ep!=NULL && this->topology_.number()==3 && ep->number()==1 && this->curved_)
  {
    // we need to check for visibility on all Faces that parent this Edge
    for (index_t k=0;k<ep->nb_parents();k++)
    {
      Entity* parent = ep->parents(k);
      if (parent->number()!=2 || !parent->tessellatable())
        continue;
      if (!visibleParameterSpace(p, x.data() , params.data() ,parent))
      {
        // revert the physical coordinates
        for (index_t d=0;d<dim;d++)
          this->topology_.points()[p][d] = x0[d];

        // revert the parameter space coordinates
        for (coord_t d=0;d<udim;d++)
          this->topology_.points().u(p,d) = params0[d];

        // signal the smoothing was not applied
        return false;
      }
    }
  }

  this->enlarge_ = false;

  // save the previous metric
  for (index_t i=0;i<dim;i++)
  for (index_t j=i;j<dim;j++)
    M0_(i,j) = metric(this->topology_.points(),p)(i,j);
  index_t elem0 = metric.value(p).elem();

  // recompute the metric at the new point
  bool success = metric.recompute( p , x.data() );

  // count this metric interpolation as "outside"
  if (!success)
    nb_interpolated_outside_++;

  if (!success)
  {
    for (index_t d=0;d<dim;d++)
      this->topology_.points()[p][d] = x0[d];

    if (ep!=NULL)
      for (coord_t d=0;d<udim;d++)
        this->topology_.points().u(p,d) = params0[d];

    metric.attachment().assign( p , M0_ , elem0 );
    return false;
  }

  // evaluate the new quality
  if (Q0>0.0)
  {
    if (worst_quality(*this,metric)<Q0)
    {
      // revert the coordinates and re-assign the metric
      for (index_t d=0;d<dim;d++)
        this->topology_.points()[p][d] = x0[d];

      if (ep!=NULL)
        for (coord_t d=0;d<udim;d++)
          this->topology_.points().u(p,d) = params0[d];

      metric.attachment().assign( p , M0_ , elem0 );
    }
  }
  nb_accepted_++;

  return true;
}

template class Collapse<Simplex>;
template class Insert<Simplex>;
template class EdgeSwap<Simplex>;
template class FacetSwap<Simplex>;
template class RidgeSwap<Simplex>;
template class Smooth<Simplex>;
template class Primitive<Simplex>;

} // luna
