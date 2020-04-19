#include "adaptation/adapt.h"
#include "adaptation/metric.h"

#include "mesh/topology.h"

#include <set>

namespace avro
{

template<typename type>
void
AdaptThread<type>::smooth_points( index_t nb_iter )
{
  std::vector<index_t> elems;
  std::vector<index_t> N;

  smoother_.nb_parameter_tests() = 0;
  smoother_.nb_parameter_rejections() = 0;
  smoother_.nb_geometry() = 0;

  real_t Q0 = worst_quality(topology_,metric_);

  // loop over smoothing iterations
  printf("-> performing vertex smoothing:\n");
  for (index_t iter=0;iter<nb_iter;iter++)
  {

    smoother_.delta() = 0.0;
    smoother_.delta_min() = 1e20;
    smoother_.delta_max() = -1.;

    smoother_.resetRejections();
    smoother_.nb_accepted() = 0;

    smoother_.objective() = 0.0;

    // loop through the points
    for (index_t k=0;k<topology_.points().nb();k++)
    {
      if (k<topology_.points().nb_ghost()) continue;
      smoother_.apply( k , metric_ , Q0 );
    }
    smoother_.objective() /= (topology_.points().nb()/topology_.points().dim());

    printf("\titer[%lu]: dx = %3.2e, min = %3.2e, max = %3.2e -> accepted %lu\n",iter,smoother_.delta()/topology_.points().nb(),smoother_.delta_min(),smoother_.delta_max(),
          smoother_.nb_accepted());
    printf("\t\tnb_visibility_rej = %lu, nb_implied_metric_rej = %lu, nb_enlarged_rej = %lu/%lu, Navg = %lu\n",
      smoother_.nb_visibility_rejections(),smoother_.nb_implied_metric_rejections(),
      smoother_.nb_enlarged_rejections(),smoother_.nb_geometry(),index_t(smoother_.Ntot()/topology_.points().nb()));
      printf("\t\tobjective = %1.12e, nb_error = %lu, nb_interp_outside = %lu\n",smoother_.objective(),smoother_.nb_zero_valency(),smoother_.nb_interpolated_outside());
  }
  printf("\tdone %lu iterations of smoothing.\n\tnb_elem = %lu, nb_vert = %lu. parameter rej = (%lu/%lu)\n",
    nb_iter,topology_.nb(),topology_.points().nb(),
    smoother_.nb_parameter_rejections(),smoother_.nb_parameter_tests());
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

  avro_assert( ep->number()==2 );
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
  this->extract_geometry( ep , {p} );
  avro_assert( this->G_.nb()>0 );
  //avro_assert_msg( this->G_.nb_real()==this->nb_ghost() , "|g| = %lu, |c| = %lu\n" , this->G_.nb_real() , this->nb_ghost() );

  for (coord_t d=0;d<2;d++)
    this->u_[this->v2u_[p]][d] = params[d];

  // the following assertion won't hold when p is on an Edge
  if (this->topology_.points().entity(p)->number()==2 && !this->G_.closed())
  {
    avro_assert_msg( this->G_.nb()==this->nb_ghost() ,
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
  avro_assert( checker.hasPositiveVolumes(this->gcavity_,mtype));

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
    if (!this->curved_)
      ep->inverse(x,params);
    else
    {
      ep->inverse_guess(x,params);

      // check the re-evaluated parameter coordinates are close to the new x coordinates
      #if 0
      real_t result[18];
      EGADS::Object* eg_entity = (EGADS::Object*)ep;
      EGADS_ENSURE_SUCCESS( EG_evaluate(*eg_entity->object(),params.data(),result));
      real_t dg = numerics::distance(result,x.data(),dim);
      real_t gtol = 0;
      EGADS_ENSURE_SUCCESS( EG_tolerance(*eg_entity->object(),&gtol) );
      avro_assert_msg( dg < gtol , "x = (%g,%g,%g), x(u) = (%g,%g,%g), distance = %g\n",
                        x[0],x[1],x[2],result[0],result[1],result[2],dg);
      #endif
    }

    // make sure the cavity is not enlarged
    std::vector<index_t> C0 = this->C_;
    this->find_geometry( x.data() , this->C_ );
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
  if (this->nb_removed_nodes()>0)
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

  // apply the parameter space coordinates if necessary
  if (ep!=NULL)
  {
    for (coord_t d=0;d<udim;d++)
      this->topology_.points().u(p,d) = params[d];
  }

  // check if all produce elements have a positive determinant of implied metric
  if (!this->positive_implied_metrics())
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
  index_t elem0 = metric.attachment()[p].elem();

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

template class Smooth<Simplex>;
template class AdaptThread<Simplex>;

} // avro
