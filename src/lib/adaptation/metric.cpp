#include "adaptation/cavity.h"
#include "adaptation/implied_metric.h"
#include "adaptation/metric.h"

#include "geometry/entity.h"

#include "library/metric.h"

#include "mesh/field.hpp"
#include "mesh/search.h"
#include "mesh/topology.h"
#include "mesh/points.h"

#include "numerics/determinant.h"
#include "numerics/functions.h" // factorial
#include "numerics/geometry.h"

#include <json/json.hpp>
#include <libmeshb/libmeshb7.h>

#include <cmath>

namespace avro
{

template<typename type>
MetricField<type>::MetricField( Topology<type>& topology , MetricAttachment& fld , FieldInterpolation<type,Metric>* interpolation ) :
	Field<type,Metric>(topology,1,CONTINUOUS),
	topology_(topology),
	number_(topology_.number()),
	attachment_(fld),
	searcher_(topology_),
	evaluator_(nullptr)
{
	avro_assert( attachment_.nb()==topology_.points().nb() );
	avro_assert( attachment_.nb()>0 );

	if (number_==2)
	{
		normalization_ = 4.*std::sqrt(3.);
	}
	else if (number_==3)
	{
		normalization_ = 36./std::pow(3.,1./3.);
	}
	else if (number_==4)
	{
		normalization_ = 10./std::pow( topology_.master().reference().vunit() ,
																	 2./topology_.number() );
	}
	else
		avro_assert_not_reached;

	// build the field
	Field<type,Metric>::build();

  // save the tensors stored in the field
  // these are static and need to be kept for interpolation
  for (index_t k=0;k<attachment_.nb();k++)
	{
		Field<type,Metric>::value(k).allocate(number_);
    Field<type,Metric>::value(k).set( attachment_[k] );
		Field<type,Metric>::value(k).calculate();
	}

	this->master().set_basis( BasisFunctionCategory_Lagrange );
	if (interpolation==nullptr)
		interpolation_ = std::make_shared< FieldInterpolation<type,Metric> >(*this);
	else
		interpolation_.reset(interpolation);


}

template<typename type>
void
MetricField<type>::reset( MetricAttachment& attachment )
{
  attachment_.clear();
  attachment_.reset( attachment );
}

template<typename type>
numerics::SymMatrixD<real_t>&
MetricField<type>::operator() ( const Points& points , index_t p )
{
  avro_assert_msg( p<attachment_.nb() ,
                   "p = %lu but field has %lu tensors",p,attachment_.nb() );
  return attachment_[p];
}

real_t
geometric_interpolation( const numerics::SymMatrixD<real_t>& m0,
	                       const numerics::SymMatrixD<real_t>& m1,
											   const numerics::VectorD<real_t>& edge )
{
	real_t l0_sqr = numpack::Transpose(edge)*m0*edge;
	real_t l1_sqr = numpack::Transpose(edge)*m1*edge;
	real_t l0 = std::sqrt( l0_sqr );
	real_t l1 = std::sqrt( l1_sqr );
	real_t lm,r;
	if (l0<1e-16) return 0.;
	if (l1<1e-16) return 0.;
	if (l0>l1)
	{
		r = l0/l1;
		lm = l0;
	}
	else
	{
		r  = l1/l0;
		lm = l1;
	}
	if (fabs(r-1.)<1e-12)
		return l0;
	if (r!=r)
	{
		m0.dump();
		m1.dump();
	}
	if (r!=r)
	{
		printf("l0_sqr = %g, l1_sqr = %g\n",l0_sqr,l1_sqr);
		printf("l0 = %g, l1 = %g\n",l0,l1);
		edge.dump();
	}
	avro_assert_msg(r>0., "r = %.20e",r);
	return lm*(r -1.)/( r*std::log(r) );
}

template<typename type>
real_t
MetricField<type>::length( index_t n0 , index_t n1 ) const
{
  avro_assert_msg( n0 < attachment_.nb() ,
                  "n0 = %lu, attachment_.nb() = %lu" , n0, attachment_.nb() );
  avro_assert_msg( n1 < attachment_.nb() ,
                  "n1 = %lu, attachment_.nb() = %lu" , n1, attachment_.nb() );

	// get the vector associated with this edge
	Entity* entity = nullptr;
	if (topology_.master().parameter())
	{
		index_t e[2] = {n0,n1};
		entity = BoundaryUtils::geometryFacet( attachment_.points() , e , 2 );
		if (entity->number()==1)
		{
			// find a parent
			for (index_t k=0;k<entity->nb_parents();k++)
			{
				if (entity->parents(k)->number()==2 and entity->parents(k)->tessellatable())
				{
					entity = entity->parents(k);
					break;
				}
			}
		}
		avro_assert( entity->number()!=1 );
	}
	coord_t dim = topology_.points().dim();
	if (topology_.master().parameter()) dim = topology_.points().udim();
	std::vector<real_t> edge0( dim );
	topology_.master().edge_vector( attachment_.points() , n0 , n1 , edge0.data() , entity );
/*
	numerics::vector( attachment_.points()[n0] , attachment_.points()[n1] , topology_.points().dim() , edge0.data() );
*/
	numerics::VectorD<real_t> edge(dim,edge0.data());
  return geometric_interpolation( attachment_[n0] , attachment_[n1] , edge );
}

template<typename type>
void
MetricField<type>::lengths( const Topology<type>& topology , std::vector<real_t>& lens ) const
{
	std::vector<index_t> edges;
	topology.get_edges(edges);
	lens.resize( edges.size()/2 );
	for (index_t k=0;k<edges.size()/2;k++)
		lens[k] = length( edges[2*k] , edges[2*k+1] );
}

template<typename type>
real_t
MetricField<type>::volume( const Topology<type>& topology , const index_t k )
{
	const index_t *V = topology(k);
	const index_t NV = topology.nv(k);
	const type& master = topology.master();

	if (topology.ghost(k)) return 0.0;

	// find the metric with the maximum determinant
	index_t jmax = 0;
	real_t d = attachment_[V[jmax]].sqdet();
	real_t dmax = d;
	for (coord_t j=1;j<NV;j++)
	{
		d = attachment_[V[j]].sqdet();
		if (d>dmax)
		{
			dmax = d;
			jmax = j;
		}
	}

	real_t sqrtdetM = dmax;
	real_t v = sqrtdetM*master.volume(topology.points(),V,NV);
	return v;
}

template<typename type>
real_t
MetricField<type>::volume( const Topology<type>& t )
{
  real_t v = 0.;
  for (index_t k=0;k<t.nb();k++)
  {
    if (t.ghost(k)) continue;
    v += volume(t,k);
  }
  return v;
}

template<typename type>
real_t
MetricField<type>::quality( const Topology<type>& topology , index_t k )
{
	const index_t *V = topology(k);
	const index_t NV = topology.nv(k);
	const Points& points = topology.points();
	const coord_t dim = points.dim();
	const coord_t num = topology.number();
	const type& master = topology.master();

	if (topology.ghost(k)) return -1.;

	// find the metric with the maximum determinant
	index_t jmax = 0, jmin = 0;
	real_t d = attachment_[V[jmax]].sqdet();
	real_t dmax = d;
	real_t dmin = d;
	for (coord_t j=1;j<NV;j++)
	{
		d = attachment_[V[j]].sqdet();
		if (d>dmax)
		{
			dmax = d;
			jmax = j;
		}
		if (d<dmin)
		{
			dmin = d;
			jmin = j;
		}
	}

	real_t sqrtdetM = dmax;
	index_t idxM = jmax;
	UNUSED(jmin);

	// tensor with maximum determinant
	const numerics::SymMatrixD<real_t>& M = attachment_[ V[idxM] ];

	Entity* entity = nullptr;
	if (topology_.master().parameter())
	{
		entity = BoundaryUtils::geometryFacet( attachment_.points() , V , NV );
		avro_assert( entity!=nullptr );
	}

	// compute the edge lengths under m
  real_t l = 0.,lj;
	numerics::VectorD<real_t> e(dim);
  for (index_t j=0;j<master.nb_edges();j++)
  {
    index_t p0 = master.edge(j,0);
    index_t p1 = master.edge(j,1);

		#if 1
		topology_.master().edge_vector( attachment_.points() , V[p0] , V[p1] , e.data() , entity );
		#else
		numerics::vector( points[V[p0]] , points[V[p1]] , dim , e.data() );
		#endif
		lj = numpack::Transpose(e)*M*e;
    l  += lj;
  }

	// compute the volume under m
	real_t v = sqrtdetM*master.volume(topology.points(),V,NV);
	if (v<0)
	{
		for (index_t j=0;j<NV;j++)
			topology.points().print(V[j],true);
	}
	avro_assert_msg( v>=0. , "v = %g, sqrtDetM = %g, v_e = %g" , v , sqrtdetM , master.volume(topology.points(),V,NV) );
	v = std::pow( v , 2./num );

	// normalize to be within [0,1]
	real_t q = normalization_*v/l;
	return q;
}

template<typename type>
bool
MetricField<type>::check( Topology<type>& topology )
{
  if (topology.points().nb()!=attachment_.nb()) return false;

	bool success = true;
	for (index_t k=0;k<attachment_.nb();k++)
	{
		if (topology_.ghost(attachment_[k].elem()))
		{
			printf("element %lu is a ghost containing point %lu\n",attachment_[k].elem(),k);
			success = false;
			break;
		}
	}
  return success;
}

template<typename type>
int
MetricField<type>::find( index_t n0 , index_t n1 , real_t* x )
{
  avro_assert_msg( n0>=topology_.points().nb_ghost() ,
									"n0 = %lu, nb_ghost = %lu" ,
									 n0 , topology_.points().nb_ghost() );
  avro_assert_msg( n1>=topology_.points().nb_ghost() ,
									"n1 = %lu, nb_ghost = %lu" ,
									 n1 , topology_.points().nb_ghost() );

  avro_assert( n0 < attachment_.points().nb() );
  avro_assert( n1 < attachment_.points().nb() );

  avro_assert_msg( n0 < attachment_.nb() ,
									 "n0 = %lu , attachment_.nb = %lu\n" ,
									  n0 , attachment_.nb() );

	// search by starting with the element containing n0
  index_t guess = attachment_[n0].elem();
  avro_assert_msg( guess < topology_.nb() ,
									"guess = %lu but nb = %lu\n", guess,topology_.nb());
  int elem = searcher_.find( x , guess );
  if (elem>=0) return elem;

	// search by starting with the element containing n1, if not found yet
  guess = attachment_[n1].elem();
  avro_assert_msg( guess < topology_.nb() ,
									"guess = %lu but nb = %lu\n", guess,topology_.nb());
  elem = searcher_.find( x , guess );
  if (elem>=0) return elem;

	// brute force the search (slow)
  elem = searcher_.brute(x);
  return elem;
}

template<typename type>
bool
MetricField<type>::add( index_t n0 , index_t n1 , real_t* x )
{

	Metric mp(number_);
	index_t ns = attachment_.points().nb()-1; // hack because edge split
	index_t g0 = attachment_[n0].elem();
	index_t g1 = attachment_[n1].elem();
	int ielem1 = interpolation_->eval( attachment_.points() , ns , {g0,g1} , mp );
	if (ielem1<0) return false;
	attachment_.add( mp , index_t(ielem1) );
	return true;


#if 1  // ultra hack for surface meshing
	attachment_.add( attachment_[n0] , attachment_[n0].elem() );
	avro_assert( attachment_.check() );
	return true;
#endif

  // find the element containing x bordering n0 and n1
  int ielem = find(n0,n1,x);

  // check if the containing element was not found (point outside domain)
  // this can only happen with curved geometries
  if (ielem<0)
  {
		numerics::SymMatrixD<real_t> tensor(number_);
		#if 0
		std::vector<real_t> alpha(2,0.5);
		std::vector<numerics::SymMatrixD<real_t>> metrics(alpha.size());
		metrics[0] = attachment_[n0];
		metrics[1] = attachment_[n1];
		bool interp_result = interp(alpha,metrics,tensor);
		if (!interp_result) return false;
    attachment_.add( tensor , attachment_[n0].elem() );
		#else
		std::vector<real_t> alpha(number_+1,0.);
		ielem = searcher_.closest( x , alpha );
		index_t elem = index_t(ielem);
		std::vector<numerics::SymMatrixD<real_t>> metrics(alpha.size());
		for (index_t j=0;j<alpha.size();j++)
		{
			//metrics[j] = Field<type,Metric>::operator()(elem,j);
			metrics[j] = Field<type,Metric>::value( topology_(elem,j) );
		}
		bool interp_result = interp(alpha,metrics,tensor);
		if (!interp_result) return false;
		attachment_.add( tensor , elem );
		#endif
		return true;
  }

  index_t elem = index_t(ielem);
  const coord_t nv = topology_.nv(elem);

  // get the points of the element
  std::vector<const real_t*> xk(nv,0);
  for (index_t j=0;j<nv;j++)
    xk[j] = topology_.points()[ topology_(elem,j) ];
  real_t V = 1./numerics::simplex_volume(xk,topology_.points().dim());

  // compute the barycentric coordinates
  std::vector<real_t> alpha( nv , 0. );
	std::vector<numerics::SymMatrixD<real_t>> metrics(nv);
  for (index_t j=0;j<nv;j++)
  {
    std::vector<const real_t*> xk0(nv);
    for (index_t i=0;i<nv;i++)
    {
      if (i==j) xk0[i] = x;
      else xk0[i] = xk[i];
    }
    alpha[j] = V*numerics::simplex_volume(xk0,topology_.points().dim());

    if (alpha[j]<0. || alpha[j]>1.)
    {
      printf("element search failed\n");
      avro_assert_not_reached;
    }

    // save the metric at this vertex
		metrics[j] = Field<type,Metric>::operator()(elem,j);
  }

  // perform the interpolation
	numerics::SymMatrixD<real_t> tensor(number_);
	bool interp_result = interp(alpha,metrics,tensor);
	if (!interp_result) return false;

  // add the tensor along with its corresponding element to the field
  attachment_.add(tensor,elem);

  // note: this check will fail if the vertex is intended to be added after
  // its corresponding tensor is computed
  avro_assert( attachment_.check() );

	return true;
}

template<typename type>
bool
MetricField<type>::recompute( index_t p , real_t* x )
{
	avro_assert( p >= attachment_.points().nb_ghost() );

	Metric mp(number_);
	int ielem1 = interpolation_->eval( attachment_.points() , p , {attachment_[p].elem()} , mp );
	if (ielem1<0) return false;
	attachment_[p].set_elem( index_t(ielem1) );
	attachment_[p].calculate();
	return true;

	// look for the element in the background topology with the searcher
  index_t guess = attachment_[p].elem();
  int ielem = searcher_.find( x , guess );
  if (ielem<0)
  {
    // point is probably outside domain
		// let's make sure by first brute forcing the check
    ielem = searcher_.brute(x);
    if (ielem<0)
    {
			std::vector<real_t> alpha(number_+1,0.);
			std::vector<numerics::SymMatrixD<real_t>> metrics(alpha.size());
			ielem = searcher_.closest( x , alpha );
			index_t elem = index_t(ielem);
			for (index_t j=0;j<alpha.size();j++)
				metrics[j] = Field<type,Metric>::value(topology_(elem,j));

			bool interp_result = interp(alpha,metrics,attachment_[p]);
			if (!interp_result) return false;
			attachment_[p].set_elem(elem);
			attachment_[p].calculate();
			return true;
    }
  }

  index_t elem = index_t(ielem);
  const coord_t nv = topology_.nv(elem);

  // get the points of the element
  std::vector<const real_t*> xk(nv,0);
  for (index_t j=0;j<nv;j++)
    xk[j] = topology_.points()[ topology_(elem,j) ];
  real_t V = 1./numerics::simplex_volume(xk,topology_.points().dim());

  // compute the barycentric coordinates
  std::vector<real_t> alpha( nv , 0. );
	std::vector<numerics::SymMatrixD<real_t>> metrics(nv);
  for (index_t j=0;j<nv;j++)
  {
    std::vector<const real_t*> xk0(nv);
    for (index_t i=0;i<nv;i++)
    {
      if (i==j) xk0[i] = x;
      else xk0[i] = xk[i];
    }
    alpha[j] = V*numerics::simplex_volume(xk0,topology_.points().dim());

    if (alpha[j]<0. || alpha[j]>1.)
    {
			printf("alpha[%lu] = %1.16e\n",j,alpha[j]);
      printf("element search failed\n");
      avro_assert_not_reached;
    }

    // retrieve the metric at this vertex
		metrics[j] = Field<type,Metric>::value(topology_(elem,j));
  }

  // perform the interpolation
	bool interp_result = interp( alpha , metrics , attachment_[p] );
	if (!interp_result) return false;
	attachment_[p].set_elem(elem);
	attachment_[p].calculate();

	return true;
}

template<typename type>
bool
MetricField<type>::check_cells()
{
  index_t nb_violate = 0;
  for (index_t k=0;k<attachment_.points().nb();k++)
  {
    if (k<attachment_.points().nb_ghost()) continue;

    // get the element in the topology this vertex is inside
    index_t elem = attachment_[k].elem();
    real_t* x = attachment_.points()[k];
    avro_assert( !topology_.ghost(elem) );

    // ensure the barcyentric coordinates are ok
    bool inside = true;
    const coord_t nv = topology_.nv(elem);

    // get the points of the element
    std::vector<const real_t*> xk(nv,0);
    for (index_t j=0;j<nv;j++)
      xk[j] = topology_.points()[ topology_(elem,j) ];
    real_t V = 1./numerics::simplex_volume(xk,topology_.points().dim());

    // compute the barycentric coordinates
    std::vector<real_t> alpha( nv , 0. );
    for (index_t j=0;j<nv;j++)
    {
      std::vector<const real_t*> xk0(nv);
      for (index_t i=0;i<nv;i++)
      {
        if (i==j) xk0[i] = x;
        else xk0[i] = xk[i];
      }
      alpha[j] = V*numerics::simplex_volume(xk0,topology_.points().dim());

      if (alpha[j]<0.0 || alpha[j]>1.0)
      {
        inside = false;
      }
    }

    if (!inside)
    {
      nb_violate++;
      printf("vertex %lu is not in cell %lu!\n",k,elem);
      print_inline(alpha,"barycentric coords: ");
    }
  }
  if (nb_violate>0) printf("there are %lu violations\n",nb_violate);
  return nb_violate==0;
}

template<typename type>
void
MetricField<type>::remove( index_t k )
{
	// removes the metric stored at vertex k
  attachment_.remove(k);
}

void
edge_vector( const index_t i , const index_t j , const Points& v ,
              std::vector<real_t>& X )
{
  avro_assert( X.size()==v.dim() );
  for (coord_t d=0;d<v.dim();d++)
    X[d] = v[j][d] -v[i][d];
}

MetricAttachment::MetricAttachment( Points& points ) :
	number_(points.dim()),
	points_(points)
{}

MetricAttachment::MetricAttachment( Points& points , const std::vector<numerics::SymMatrixD<real_t>>& metrics ) :
  number_(metrics[0].n()),
	points_(points)
{
  avro_assert_msg( metrics.size()==points.nb() ,
									"number of metrics = %lu , number of points = %lu" ,
									metrics.size() , points.nb() );

	for (index_t k=0;k<points_.nb();k++)
  {
    if (points.ghost(k))
    {
			Metric mk(number_);
      for (index_t i=0;i<number_;i++)
        mk(i,i) = 1.;
      mk.calculate();
      Array<Metric>::add( mk );
    }
    else
		{
      Metric mk(number_);
			mk.set( metrics[k] );
      mk.calculate();
		  Array<Metric>::add( mk );
		}
  }
}

void
MetricAttachment::from_solb( const std::string& filename )
{
	avro_implement;
}

template<typename type>
void
MetricAttachment::set_cells( const Topology<type>& topology )
{
  std::vector<index_t> vertex(1);
  int offset = topology.points().nb_ghost() -points_.nb_ghost();
  avro_assert( points_.nb()+offset==topology.points().nb() );

	std::vector<bool> visited( points_.nb() , false );
	index_t counted = points_.nb_ghost();
	for (index_t k=0;k<points_.nb_ghost();k++)
		Array<Metric>::data_[k].set_elem(0); // doesn't matter, won't be used for interpolation

	for (index_t k=0;k<topology.nb();k++)
	{
		if (topology.ghost(k)) continue;
		for (index_t j=0;j<topology.nv(k);j++)
		{
			if (visited[topology(k,j)]) continue;

			counted++;
			Array<Metric>::data_[ topology(k,j) ].set_elem(k);
			visited[topology(k,j)] = true;
		}
	}
	avro_assert_msg(counted==points_.nb(),
									"counted = %lu, nb_points = %lu",counted,points_.nb());
}

template<typename type>
void
MetricAttachment::limit( const Topology<type>& topology , real_t href )
{
	const coord_t dim = topology.number();
	avro_assert_msg( topology.points().dim() == dim , "dim = %u , num = %u" , topology.points().dim() , dim );

	// the points should be associated with each other
	avro_assert_msg( topology.points().nb() == points_.nb() , "topology nb_points = %lu, points_.nb() = %lu" ,
									 topology.points().nb() , points_.nb() );

	// compute the implied metric of the input topology
	MeshImpliedMetric<type> implied( topology );
	implied.initialize();
	implied.optimize();

	// go through the entries of the current field and limit them using the step
	index_t nb_limited = 0;
	for (index_t k=0;k<nb();k++)
	{
		if (k < topology.points().nb_ghost()) continue;

		// compute the step from the implied metric to the current metric
		numerics::SymMatrixD<real_t> mi = implied[k];
		numerics::SymMatrixD<real_t> mt = this->operator[](k);

		numerics::SymMatrixD<real_t> invsqrt_M0 = numerics::powm(mi,-0.5);
		numerics::SymMatrixD<real_t> expS = invsqrt_M0*mt*invsqrt_M0;
		numerics::SymMatrixD<real_t> s = numerics::logm(expS);

		// limit the step
		bool limited = false;
		for (index_t i=0;i<dim;i++)
		for (index_t j=0;j<=i;j++)
		{
			if (s(i,j) > 2*std::log(href))
			{
				s(i,j) = 2*std::log(href);
				limited = true;
			}
			else if (s(i,j) < -2*std::log(href))
			{
				s(i,j) = -2*std::log(href);
				limited = true;
			}
		}
		if (limited) nb_limited++;
		else continue;

		// compute the new metric from the step
		numerics::SymMatrixD<real_t> sqrt_M0 = numerics::sqrtm(mi);
		numerics::SymMatrixD<real_t> mk = sqrt_M0*numerics::expm(s)*sqrt_M0;

		this->operator[](k).set(mk);
		this->operator[](k).calculate();
	}
	printf("limited %lu metrics out of %lu\n",nb_limited,this->nb());
}

void
MetricAttachment::reset( MetricAttachment& fld )
{
  Array<Metric>::clear();
  for (index_t k=0;k<fld.nb();k++)
    Array<Metric>::add(fld[k]);
}

void
MetricAttachment::add( numerics::SymMatrixD<real_t>& tensor , index_t elem )
{
	Metric mk(number_);
	mk.set(tensor);
	mk.set_elem(elem);
	mk.calculate();
	Array<Metric>::add(mk);
}

void
MetricAttachment::assign( index_t p , const numerics::SymMatrixD<real_t>& m0 , index_t elem )
{
	for (index_t j=0;j<number_;j++)
	for (index_t i=j;i<number_;i++)
		Array<Metric>::data_[p](i,j) = m0(i,j);
	Array<Metric>::data_[p].set_elem(elem);
	Array<Metric>::data_[p].calculate();
}

void
MetricAttachment::remove( index_t k , bool recheck )
{
	Array<Metric>::remove(k);
	if (recheck)
	{
	  avro_assert_msg( check() ,
	  "nb_points = %lu, nb_metrics = %lu" ,
	  points_.nb(),Array<Metric>::nb() );
	}
}

bool
MetricAttachment::check() const
{
	if (points_.nb()!=Array<Metric>::nb()) return false;
  return true;
}

void
MetricAttachment::to_json( json& J ) const
{
  J["name"] = "metric";
  J["membertype"] = "SPDT";
  J["evaltype"] = "vertex";
  J["numbertype"] = "real";
  J["nb"] = this->nb();
  index_t nb_rank = number_*(number_+1)/2;
  J["nb_rank"] = nb_rank;
  J["number"] = number_;
  std::vector<real_t> data(this->nb()*nb_rank);
  index_t n = 0;
  for (index_t k=0;k<this->nb();k++)
  for (index_t j=0;j<number_;j++)
	for (index_t i=j+1;i<number_;i++)
    data[n++] = (*this)[k](i,j);
  J["data"] = data;
}

#if 0
void
MetricAttachment::to_solb( const std::string& filename ) const
{
	int64_t fid;

	double buf[GmfMaxTyp];
	int dim = points_.dim();
	fid = GmfOpenMesh(filename.c_str(),GmfWrite,GmfDouble,dim);
	avro_assert( fid );

	int TypTab[GmfMaxTyp];
	TypTab[0] = GmfSymMat;
	GmfSetKwd( fid , GmfSolAtVertices , this->nb() , 1 , TypTab );

	for (index_t k=0;k<this->nb();k++)
	{
		const numerics::SPDT<real_t>& m = this->operator[](k);
		for (index_t j=0;j<m.nb();j++)
			buf[j] = m.data(j);
		GmfSetLin( fid , GmfSolAtVertices , buf );
	}
	GmfCloseMesh(fid);
}

void
MetricAttachment::from_solb( const std::string& filename )
{
	// open the file
	int dim,status;
	int nb_sol,numberType,solSize , TypTab[GmfMaxTyp];
	float fvalues[GmfMaxTyp];
	real_t dvalues[GmfMaxTyp];
	int version; // 2 for 32-bit int, 64-bit real
	int64_t fid;

	fid = GmfOpenMesh(filename.c_str(),GmfRead,&version,&dim);
	avro_assert_msg( fid , "could not open sol file %s ",filename.c_str() );

	printf("version = %d, dimension = %d\n",version,dim);

	// create a field whether this is attached at points or cells

	nb_sol = GmfStatKwd( fid , GmfSolAtVertices , &numberType , &solSize , TypTab );
	printf("nb_sol = %d, numberType = %d, solSize = %d\n",nb_sol,numberType,solSize);

	avro_assert( nb_sol == int(this->points_.nb()) );
	avro_assert( solSize == int(dim*(dim+1)/2) );

	avro_assert( GmfGotoKwd( fid , GmfSolAtVertices ) > 0 );
	for (int k=0;k<nb_sol;k++)
	{

		// read the metric
		if (version==1)
		{
			status = GmfGetLin( fid , GmfSolAtVertices , fvalues );
			for (index_t j=0;j<6;j++)
				dvalues[j] = real(fvalues[j]);
		}
		else
			status = GmfGetLin( fid , GmfSolAtVertices , dvalues );

		avro_assert( status==1 );

		std::vector<real_t> data(dvalues,dvalues+solSize);
		numerics::SPDT<real_t> m(data);
		this->operator[](k) = m;
	}
}
#endif

template class MetricField<Simplex>;
template void MetricAttachment::set_cells( const Topology<Simplex>& );
template void MetricAttachment::limit(const Topology<Simplex>&,real_t);

} // avro
