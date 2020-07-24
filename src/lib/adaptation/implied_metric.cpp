//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "adaptation/implied_metric.h"
#include "adaptation/metric.h"

#include "common/tools.h"

#include "element/simplex.h"

#include "mesh/points.h"
#include "mesh/topology.h"

#include "numerics/linear_algebra.h"
#include "numerics/nlopt_result.h"

#include <numpack/types/SurrealS.h>
#include <numpack/types/PromoteSurreal.h>

#include <nlopt.hpp>

#include <cmath>

namespace avro
{

template<typename type>
class JacobianEquilateral : public numerics::MatrixD<real_t>
{
public:
  JacobianEquilateral( coord_t number ) :
    numerics::MatrixD<real_t>(number,number)
  {
    numerics::MatrixD<real_t>& J = *this;
    if (number==1) J(0,0) = 1.0;
    else if (number==2)
    {
      J(0,0) = 1.0;
      J(0,1) = -1.0/sqrt(3.0);
      J(1,0) = 0.0;
      J(1,1) = 2.0/sqrt(3.0);
    }
    else if (number==3)
    {
      J(0,0) = 1.0; J(0,1) = -1.0/sqrt(3.0); J(0,2) = -1.0/sqrt(6.0);
      J(1,0) = 0.0; J(1,1) = 2.0/sqrt(3.0) ; J(1,2) = -1.0/sqrt(6.0);
      J(2,0) = 0.0; J(2,1) = 0.0           ; J(2,2) = 3.0/sqrt(6.0);
    }
    else if (number==4)
    {
      J(0,0) = (sqrt(2.0)*sqrt(3.0)*sqrt(5.0)*-2.0)/(sqrt(5.0)*sqrt(1.5E1)*3.0+sqrt(3.0)*5.0);
      J(0,1) = (sqrt(2.0)*3.0E1)/(sqrt(5.0)*sqrt(1.5E1)*3.0+sqrt(3.0)*5.0);
      J(0,2) = 0.0;
      J(0,3) = 0.0;
      J(1,0) = (sqrt(2.0)*sqrt(1.5E1)*-2.0)/(sqrt(5.0)*sqrt(1.5E1)*3.0+sqrt(3.0)*5.0);
      J(1,1) = (sqrt(2.0)*-1.0E1)/(sqrt(5.0)*sqrt(1.5E1)*3.0+sqrt(3.0)*5.0);
      J(1,2) = sqrt(2.0)*sqrt(6.0)*(1.0/3.0);
      J(1,3) = 0.0;
      J(2,0) = (sqrt(2.0)*sqrt(1.5E1)*-2.0)/(sqrt(5.0)*sqrt(1.5E1)*3.0+sqrt(3.0)*5.0);
      J(2,1) = (sqrt(2.0)*-1.0E1)/(sqrt(5.0)*sqrt(1.5E1)*3.0+sqrt(3.0)*5.0);
      J(2,2) = sqrt(2.0)*sqrt(6.0)*(-1.0/6.0);
      J(2,3) = 1.0;
      J(3,0) = (sqrt(2.0)*sqrt(1.5E1)*-2.0)/(sqrt(5.0)*sqrt(1.5E1)*3.0+sqrt(3.0)*5.0);
      J(3,1) = (sqrt(2.0)*-1.0E1)/(sqrt(5.0)*sqrt(1.5E1)*3.0+sqrt(3.0)*5.0);
      J(3,2) = sqrt(2.0)*sqrt(6.0)*(-1.0/6.0);
      J(3,3) = -1.0;
    }
    else
      avro_implement;
  }
};

template<typename type>
ElementImpliedMetric<type>::ElementImpliedMetric( const type& element ) :
  numerics::SymMatrixD<real_t>(element.number()),
  element_(element),
  J_( element.number() , element.number() ),
  J0_( element.number() , element.number() ),
  Jeq_(JacobianEquilateral<type>(element.number())),
  M_( element.number() , element.number() )
{
  detJeq_ = numerics::determinant(Jeq_);
}

template<typename type>
void
ElementImpliedMetric<type>::compute( const std::vector<const real_t*>& xk )
{
  element_.jacobian( xk , J0_ );
  J_ = J0_*Jeq_;
  numerics::SymMatrixD<real_t> JJt = J_*numpack::Transpose(J_);
  M_ = numerics::inverse( JJt );
  for (index_t i=0;i<element_.number();i++)
  for (index_t j=i;j<element_.number();j++)
    this->operator()(i,j) = M_(i,j);
}

template<typename type>
void
ElementImpliedMetric<type>::compute( const Points& points , const index_t* v , index_t nv )
{
  element_.jacobian( v , nv , points , J0_ );
  J_ = J0_*Jeq_;
  numerics::SymMatrixD<real_t> JJt = J_*numpack::Transpose(J_);
  M_ = numerics::inverse( JJt );
  for (index_t i=0;i<element_.number();i++)
  for (index_t j=i;j<element_.number();j++)
    this->operator()(i,j) = M_(i,j);
}


template<typename type>
void
ElementImpliedMetric<type>::inverse( const Points& points , const index_t *v , index_t nv )
{
  element_.jacobian( v , nv , points , J0_ );
  J_ = J0_*Jeq_;
  M_ = J_*numpack::Transpose(J_); // no inverse
  for (index_t i=0;i<element_.number();i++)
  for (index_t j=i;j<element_.number();j++)
    this->operator()(i,j) = M_(i,j);
}

template<typename type>
real_t
ElementImpliedMetric<type>::determinant( const std::vector<const real_t*>& xk )
{
  element_.jacobian( xk , J0_ );
  real_t detJ0 = numerics::determinant(J0_);
  if (detJ0==0.0) return 0.0;
  return 1./( detJ0*detJ0*detJeq_*detJeq_ );
}

template<typename type>
real_t
ElementImpliedMetric<type>::determinant( const Points& points , const index_t* v , const index_t nv )
{
  element_.jacobian( v, nv, points , J0_ );
  real_t detJ0 = numerics::determinant(J0_);
  if (detJ0==0.0) return 0.0;
  return 1./( detJ0*detJ0*detJeq_*detJeq_ );
}

template<typename type>
MeshImpliedMetric<type>::MeshImpliedMetric( const Topology<type>& topology ) :
  topology_(topology)
{
  numerics::SymMatrixD<real_t> zero( topology_.element().number() );
  this->resize( topology.points().nb() , zero );
  nodalMetricSqrt_.resize( topology_.points().nb() , zero );
  nodalMetricSqrtDet_.resize( topology_.points().nb() , 0. );
  edges_.clear();
  topology_.get_edges(edges_);

  // keep track of the vertex-to-element relations
  v2e_.resize( topology_.points().nb() );
  for (index_t k=0;k<topology_.nb();k++)
  {
    if (topology_.ghost(k)) continue;
    for (index_t j=0;j<topology_.nv(k);j++)
    {
      v2e_[ topology_(k,j) ].push_back(k);
    }
  }
}

template<typename type>
const Topology<type>&
MeshImpliedMetric<type>::topology() const
{
  return topology_;
}

template<typename type>
void
MeshImpliedMetric<type>::initialize()
{
  std::vector<real_t> alpha;
  std::vector<const real_t*> xj;
  std::vector<real_t> volk(topology_.nb());
  std::vector<numerics::SymMatrixD<real_t>> metrics( topology_.nb() , numerics::SymMatrixD<real_t>(topology_.number()) );
  std::vector<numerics::SymMatrixD<real_t>> mb;

  // v0 is not used here, but for the local volume sum in the vertex loop
  topology_.get_volumes( volk );

  // pre-compute the implied metric of each element
  for (index_t k=0;k<topology_.nb();k++)
  {
    if (topology_.ghost(k))
    {
      numerics::SymMatrixD<real_t> zero( topology_.element().number() );
      metrics[k] = zero;
      continue;
    }
    topology_.get_elem( k , xj );
    ElementImpliedMetric<type> mk( topology_.element() );
    //mk.compute( xj );
    mk.compute( topology_.points() , topology_(k) , topology_.nv(k) );
    metrics[k] = mk;
  }

  // initialize the vertex implied metric as the volume-weighted average
  // of the element metrics in the ball of the vertex
  for (index_t k=0;k<topology_.points().nb();k++)
  {
    if (k<topology_.points().nb_ghost()) continue;

    // retrieve the ball of the vertex
    const std::vector<index_t>& B = v2e_[k];

    // compute the weights (from the volumes)
    alpha.resize( B.size() , 0.0 );
    real_t v0 = 0.0;
    for (index_t j=0;j<B.size();j++)
    {
      if (topology_.ghost(B[j])) continue;
      alpha[j] = volk[B[j]];
      v0 += alpha[j];
    }
    avro_assert_msg( v0 > 0 , "topology may need to be oriented" );

    // normalize the average and retrieve the element metrics
    mb.resize( B.size() );
    for (index_t j=0;j<B.size();j++)
    {
      alpha[j] /= v0;
      mb[j] = metrics[ B[j] ];
    }

    // compute the weighted average
    try
    {
      interp( alpha , mb , this->data_[k] );
    }
    catch(...)
    {
      print_inline(alpha);
      printf("v0 = %g\n",v0);
      avro_implement;
    }

    nodalMetricSqrt_[k]    = numerics::sqrtm(this->data_[k]);
    nodalMetricSqrtDet_[k] = numerics::determinant(nodalMetricSqrt_[k]);
  }
}

template<typename type>
template<int DIM>
real_t
MeshImpliedMetric<type>::cost( const std::vector<numerics::SymMatrixD<real_t>>& sv ,
                               std::vector<numerics::SymMatrixD<real_t>>& dc_dS ,
                               real_t& complexity0 ) const
{
  // precompute volume
  std::vector<real_t> volk(topology_.nb());
  topology_.get_volumes( volk );

  index_t nrank = DIM*(DIM+1)/2;
  real_t  one_over_nv = 1./real_t(DIM+1);

  // reference complexity
  complexity0 = topology_.nb_real()*topology_.element().reference().vunit();

  avro_assert( sv.size()==topology_.points().nb() );
  if (dc_dS.size()!=0)
    avro_assert( dc_dS.size()==topology_.points().nb() );

  // compute elemental step matrices and elemental costs
  std::vector< numerics::SymMatrixD<real_t> > sk( topology_.nb() , numerics::SymMatrixD<real_t>(DIM) );
  std::vector<real_t> ck( topology_.nb() , 0. );
  real_t complexity = 0.0;
  for (index_t k=0;k<topology_.nb();k++)
  {
    if (topology_.ghost(k)) continue;

    real_t sqrtDetM = -1;
    for (index_t j=0;j<topology_.nv(k);j++)
    {
      // use the maximum nodal metric determinant to estimate the original cost
      if (nodalMetricSqrtDet_[ topology_(k,j) ] > sqrtDetM)
        sqrtDetM = nodalMetricSqrtDet_[ topology_(k,j) ];

      // arithmetic average sv to form sk
      for (index_t i=0;i<nrank;i++)
        sk[k].data[i] += sv[ topology_(k,j) ].data[i]*one_over_nv;
    }

    // compute the contribution from this element
    // the first term is caused by m0 and the second term is caused by the step
    for (index_t i=0;i<DIM;i++)
      ck[k] += sk[k](i,i);
    ck[k] = sqrtDetM*volk[k] * exp(0.5*ck[k]);
    //ck[k] = exp(0.5*ck[k]);

    // add the contribution to the complexity
    complexity += ck[k];
  }

  // option to compute the gradients
  if (dc_dS.size()!=0)
  {
    for (index_t k=0;k<topology_.points().nb();k++)
    {
      // zero out the derivative
      for (index_t i=0;i<nrank;i++)
        dc_dS[k].data[i] = 0.0;

      // add the contribution from the surrounding elements
      for (index_t j=0;j<v2e_[k].size();j++)
      {
        for (index_t i=0;i<DIM;i++)
          dc_dS[k](i,i) += ck[ v2e_[k][j] ]*0.5*one_over_nv;
      }
    }
  }

  return complexity;
}

template<typename T>
static T
power( const T& x , index_t p )
{
  T result = 1;
  for (index_t k=1;k<=p;k++)
    result *= x;
  return result;
}

// Smooth max with C2 continuity that matches at x = +- eps/2
template<typename T>
static T
smoothmaxC2(const T x, const T y, const real_t eps = 1e-2 )
{
  const T m = x > y ? y : x;
  const T M = x > y ? x : y;
  const T d = M-m; // max(x,y) = y + max(0,x-y);

  // from solving linear system with 6 boundary conditions
  return ( d > eps / 2 ) ? M : m + (power(eps+2*d,3)*(3*eps - 2*d))/(32*power(eps,3));
}

template<typename type>
template<int DIM>
real_t
MeshImpliedMetric<type>::deviation( const std::vector<numerics::SymMatrixD<real_t>>& Svec ,
                                           std::vector<numerics::SymMatrixD<real_t>>& df_dS ) const
{
  typedef SurrealS<DIM*(DIM+1)/2> SurrealClassVertex;
  typedef SurrealS<2*DIM*(DIM+1)/2> SurrealClassEdge;
  typedef numerics::SymMatrixD<SurrealClassVertex> MatrixSymSurrealVertex;

  index_t nrank = DIM*(DIM+1)/2;

  real_t delta = 0.0;

  avro_assert( Svec.size()==topology_.points().nb() );
  if (df_dS.size()!=0)
    avro_assert( df_dS.size()==topology_.points().nb() );

  // set the nodal metric derivatives
  std::vector<MatrixSymSurrealVertex> nodalMetric( topology_.points().nb() );
  for (index_t k=0;k<topology_.points().nb();k++)
  {
    MatrixSymSurrealVertex S = Svec[k];
    for (index_t i=0;i<nrank;i++)
      S.data[i].deriv(i) = 1.0; // derivative with respect to its own value is 1

    MatrixSymSurrealVertex sqrtM0 = nodalMetricSqrt_[k];
    nodalMetric[k] = sqrtM0*numerics::expm(S)*sqrtM0;

  }

  // compute the deviation for every edge
  std::vector<real_t> dx( topology_.points().dim() );
  index_t nb_edges = edges_.size()/2;
  index_t nb_violate = 0;
  for (index_t k=0;k<nb_edges;k++)
  {
    index_t p = edges_[2*k];
    index_t q = edges_[2*k+1];

    if (p < topology_.points().nb_ghost() || q < topology_.points().nb_ghost())
      continue;

    // get the edge vector and copy to surreal_t version
    numerics::vector( topology_.points()[p] , topology_.points()[q] ,
                      topology_.points().dim() , dx.data() );

    Entity* entity = BoundaryUtils::geometryFacet( topology_.points() , edges_.data()+2*k , 2 );
    topology_.element().edge_vector( topology_.points() , p , q , dx.data() , entity );

    // get the edge length squared
    numerics::VectorD<real_t> e( DIM , dx.data() );
    SurrealClassVertex lni = quadratic_form(nodalMetric[p],e);//numpack::Transpose(e)*nodalMetric[p]*e;
    SurrealClassVertex lnj = quadratic_form(nodalMetric[q],e);//numpack::Transpose(e)*nodalMetric[q]*e;
    lni = sqrt(lni);
    lnj = sqrt(lnj);

    // edge length under each nodal metric
    SurrealClassEdge lmi = lni.value();
    SurrealClassEdge lmj = lnj.value();

    // save the derivatives from each nodal metric
    for (index_t d=0;d<nrank;d++)
    {
      lmi.deriv(d)       = lni.deriv(d);
      lmj.deriv(nrank+d) = lnj.deriv(d);
    }

    SurrealClassEdge edgeLen = 0;

    if (fabs(lmi-lmj) <= 1.0e-6*lmi )
      edgeLen = lmi;
    else
      edgeLen = (lmi -lmj)/(log(lmi/lmj));

    // compute the contribution to the objective function
    SurrealClassEdge del = 0;
    if (edgeLen>1)
      del = edgeLen -sqrt(2.);
    else
      del = 1./sqrt(2.) -edgeLen;
    if (del.value()>0) nb_violate++;
    del = smoothmaxC2( SurrealClassEdge(0.0) , del , 1e-2 );
    SurrealClassEdge delUnity = power( del , topology_.number() );

    // add the contribution to the objective function
    delta += delUnity.value()/nb_edges;

    // no gradients requested
    if (df_dS.size()==0) continue;

    // store the gradients
    for (index_t i=0;i<nrank;i++)
    {
      df_dS[p].data[i] += delUnity.deriv(i)/nb_edges;
      df_dS[q].data[i] += delUnity.deriv(nrank+i)/nb_edges;
    }
  }
  return delta;
}

template<typename type>
struct nlopt_data
{
	MeshImpliedMetric<type>& metric;
	index_t eval_count;
	real_t objective;
  real_t complexity;
  real_t deviation;
  bool include_complexity;
};

template<int DIM,typename type>
double
impliedMetric_objective( unsigned n , const double* x , double* grad, void* data0 )
{
	nlopt_data<type>* data = static_cast<nlopt_data<type>*>(data0);
	MeshImpliedMetric<type>& metric = data->metric;
  const Topology<type>& topology = metric.topology();

  // increase the evaluation count
	data->eval_count++;

  // rank of each step matrix
	const index_t nrank = DIM*(DIM+1)/2;

  // fill in the step matrices
	std::vector<numerics::SymMatrixD<real_t>> S( topology.points().nb() , numerics::SymMatrixD<real_t>(DIM) );
	for (index_t k=0;k<topology.points().nb();k++)
	{
		for (index_t i=0;i<nrank;i++)
			S[k].data[i] = x[k*nrank+i];
	}

  // size the gradients if necessary
	std::vector<numerics::SymMatrixD<real_t>> dl_dS;
  std::vector<numerics::SymMatrixD<real_t>> dc_dS;
	if (grad)
	{
    numerics::SymMatrixD<real_t> zero(DIM);
    zero = 0;
		dl_dS.resize( topology.points().nb() , zero );
    dc_dS.resize( topology.points().nb() , zero );
	}

  // compute the edge length deviation
	real_t deviation = metric.template deviation<DIM>( S , dl_dS );

  // compute the complexity
  real_t complexity0;
  real_t complexity = metric.template cost<DIM>( S , dc_dS , complexity0 );

  real_t factor = 1;

  // fill the gradients
  real_t gradnorm = 0.0;
	if (grad)
	{
		for (index_t k=0;k<topology.points().nb();k++)
		{
			for (index_t i=0;i<nrank;i++)
			{
				grad[k*nrank+i] = dl_dS[k].data[i];
        if (data->include_complexity)
          grad[k*nrank+i] += 2*factor*(complexity/complexity0-1.0)*dc_dS[k].data[i]*(1./complexity0);
				gradnorm += pow(grad[k*nrank+i],2.);
			}
		}
		gradnorm = sqrt(gradnorm/n);
	}

  // compute and save the objective function
  real_t obj = 0.0;
  obj = deviation;
  if (data->include_complexity)
    obj += factor*pow( complexity/complexity0 - 1.0 , 2.0 );
	data->objective = obj;
  data->complexity = complexity;
  data->deviation  = deviation;

  //printf("|g| = %e, obj = %e (d = %g, c = %g)\n",gradnorm,obj,data->deviation,data->complexity/complexity0);

  return obj;
}

template<typename type>
void
MeshImpliedMetric<type>::optimize()
{
  // compute the number of variables
  coord_t n = topology_.points().dim();
  index_t nrank = n*(n+1)/2;
  index_t N = topology_.points().nb()*nrank;

  // setup the optimizer
  nlopt::opt opt( nlopt::LD_MMA , N );

  // initialize the step matrices to zero
	std::vector<real_t> x( N );
	for (index_t k=0;k<topology_.points().nb();k++)
	for (index_t i=0;i<nrank;i++)
		x[k*nrank+i] = 0.0;

	// assign the data used by the nlopt objective function
	nlopt_data<type> data = {*this,0,1,1,1,true};

	// set the objective function
	if (topology_.number()==1)
		opt.set_min_objective( &impliedMetric_objective<1,type> , static_cast<void*>(&data) );
	else if (topology_.number()==2)
		opt.set_min_objective( &impliedMetric_objective<2,type> , static_cast<void*>(&data) );
	else if (topology_.number()==3)
		opt.set_min_objective( &impliedMetric_objective<3,type> , static_cast<void*>(&data) );
	else if (topology_.number()==4)
		opt.set_min_objective( &impliedMetric_objective<4,type> , static_cast<void*>(&data) );
	else
		avro_implement;

	// set the lower and upper bounds on the entries of the step matrix
	std::vector<real_t> lower_bound( N , -2*log(2.) );
	std::vector<real_t> upper_bound( N ,  2*log(2.) );
	opt.set_lower_bounds(lower_bound);
	opt.set_upper_bounds(upper_bound);

	// set some optimization parameters
  opt.set_stopval(1e-7);
  opt.set_xtol_rel(1e-7);
  opt.set_ftol_rel(1e-7);
  opt.set_maxeval(200);

	// optimize the step matrices
  printf("optimizing implied metric...\n");
	real_t f_opt;
	nlopt::result result = nlopt::result::SUCCESS;
  result = opt.optimize( x , f_opt );
  printf("done!\n");
  printf("summary:\n");
	printf("\tresult = %s\n",nloptResultDescription(result).c_str());
	printf("\teval_count = %lu, objective = %3.6e\n",data.eval_count,data.objective);
  printf("\tcomplexity = %g\n",data.complexity);
  printf("\tdeviation = %g\n",data.deviation);

	// store the result
	for (index_t k=0;k<topology_.points().nb();k++)
	{

		numerics::SymMatrixD<real_t> S( topology_.number() );
		for (index_t i=0;i<nrank;i++)
			S.data[i] = x[k*nrank+i];

    // assign the implied metric
    this->data_[k] = nodalMetricSqrt_[k]*numerics::expm(S)*nodalMetricSqrt_[k];

    if (k<topology_.points().nb_ghost())
    {
      this->data_[k] = 0;
      continue;
    }

    real_t detm = numerics::determinant(this->data_[k]);
    if (detm <= 0.0)
    {
      // forget about the step
      //this->data_[k].dump();
      this->data_[k] = nodalMetricSqrt_[k]*nodalMetricSqrt_[k];

      std::vector<index_t> ball;
      topology_.inverse().ball( k , ball );
      for (index_t j=0;j<ball.size();j++)
      for (index_t i=0;i<topology_.nv(ball[j]);i++)
        topology_.points().set_fixed( topology_(ball[j],i) , true );

      printf("detm = %g! -> fixing ball of vertex %lu",detm,k);
    }
    //avro_assert_msg( detm > 0. , "detm = %g for vertex %lu",detm,k );
	}
}

template class ElementImpliedMetric<Simplex>;
template class MeshImpliedMetric<Simplex>;

} // avro
