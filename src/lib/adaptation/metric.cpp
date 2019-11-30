#include "adaptation/metric.h"

#include "mesh/search.h"
#include "mesh/topology.h"
#include "mesh/points.h"

#include "numerics/determinant.h"
#include "numerics/functions.h" // factorial
#include "numerics/geometry.h"

#if 0
#include <json/json.hpp>
#include <libMeshb/libMeshb7/libmeshb.h>
#endif

#include <cmath>

namespace luna
{

#if 0
template<typename type>
MetricField<type>::MetricField( coord_t _n , Topology<type>& _topology , MetricAttachment& fld ) :
	topology_(_topology),
	field_(_field),
  searcher_(topology_),
{
	luna_assert( field_.nb()==topology_.vertices().nb() );
	luna_assert( field_.nb()>0 );
	this->n_ = field_[0].n();

	if (this->n_==2)
	{
		normalization_ = 4.*std::sqrt(3.);
	}
	else if (this->n_==3)
	{
		normalization_ = 36./std::pow(3.,1./3.);
	}
	else if (this->n_==4)
	{
		normalization_ = 10./std::pow( topology_.master().unitVolume() ,
																	 2./topology_.number() );
	}
	else
		luna_assert_not_reached;

	if (topology_.number()==4)
		normalization_ = 10./std::pow( topology_.master().unitVolume() , 2./4. );

  // save the tensors stored in the field
  // these are static and need to be kept for interpolation
  for (index_t k=0;k<field_.nb();k++)
    tensors_.add( field_[k] );
}

template<typename type>
MetricField<type>::MetricField( Topology<type>& _topology , MetricAttachment& _field ) :
	MetricField<type>( _topology.number() , _topology , _field )
{}

template<typename type>
void
MetricField<type>::reset( MetricAttachment& fld )
{
  tensors_.clear();
  field_.reset( fld );

  for (index_t k=0;k<fld.nb();k++)
    tensors_.add( fld[k] );
}

template<typename type>
SPDT<real>&
MetricField<type>::operator() ( const Vertices& x , index_t v )
{
  luna_assert_msg( v<field_.nb() ,
                   "v = %lu but field has %lu tensors",v,field_.nb() );
  return field_[v];
}

template<typename type>
real
MetricField<type>::length( index_t n0 , index_t n1 )
{
  luna_assert_msg( n0 < field_.nb() ,
                  "n0 = %lu, field_.nb() = %lu" , n0, field_.nb() );
  luna_assert_msg( n1 < field_.nb() ,
                  "n1 = %lu, field_.nb() = %lu" , n1, field_.nb() );
  real l = this->geometricInterpolation( field_.vertices() , n0 , n1 );

  // limit the edge length change
  if (limits_[0]>0 && l < limits_[0]) return limits_[0];
  if (limits_[1]>0 && l > limits_[1]) return limits_[1];
  return l;
}

template<typename type>
real
MetricField<type>::volume( const Topology<type>& topology , const index_t k )
{
	const index_t *V = topology(k);
	const index_t NV = topology.nv(k);
	const type& master = topology.master();

	if (topology.ghost(k)) return 0.0;

	// find the metric with the maximum determinant
	index_t jmax = 0;
	real d = field_.sqrtDetM(V[jmax]);
	real dmax = d;
	for (coord_t j=1;j<NV;j++)
	{
		d = field_.sqrtDetM(V[j]);
		if (d>dmax)
		{
			dmax = d;
			jmax = j;
		}
	}

	real sqrtdetM = dmax;
	real v = sqrtdetM*master.volume(topology.vertices(),V,NV);
	return v;
}

template<typename type>
real
MetricField<type>::volume( const Topology<type>& t )
{
  real v = 0.;
  for (index_t k=0;k<t.nb();k++)
  {
    if (t.ghost(k)) continue;
    v += volume(t,k);
  }
  return v;
}

template<typename type>
real
MetricField<type>::volume0()
{
  return volume(topology_);
}

template<typename type>
real
MetricField<type>::quality( const Topology<type>& topology , index_t k )
{
	const index_t *V = topology(k);
	const index_t NV = topology.nv(k);
	const Vertices& vertices = topology.vertices();
	const coord_t dim = vertices.dim();
	const coord_t num = topology.number();
	const type& master = topology.master();

	if (topology.ghost(k)) return -1.;

	// find the metric with the maximum determinant
	index_t jmax = 0, jmin = 0;
	real d = field_.sqrtDetM(V[jmax]);
	real dmax = d;
	real dmin = d;
	for (coord_t j=1;j<NV;j++)
	{
		d = field_.sqrtDetM(V[j]);
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

	real sqrtdetM = dmax;
	index_t idxM = jmax;
	UNUSED(jmin);

	// tensor with maximum determinant
	const SPDT<real>& M = field_[ V[idxM] ];

	// compute the edge lengths under m
  real l = 0.,lj;
	std::vector<real> e( dim , 0. );
  for (index_t j=0;j<master.nb_Edges;j++)
  {
    index_t p0 = master.edge(j,0);
    index_t p1 = master.edge(j,1);

		geometrics::vector( vertices[V[p0]] , vertices[V[p1]] , dim , e.data() );
		lj = M.quadraticForm( e.data() );
    l += lj;
  }

	// compute the volume under m
	real v = sqrtdetM*master.volume(topology.vertices(),V,NV);
	luna_assert_msg( v>=0. , "v = %g, sqrtDetM = %g, v_e = %g" , v , sqrtdetM , master.volume(topology.vertices(),V,NV) );
	v = std::pow( v , 2./num );

	// normalize to be within [0,1]
	real q = normalization_*v/l;
	return q;
}

template<typename type>
void
MetricField<type>::initializeCells()
{
  // assume the vertices in the topology are associated with the metric
  // with the VertexField
  // this means we only need to find a cell containing each vertex
  std::vector<index_t> vertex(1);

  for (index_t k=0;k<topology_.vertices().nb();k++)
  {
    vertex[0] = k;
    std::vector<index_t> elems;
    topology_.allWithSubset(vertex,elems);
    luna_assert( elems.size()>0 );

    if (k<topology_.vertices().nb_ghost())
    {
      field_.cell(k) = elems[0];
      continue;
    }

    bool found = false;
    for (index_t j=0;j<elems.size();j++)
    {
      if (topology_.ghost(elems[j])) continue;
      field_.cell(k) = elems[j];
      found = true;
      break;
    }
    luna_assert(found);
  }
}

template<typename type>
template<typename shape>
bool
MetricField<type>::check( Topology<shape>& topology )
{
  if (topology.vertices().nb()!=field_.nb()) return false;
  if (topology.vertices().nb()!=field_.cell().size()) return false;
  if (topology.vertices().nb()!=field_.nb()) return false;
  return true;
}

template<typename type>
int
MetricField<type>::find( index_t n0 , index_t n1 , real* x )
{
  luna_assert_msg( n0>=topology_.vertices().nb_ghost() , "n0 = %lu, nb_ghost = %lu" , n0 , topology_.vertices().nb_ghost() );
  luna_assert_msg( n1>=topology_.vertices().nb_ghost() , "n1 = %lu, nb_ghost = %lu" , n1 , topology_.vertices().nb_ghost() );

  luna_assert( n0 < field_.vertices().nb() );
  luna_assert( n1 < field_.vertices().nb() );

  luna_assert_msg( n0 < field_.nb() , "n0 = %lu , field_.nb = %lu\n" , n0 , field_.nb() );
  luna_assert_msg( n0 < field_.cell().size() , "n0 = %lu , field_nb = %lu\n" , n0 , field_.nb_cell() );
  luna_assert( n1 < field_.cell().size() );

  index_t guess = field_.cell(n0);
  luna_assert_msg( guess < topology_.nb() , "guess = %lu but nb = %lu\n",guess,topology_.nb());
  int elem = searcher_.find( x , guess );
  if (elem>=0) return elem;
  guess = field_.cell(n1);
  luna_assert_msg( guess < topology_.nb() , "guess = %lu but nb = %lu\n",guess,topology_.nb());
  elem = searcher_.find( x , guess );
  if (elem>=0) return elem;
  elem = searcher_.brute(x);
  return elem;
}

template<typename type>
void
MetricField<type>::interpolate( real* x , index_t elem , SPDT<real>& tensor , bool STFU )
{
  const coord_t nv = topology_.nv(elem);

  // get the vertices of the element
  std::vector<real*> xk(nv,0);
  for (index_t j=0;j<nv;j++)
    xk[j] = topology_.vertices()[ topology_(elem,j) ];
  real V = 1./geometrics::signedVolume(xk,topology_.vertices().dim());

  // compute the barycentric coordinates
  std::vector<real> alpha( nv , 0. );
  for (index_t j=0;j<nv;j++)
  {
    std::vector<real*> xk0(nv);
    for (index_t i=0;i<nv;i++)
    {
      if (i==j) xk0[i] = x;
      else xk0[i] = xk[i];
    }
    alpha[j] = V*geometrics::signedVolume(xk0,topology_.vertices().dim());

    if (alpha[j]<0. || alpha[j]>1.)
    {
      if (!STFU)
      {
        printf("element search failed!, alpha[%lu] = %g\n",j,alpha[j]);
        luna_assert_not_reached;
      }
    }

    // save the metric at this vertex
    this->Melem_[j] = tensors_[ topology_(elem,j) ];
  }

  // perform the interpolation
  tensor.interpolate<LogEuclidean<real>>(alpha,this->Melem_);
}

template<typename type>
void
MetricField<type>::add( index_t n0 , index_t n1 , real* x )
{
  // find the element containing x bordering n0 and n1
  int ielem = find(n0,n1,x);

  // check if the containing element was not found (point outside domain)
  // this can only happen with curved geometries
  if (ielem<0)
  {
		//luna_assert_not_reached; // pcaplan temporary
		SPDT<real> tensor( this->n_ );
		#if 0
		std::vector<real> alpha(2,0.5);
		tensor.interpolate<LogEuclidean<real>>(alpha, {field_[n0],field_[n1]} );
    field_.add( tensor , field_.cell(n0) );
		#else
		std::vector<real> alpha(this->n_+1,0.);
		ielem = searcher_.closest( x , alpha );
		index_t elem = index_t(ielem);
		for (index_t j=0;j<alpha.size();j++)
			this->Melem_[j] = tensors_[topology_(elem,j)];
		tensor.interpolate<LogEuclidean<real>>(alpha,this->Melem_);
		field_.add( tensor , elem );
		#endif
		return;
  }

  index_t elem = index_t(ielem);
  const coord_t nv = topology_.nv(elem);

  // get the vertices of the element
  std::vector<real*> xk(nv,0);
  for (index_t j=0;j<nv;j++)
    xk[j] = topology_.vertices()[ topology_(elem,j) ];
  real V = 1./geometrics::signedVolume(xk,topology_.vertices().dim());

  // compute the barycentric coordinates
  std::vector<real> alpha( nv , 0. );
  for (index_t j=0;j<nv;j++)
  {
    std::vector<real*> xk0(nv);
    for (index_t i=0;i<nv;i++)
    {
      if (i==j) xk0[i] = x;
      else xk0[i] = xk[i];
    }
    alpha[j] = V*geometrics::signedVolume(xk0,topology_.vertices().dim());

    if (alpha[j]<0. || alpha[j]>1.)
    {
      printf("element search failed\n");
      luna_assert_not_reached;
    }

    // save the metric at this vertex
    this->Melem_[j] = tensors_[ topology_(elem,j) ];
  }

  // perform the interpolation
  SPDT<real> tensor( this->n_ );
  tensor.interpolate<LogEuclidean<real>>(alpha,this->Melem_);

  // add the tensor along with its corresponding element to the field
  field_.add(tensor,elem);

  // note: this check will fail if the vertex is intended to be added after
  // its corresponding tensor is computed
  luna_assert( field_.check() );
}

template<typename type>
void
MetricField<type>::recompute( index_t p , real* x ,
                                 const std::vector<index_t>& N )
{
  // the vertex at p has new coordinates so the element containing it
  // (within the background topology) may have changed
  // check the elements attached to the ball (the vertex needs to be in one of
  // these cells)
  index_t elem;
  std::vector<index_t> elems;

  printf("recomputing metric for vertex %lu\n",p);
  printInline( std::vector<real>(x,x+topology_.vertices().dim()), "coords: ");

  // list all the elements to try
  for (index_t k=0;k<N.size();k++)
    elems.push_back( field_.cell( N[k] ) );

  bool inside = true;
  for (index_t k=0;k<elems.size();k++)
  {
    inside = true;
    elem = elems[k];
    if (topology_.ghost(elem)) continue;
    const coord_t nv = topology_.nv(elem);

    // get the vertices of the element
    std::vector<real*> xk(nv,0);
    for (index_t j=0;j<nv;j++)
      xk[j] = topology_.vertices()[ topology_(elem,j) ];
    real V = 1./geometrics::signedVolume(xk,topology_.vertices().dim());

    // compute the barycentric coordinates
    std::vector<real> alpha( nv , 0. );
    for (index_t j=0;j<nv;j++)
    {
      std::vector<real*> xk0(nv);
      for (index_t i=0;i<nv;i++)
      {
        if (i==j) xk0[i] = x;
        else xk0[i] = xk[i];
      }
      alpha[j] = V*geometrics::signedVolume(xk0,topology_.vertices().dim());

      if (alpha[j]==0.0 || alpha[j]==1.0)
      {}
      else if (alpha[j]<1e-8 || alpha[j]>1.-1e-8)
      {
        printf("alpha[%lu] = %3.16e\n",j,alpha[j]);
        // the point is not inside this element;
        inside = false;
        break;
      }

      // save the metric at this vertex
      this->Melem_[j] = tensors_[ topology_(elem,j) ];
    }
    printInline( alpha , "\talpha: ");

    if (!inside) continue;

    // perform the interpolation
    field_[p].template interpolate<LogEuclidean<real>>(alpha,this->Melem_);
    field_.cell(p) = elem;
		field_.sqrtDetM(p) = sqrt( field_[p].determinant() );

    break; // we found the element so we're done
  }
  luna_assert(inside); // we should have found something
}

template<typename type>
bool
MetricField<type>::recompute( index_t p , real* x )
{
	// look for the element in the background topology with the searcher
  index_t guess = field_.cell(p);
  int ielem = searcher_.find( x , guess );
  if (ielem<0)
  {
    // point is probably outside domain
		// let's make sure by first brute forcing the check
    ielem = searcher_.brute(x);
    if (ielem<0)
    {
			#if 1
			std::vector<real> alpha(this->n_+1,0.);
			ielem = searcher_.closest( x , alpha );
			index_t elem = index_t(ielem);
			for (index_t j=0;j<alpha.size();j++)
				this->Melem_[j] = tensors_[topology_(elem,j)];
			field_[p].template interpolate<LogEuclidean<real>>(alpha,this->Melem_);
			field_.cell(p) = elem;
			field_.sqrtDetM(p) = sqrt( field_[p].determinant() );
			return true;
			#endif

			// point is definitely outside the domain
			// there might be a projection to a geometry entity
			// and some tolerances we need to account for
      return false;
    }
  }

  index_t elem = index_t(ielem);
  const coord_t nv = topology_.nv(elem);

  // get the vertices of the element
  std::vector<real*> xk(nv,0);
  for (index_t j=0;j<nv;j++)
    xk[j] = topology_.vertices()[ topology_(elem,j) ];
  real V = 1./geometrics::signedVolume(xk,topology_.vertices().dim());

  // compute the barycentric coordinates
  std::vector<real> alpha( nv , 0. );
	bool problem = false;
  for (index_t j=0;j<nv;j++)
  {
    std::vector<real*> xk0(nv);
    for (index_t i=0;i<nv;i++)
    {
      if (i==j) xk0[i] = x;
      else xk0[i] = xk[i];
    }
    alpha[j] = V*geometrics::signedVolume(xk0,topology_.vertices().dim());

		//if (fabs(alpha[j]-1.0)>0.0 && fabs(alpha[j]-1.0)<1e-12)
		// 	alpha[j] = 1.0; // round-off due to above multiplication can exceed machine precision

    if (alpha[j]<0. || alpha[j]>1.)
    {
			printf("alpha[%lu] = %1.16e\n",j,alpha[j]);
      printf("element search failed\n");
      luna_assert_not_reached;
    }

		// signal a problem
		//if (alpha[j]<1e-12) problem = true;

    // save the metric at this vertex
    this->Melem_[j] = tensors_[ topology_(elem,j) ];
  }

	if (false and problem)
	{
		index_t nb_geometry = 0;
		for (index_t j=0;j<nv;j++)
			if (topology_.vertices().entity(j))
				nb_geometry++;
		if (nb_geometry==nv) return false;

		// determine how many nonzero barycentric coordinates there are
		std::vector<const real*> X;
		std::vector<SPDT<real>> minterp;
		for (index_t j=0;j<alpha.size();j++)
		{
			if (alpha[j]<1e-8) continue;
			X.push_back( xk[j] );
			minterp.push_back( this->Melem_[j] );
		}

		if (X.size()<alpha.size())
		{
			// compute the barycentric coordinates w.r.t this simplex facet
			std::vector<real> w(X.size());
			geometrics::barycentric( x , X , topology_.vertices().dim() , w );

			// the barycentric coordinate calculation is sensitive to numerical precision
			real wsum = 0.;
			for (coord_t j=0;j<w.size();j++)
				wsum += w[j];
			if (fabs(wsum-1.)>1e-20) return false; // signal an error

  		field_[p].template interpolate<LogEuclidean<real>>(w,minterp);
  		field_.cell(p) = elem;
			field_.sqrtDetM(p) = sqrt( field_[p].determinant() );
			return true;
		}

	} // problem

  // perform the interpolation
  field_[p].template interpolate<LogEuclidean<real>>(alpha,this->Melem_);
  field_.cell(p) = elem;
	field_.sqrtDetM(p) = sqrt( field_[p].determinant() );
	return true;
}

template<typename type>
bool
MetricField<type>::checkCells()
{
  index_t nb_violate = 0;
  for (index_t k=0;k<field_.vertices().nb();k++)
  {
    if (k<field_.vertices().nb_ghost()) continue;

    // get the element in the topology this vertex is inside
    index_t elem = field_.cell(k);
    real* x = field_.vertices()[k];
    luna_assert( !topology_.ghost(elem) );

    // ensure the barcyentric coordinates are ok
    bool inside = true;
    const coord_t nv = topology_.nv(elem);

    // get the vertices of the element
    std::vector<real*> xk(nv,0);
    for (index_t j=0;j<nv;j++)
      xk[j] = topology_.vertices()[ topology_(elem,j) ];
    real V = 1./geometrics::signedVolume(xk,topology_.vertices().dim());

    // compute the barycentric coordinates
    std::vector<real> alpha( nv , 0. );
    for (index_t j=0;j<nv;j++)
    {
      std::vector<real*> xk0(nv);
      for (index_t i=0;i<nv;i++)
      {
        if (i==j) xk0[i] = x;
        else xk0[i] = xk[i];
      }
      alpha[j] = V*geometrics::signedVolume(xk0,topology_.vertices().dim());

      if (alpha[j]<0.0 || alpha[j]>1.0)
      {
        inside = false;
      }
    }

    if (!inside)
    {
      nb_violate++;
      printf("vertex %lu is not in cell %lu!\n",k,elem);
      printInline(alpha,"barycentric coords: ");
    }
  }
  if (nb_violate>0) printf("there are %lu violations\n",nb_violate);
  return nb_violate==0;
}

template<typename type>
void
MetricField<type>::remove( index_t k )
{
  field_.remove(k);
}

void
edge_vector( const index_t i , const index_t j , const Vertices& v ,
              std::vector<real>& X )
{
  luna_assert( X.size()==v.dim() );
  for (coord_t d=0;d<v.dim();d++)
    X[d] = v[j][d] -v[i][d];
}

MetricAttachment::MetricAttachment( VertexField<SPDT<real>>& fld , Vertices& v ) :
  n_(fld[0].n()), vertices_(v)
{
  luna_assert_msg( fld.nb()==v.nb() , "fld.nb = %lu , v.nb = %lu" , fld.nb() , v.nb()  );
  cell_.resize(v.nb());
  for (index_t k=0;k<fld.nb();k++)
  {
    if (v.ghost(k))
    {
      SPDT<real> Mghost(n_);
      for (index_t i=0;i<n_;i++)
        Mghost(i,i) = 1.;
      VertexField<SPDT<real>>::add( Mghost );
			logM_.push_back( Mghost.log() );
			sqrtDetM_.push_back( 0. );
    }
    else
		{
      VertexField<SPDT<real>>::add( fld[k] );
			logM_.push_back( fld[k].log() );
			sqrtDetM_.push_back( sqrt( fld[k].determinant() ) );
		}
  }
}

MetricAttachment::MetricAttachment( Vertices& v ) :
	vertices_(v)
{
	cell_.resize( v.nb() );
	VertexField<SPDT<real>>::data_.resize( v.nb() );
}

template<typename type>
void
MetricAttachment::setCells( Topology<type>& topology )
{
  std::vector<index_t> vertex(1);
  int offset = topology.vertices().nb_ghost() -vertices_.nb_ghost();
  luna_assert( vertices_.nb()+offset==topology.vertices().nb() );

	#if 0
  for (index_t k=0;k<vertices_.nb();k++)
  {

    if (k<vertices_.nb_ghost())
    {
      cell_[k] = 0;
      continue;
    }

    vertex[0] = k +offset;
    std::vector<index_t> elems;
    topology.allWithSubset(vertex,elems);
    if (elems.size()==0)
    {
      printf("vertex %lu touches no elements?\n",k);
      vertices_.print(k,true);
      topology.printData();
    }
    luna_assert( elems.size()>0 );

    bool found = false;
    for (index_t j=0;j<elems.size();j++)
    {
      if (topology.ghost(elems[j])) continue;
      cell_[k] = elems[j];
      found = true;
      break;
    }
    if (!found)
    {
      vertices_.print("v",true);
      topology.printData();
      cell_[k] = topology.nb();
    }
    //luna_assert_msg(found,"vertex %lu not found in topology, offset = %d",k,offset);
  }
	#else
	std::vector<bool> visited( vertices_.nb() , false );
	index_t counted = vertices_.nb_ghost();
	for (index_t k=0;k<vertices_.nb_ghost();k++)
		cell_[k] = 0; // doesn't matter, won't be used for interpolation

	for (index_t k=0;k<topology.nb();k++)
	{
		if (topology.ghost(k)) continue;
		for (index_t j=0;j<topology.nv(k);j++)
		{
			if (visited[topology(k,j)]) continue;

			counted++;
			cell_[ topology(k,j) ] = k;
			visited[topology(k,j)] = true;
		}
	}
	luna_assert_msg(counted==vertices_.nb(),
									"counted = %lu, nb_vertices = %lu",counted,vertices_.nb());
	#endif
}

template<typename type>
void
MetricAttachment::limit( const Topology<type>& topology , real href )
{
	const coord_t dim = topology.number();
	luna_assert_msg( topology.vertices().dim() == dim , "dim = %u , num = %u" , topology.vertices().dim() , dim );

	// the vertices should be associated with each other
	luna_assert_msg( topology.vertices().nb() == vertices_.nb() , "topology nb_vertices = %lu, vertices_.nb() = %lu" ,
										topology.vertices().nb() , vertices_.nb() );

	// compute the implied metric of the input topology
	MeshImpliedMetric<type> implied( topology );
	implied.initialize();
	implied.optimize();

	// go through the entries of the current field and limit them using the step
	index_t nb_limited = 0;
	for (index_t k=0;k<nb();k++)
	{
		if (k < topology.vertices().nb_ghost()) continue;

		// compute the step from the implied metric to the current metric
		SPDT<real> mi = implied[k];
		SPDT<real> mt = this->operator[](k);

		SPDT<real> invsqrt_M0 = mi.pow(-0.5);
		SPDT<real> expS = mt.sandwich(invsqrt_M0);
		SPDT<real> s = expS.log();

		// limit the step
		bool limited = false;
		for (index_t i=0;i<dim;i++)
		for (index_t j=0;j<=i;j++)
		{
			if (s(i,j) > 2*log(href))
			{
				s(i,j) = 2*log(href);
				limited = true;
			}
			else if (s(i,j) < -2*log(href))
			{
				s(i,j) = -2*log(href);
				limited = true;
			}
		}
		if (limited) nb_limited++;
		else continue;

		// compute the new metric from the step
		SPDT<real> sqrt_M0 = mi.sqrt();
		SPDT<real> mk = (s.exp()).sandwich(sqrt_M0);
		this->operator[](k) = mk;

		logM_[k]     = mk.log();
		sqrtDetM_[k] = sqrt( mk.determinant() );

		// check positive-definiteness
		std::pair< std::vector<real> , densMat<real> > eig = mk.eig();
		for (index_t j=0;j<eig.first.size();j++)
			luna_assert_msg( eig.first[j] > 0 , "lambda(%lu) = %g" , j , eig.first[j] );
	}
	printf("limited %lu metrics out of %lu\n",nb_limited,this->nb());
}

void
MetricAttachment::reset( MetricAttachment& fld )
{
  cell_.resize( fld.nb() );
  VertexField<SPDT<real>>::clear();
  for (index_t k=0;k<fld.nb();k++)
    VertexField<SPDT<real>>::add(fld[k]);
}

void
MetricAttachment::setField( VertexField<SPDT<real>>& fld )
{
  VertexField<SPDT<real>>::clear();
  for (index_t k=0;k<fld.nb();k++)
    VertexField<SPDT<real>>::add( fld[k] );
}

index_t
MetricAttachment::rank() const
{
  return this->n_*(this->n_+1)/2;
}

std::vector<std::string>
MetricAttachment::ranknames() const
{
  std::vector<std::string> result;
  for (index_t i=0;i<this->n_;i++)
  for (index_t j=i;j<this->n_;j++)
    result.push_back("m("+stringify(i)+","+stringify(j)+")");
  return result;
}

std::vector<real>
MetricAttachment::lims() const
{
  std::vector<real> l;// = {0,1e3};
  return l;
}

void
MetricAttachment::add( SPDT<real>& T , index_t elem )
{
  Field<SPDT<real>>::add(T);
  cell_.push_back(elem);
	logM_.push_back( T.log() );
	sqrtDetM_.push_back( sqrt(T.determinant()) );
}

void
MetricAttachment::assign( index_t p , const numerics::SPDT<real>& M0 , index_t elem0 )
{
	for (index_t j=0;j<M0.nb();j++)
		Field<SPDT<real>>::operator[](p).data(j) = M0.data(j);
	cell_[p] = elem0;
	logM_[p] = M0.log();
	sqrtDetM_[p] = sqrt(M0.determinant());
}

void
MetricAttachment::remove( index_t k , bool recheck )
{
  Field<SPDT<real>>::remove(k);
  cell_.erase( cell_.begin() +k );
	logM_.erase( logM_.begin() +k );
	sqrtDetM_.erase( sqrtDetM_.begin() +k );
	if (recheck)
	{
	  luna_assert_msg( check() ,
	  "nb_vertices = %lu, nb_cell = %lu, Field<SPDT<real>>nb = %lu" ,
	  vertices_.nb(),cell_.size(),Field<SPDT<real>>::nb() );
	}
}

bool
MetricAttachment::check() const
{
  if (vertices_.nb()!=cell_.size()) return false;
  if (vertices_.nb()!=Field<SPDT<real>>::nb()) return false;
	if (vertices_.nb()!=logM_.size()) return false;
	if (vertices_.nb()!=sqrtDetM_.size()) return false;
  return true;
}

void
MetricAttachment::toJSON( json& J ) const
{
  J["name"] = "metric";
  J["membertype"] = "SPDT";
  J["evaltype"] = "vertex";
  J["numbertype"] = "real";
  J["nb"] = this->nb();
  index_t nb_rank = n_*(n_+1)/2;
  J["nb_rank"] = nb_rank;
  J["dim"] = n_;
  std::vector<real> data(this->nb()*nb_rank);
  index_t i = 0;
  for (index_t k=0;k<this->nb();k++)
  for (index_t j=0;j<nb_rank;j++)
    data[i++] = (*this)[k].data(j);
  J["data"] = data;
}

#if 0
void
MetricAttachment::to_solb( const std::string& filename ) const
{
	int64_t fid;

	double buf[GmfMaxTyp];
	int dim = vertices_.dim();
	fid = GmfOpenMesh(filename.c_str(),GmfWrite,GmfDouble,dim);
	luna_assert( fid );

	int TypTab[GmfMaxTyp];
	TypTab[0] = GmfSymMat;
	GmfSetKwd( fid , GmfSolAtVertices , this->nb() , 1 , TypTab );

	for (index_t k=0;k<this->nb();k++)
	{
		const numerics::SPDT<real>& m = this->operator[](k);
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
	real dvalues[GmfMaxTyp];
	int version; // 2 for 32-bit int, 64-bit real
	int64_t fid;

	fid = GmfOpenMesh(filename.c_str(),GmfRead,&version,&dim);
	luna_assert_msg( fid , "could not open sol file %s ",filename.c_str() );

	printf("version = %d, dimension = %d\n",version,dim);

	// create a field whether this is attached at vertices or cells

	nb_sol = GmfStatKwd( fid , GmfSolAtVertices , &numberType , &solSize , TypTab );
	printf("nb_sol = %d, numberType = %d, solSize = %d\n",nb_sol,numberType,solSize);

	luna_assert( nb_sol == int(this->vertices_.nb()) );
	luna_assert( solSize == int(dim*(dim+1)/2) );

	luna_assert( GmfGotoKwd( fid , GmfSolAtVertices ) > 0 );
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

		luna_assert( status==1 );

		std::vector<real> data(dvalues,dvalues+solSize);
		numerics::SPDT<real> m(data);
		this->operator[](k) = m;
	}
}
#endif

// discrete metric and embedding instantiations
template class MetricField<Simplex>;
template bool MetricField<Simplex>::check( Topology<Simplex>& );
template void MetricAttachment::setCells( Topology<Simplex>& );
template void MetricAttachment::limit( const Topology<Simplex>& , real );

#endif

} // luna
