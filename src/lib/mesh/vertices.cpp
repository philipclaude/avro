#include "common/tools.h"

//#include "geometry/body.h"
#include "geometrics/primitive.h"
//#include "geometry/model.h"

#include "mesh/topology.h"
#include "mesh/vertices.h"

#include "numerics/geometry.h"

#include <egads.h>
#include <json/json.hpp>

#include <math.h>
#include <unordered_map>

namespace ursa
{

Vertices::Vertices( const coord_t _dim ) :
	dim_(_dim),
	udim_(dim_-1), // default to assuming parameter space is dim-1
	ghost_(0)
{
	x_.clear();
	u_.clear();
	body_.clear();
	primitive_.clear();
	fixed_.clear();
}

Vertices::Vertices( const coord_t _dim , const coord_t _udim ) :
	dim_(_dim),
	udim_(_udim),
	ghost_(0)
{
	x_.clear();
	u_.clear();
	body_.clear();
	primitive_.clear();
	fixed_.clear();
}

Vertices::~Vertices()
{
	clear();
}

void
Vertices::create( const std::vector<real_t>& x )
{
	ursa_assert(x.size()==dim_);
	create(x.data());
}

void
Vertices::create( const real_t* x )
{
  for (index_t i=0;i<dim_;i++)
    x_.push_back( x[i] );
  for (index_t i=0;i<udim_;i++)
    u_.push_back(INFTY);
  body_.push_back(0);
  primitive_.push_back(NULL);
  const bool f = false;
  fixed_.push_back(f);
}

void
Vertices::copy( Vertices& v , const bool erase , const bool ghosts) const
{
  if (erase)
    v.clear();

  v.setDimension( dim_ );
	v.setParameterDimension( udim_ );

  // copy the coordinates
  for (index_t k=nb_ghost();k<nb();k++)
    v.create( operator[](k) );

  // copy the ghost vertices
  if (ghosts)
  {
    for (index_t k=0;k<nb_ghost();k++)
      v.createGhost();
  }

  // copy the meta data
	// TODO account for ghost offset...
  for (index_t k=0;k<nb();k++)
  {
    v.body(k)   = body_[k];
    //v.setPrimitive( k , primitive_[k] );
    v.setFixed( k , fixed_[k] );
    //v.setParam( k , u(k) );
  }
}

void
Vertices::createGhost()
{
	// add some arbitrary ghost vertices
	std::vector<real_t> x0( dim_ , INFTY );
	x_.insert( x_.begin()+ghost_ , x0.begin() , x0.end() );

	std::vector<real_t> u0( udim_ , INFTY );
	u_.insert( u_.begin()+ghost_ , u0.begin() , u0.end() );

	body_.insert( body_.begin() +ghost_ , ghost_ );
	primitive_.insert( primitive_.begin() +ghost_ , NULL );
	fixed_.insert( fixed_.begin() + ghost_ , false );

	// increment the ghost counter
	ghost_++;
}

int&
Vertices::body( const index_t k )
{
	ursa_assert_msg( k<nb() , "k = %lu , nb = %lu" , k , nb() );
	return body_[k];
}

#if 1

void
Vertices::setPrimitive( const index_t k , geometrics::Primitive* e )
{
	ursa_assert( k<nb() );
	primitive_[k] = e;
}

void
Vertices::setParam( const index_t k , const std::vector<real_t>& u )
{
	ursa_assert(u.size()==udim_);
	setParam( k , u.data() );
}

void
Vertices::setParam( const index_t k , const real_t* u )
{
  ursa_assert( k<nb() );
  for (index_t i=0;i<udim_;i++)
    u_[udim_*k+i] = u[i];
}

bool
Vertices::boundary( const index_t k ) const
{
	ursa_assert_msg( k<nb() , "k = %lu , nb = %lu" , k , nb() );
	return body_[k]>0;
}

void
Vertices::print( std::string pre , bool info ) const
{
	printf("Vertices (%p):\n",(void*)this);
	if (pre=="\0") pre = "v";
	for (index_t k=0;k<nb();k++)
	{
		printf("%s[%4d]: (",pre.c_str(),int(k));
		for (index_t d=0;d<dim_;d++)
			printf(" %12.4e ", operator[](k)[d] );
		printf(")");
		if (info)
		{
      std::string geo;
  		if (primitive_[k]!=NULL)
  		{
  			geo = "-"+primitive_[k]->name();
  		}
  		else geo = "";
  		printf(" : %s , b[ %3d ] , g[ %p%s ] , u = (",(k<ghost_)? "GHST":"REAL",
  						body_[k],(void*)primitive_[k],geo.c_str());
			for (index_t d=0;d<udim_;d++)
				printf(" %12.4e ",u(k)[d]);
			printf(")");
		}
		printf("\n");
	}
}

void
Vertices::print( const index_t k , bool info ) const
{
	printf("vertex[%4d]: (",int(k));
	for (index_t d=0;d<dim_;d++)
		printf(" %12.4e ", operator[](k)[d] );
	printf(")");
	if (info)
	{
		std::string geo;
		if (primitive_[k]!=NULL)
		{
			geo = "-"+primitive_[k]->name();
		}
		else geo = "";
		printf(" : %s , b[ %3d ] , g[ %p%s ] , u = (",(k<ghost_)? "GHST":"REAL",
						body_[k],(void*)primitive_[k],geo.c_str());
		for (index_t d=0;d<udim_;d++)
			printf("%12.4e ",u(k)[d]);
		printf(")");
	}
	printf("\n");
}

void
Vertices::remove( const index_t k )
{
	ursa_assert_msg( k < nb() ,
		"k = %lu , nb = %lu , |x| = %lu" , k , nb() , x_.size() );
	x_.erase( x_.begin()+dim_*k , x_.begin()+dim_*(k+1) );
	u_.erase( u_.begin()+udim_*k , u_.begin()+udim_*(k+1) );
  if (k<ghost_) ghost_--;
	body_.erase( body_.begin() +k );
	primitive_.erase( primitive_.begin() + k );
	fixed_.erase( fixed_.begin() + k );
}

void
Vertices::duplicates( std::vector<index_t>& idx ,real_t tol ) const
{
	// initialize the map
	idx.resize(nb());
	for (index_t k=0;k<nb();k++)
		idx[k] = k;

	// TODO use parallel scheduler #pragma omp parallel for
	for (index_t k=0;k<nb();k++)
	{
		// look for any vertices up to this one which are too close
		for (index_t j=0;j<k;j++)
		{
			if (numerics::distance(operator[](k),operator[](j),dim_)<tol)
			{
				idx[k] = j;
				break;
			}
		}
	}
}

void
Vertices::duplicates( std::vector<index_t>& idx , const Data<int>& F ) const
{
	ursa_assert( F.nb() == nb() );

	// get the symbolic information of every vertex
	std::unordered_map<std::string,index_t> symbolic;

	// loop through the vertices
	for (index_t k=0;k<nb();k++)
	{
		std::vector<int> B = F.get(k);
		std::string s = unique_label(B);

		std::unordered_map<std::string,index_t>::const_iterator it;
		it = symbolic.find(s);
		if (it==symbolic.end())
		{
			// the symbolic vertex does not exist so add it
			symbolic.insert( {s,k} );
			idx[k] = k;
		}
		else
		{
			// the symbolic information exists
			idx[k] = it->second;
			printf("symbolic vertex exists! %lu -> %lu\n",k,it->second);

			real_t d = numerics::distance( (*this)[k] , (*this)[it->second] , dim_ );
			ursa_assert( d < 1e-6 ); // pretty loose tolerance
		}
	}
	printf("there are %lu unique vertices\n",symbolic.size());
	ursa_implement;
}

void
Vertices::dump( const std::string& filename ) const
{
  FILE* fid = fopen(filename.c_str(),"w");
  for (index_t k=0;k<nb();k++)
  {
    for (coord_t d=0;d<dim_;d++)
      fprintf(fid,"%.12e ",operator[](k)[d]);
    fprintf(fid,"\n");
  }
  fclose(fid);
}

/*
void
Vertices::findGeometry( const Body& body , index_t ibody ,real_t tol )
{
  // useful for when a mesh is read but it does not have a geometry
  // even though we know what it should be
  std::vector<real_t> x(4,0.);
  coord_t dlim = (dim_<=3) ? dim_ : 4;

  Primitive *e;
  std::vector<Primitive*> e_candidates;
  std::vector<std::vector<real_t>> u_candidates;

  bool recheck;

  // get the full list of primitives
  std::vector<Primitive*> primitives;
  body.listTessellatableEntities(primitives);
  std::vector<real_t> distances( primitives.size() , 0. );
  std::vector<real_t> u( udim_*primitives.size() );

  if (dim_ != 4) ursa_assert_msg(udim_ > 0, "udim = %u" , udim_ );
  std::vector<real_t> uk(udim_, INFTY);

  for (index_t k=0;k<nb();k++)
  {

    recheck = true;

    const real_t* xk = this->operator[](k);

    e_candidates.clear();
    u_candidates.clear();

    // compute the distance to each primitive
    for (index_t j=0;j<primitives.size();j++)
    {
      // re-assign the original coordinates
      for (coord_t d=0;d<dlim;d++)
        x[d] = xk[d];

      e = primitives[j];

      // project to the primitive and compute the distance
      primitives[j]->projectToGeometry( x.data() , uk.data() );

      distances[j] = numerics::distance2( xk , x.data() , dim_ );

      #if 0
      if (primitives[j]->number()==2 && !primitives[j]->tesseractGeometry())
      {
        // better to use inTopology for surface fit
        int icode = EG_inTopology( primitives[j]->object() , xk );
        if (icode==EGADS_SUCCESS)
        {
          recheck = false;
          distances[j] = 0.;
        }
      }
      #endif

      // add as candidate if the tolerance is satisfied
      if (distances[j]<tol)
      {
        e_candidates.push_back(e);
        u_candidates.push_back(uk);
      }
    }

    // no candidates, stay null (interior)
    if (e_candidates.size()==0)
    {
      ursa_assert(primitive_[k] == NULL);
			printf("no candidiates for vertex %lu\n",k);
      continue;
    }

    // find the candidate with the lowest topological number
    Primitive* ek = e_candidates[0];
    uk = u_candidates[0];
    for (index_t j=1;j<e_candidates.size();j++)
    {
      if (e_candidates[j]->number()<ek->number())
      {
        ek = e_candidates[j];
        uk = u_candidates[j];
      }
    }
    if (primitive_[k] != NULL)
    {
      ursa_assert(ek == primitive_[k]);
    }
    else
      setPrimitive(k,ek);

    setParam(k,uk);
    body_[k] = ibody;

    // get the primitive coordinates
    if (ek->tesseractGeometry())
    {
      // re-assign the original coordinates
      for (coord_t d=0;d<dlim;d++)
        x[d] = xk[d];

      ek->projectToGeometry(x.data());
    }
    else
      ek->evaluate( uk.data() , x.data() );

    // check the distance is lower than the tolerance
   real_t d = numerics::distance2(xk,x.data(),dim_);
    if (recheck) ursa_assert_msg( d < tol , "d = %1.16e" , d );
  }
}

void
Vertices::findGeometry( const Model& model ,real_t tol )
{
  for (index_t k=0;k<model.nb_bodies();k++)
    findGeometry( *model.body(k) , k+1 , tol );
}

void
Vertices::projectToGeometry( Body& body )
{
  // project the vertices to the egads geometry
  std::vector<real_t> x(4,0.);
  coord_t dlim = (dim_<=3) ? dim_ : 4;

  std::vector<Primitive*> e_candidates;
  std::vector<std::vector<real_t>> u_candidates;

  // get the full list of primitives
  std::vector<Primitive*> primitives;
  body.listTessellatableEntities(primitives);
  std::vector<real_t> distances( primitives.size() , 0. );
  std::vector<real_t> u( udim_*primitives.size() );

  if (dim_ != 4) ursa_assert(udim_ > 0);
  std::vector<real_t> uk(udim_);

  for (index_t k=0;k<nb();k++)
  {
    // find the closest primitive with lowest topological dimension
    const real_t* xk = this->operator[](k);

    // compute the distance to each primitive
    for (index_t j=0;j<primitives.size();j++)
    {
      // re-assign the original coordinates
      for (coord_t d=0;d<dlim;d++)
        x[d] = xk[d];

      // project to the primitive and compute the distance
      primitives[j]->projectToGeometry(x.data(), u.data() + udim_*j);
      distances[j] = numerics::distance2( xk , x.data() , dim_ );

    }

   real_t dmin = *std::min_element( distances.begin() , distances.end() );

    e_candidates.clear();
    u_candidates.clear();
    for (index_t j=0;j<distances.size();j++)
    {
      if (distances[j]<dmin+1e-12)
      {
        e_candidates.push_back( primitives[j] );
        u_candidates.push_back( std::vector<real_t>(u.begin() + udim_*j, u.begin() + udim_*(j+1)) );
      }
    }

    if (e_candidates.size()==0)
		{
			//printf("no candidates for vertex %lu\n",k);
			continue;
		}

    // find the candidate with the lowest topological number
    Primitive* ek = e_candidates[0];
    uk = u_candidates[0];
    for (index_t j=1;j<e_candidates.size();j++)
    {
      if (e_candidates[j]->number()<ek->number())
      {
        ek = e_candidates[j];
        uk = u_candidates[j];
      }
    }

    //print(k,true);
    setPrimitive(k,ek);
    setParam(k,uk);

    // recall the projection (this is only so we don't have to save all the coordinates)
   real_t xp[3];
    if (ek->tesseractGeometry())
    {
      for (coord_t d=0;d<dlim;d++)
        xp[d] = xk[d];
      ek->projectToGeometry( xp );
    }
    else
      ek->evaluate( uk.data() , xp );

    for (coord_t d=0;d<dlim;d++)
      this->operator[](k)[d] = xp[d];
  }

}

void
Vertices::computePartition( Data<index_t>& data ,
                      std::vector<index_t>& cell_partition , index_t nparts )
{
  // data is technically a topology, construct the adjacency graph
  graph::Graph adjacency(nb());

  std::vector<index_t> vertex(1);

  // vertices are connected if they share a cell
  for (index_t k=0;k<data.nb();k++)
  {

    // get all cells with this vertex
    std::vector<index_t> cells;
    vertex[0] = k;
    data.allWithSubset(vertex,cells);

    // get uniquee list of all vertices in touching cells
    std::vector<index_t> N;
    for (index_t j=0;j<cells.size();j++)
    for (index_t i=0;i<data.nv(cells[j]);i++)
      N.push_back( data(cells[j],i) );
    uniquify(N);

    // add an edge for each connected vertex
    for (index_t j=0;j<N.size();j++)
    {
      if (N[j]>=k) continue; // uniqueness of edges
      adjacency.addEdge( k , N[j] );
    }
  }

  // convert to csr format
  adjacency.tocsr();

  // compute the partition
  adjacency.partition( nparts , partition_ );

  // revisit the vertices and try to assign the cell partition
  for (index_t k=0;k<nb();k++)
  {
    // get all cells with this vertex
    std::vector<index_t> cells;
    vertex[0] = k;
    data.allWithSubset(vertex,cells);

    std::vector<index_t> pc(cells.size());
    for (index_t i=0;i<cells.size();i++)
      pc[i] = cell_partition[cells[i]];
    uniquify(pc);

    if (pc.size()==1)
      partition_[k] = pc[0];
    else
    {
      // pick the partition closest to this one
     real_t d = fabs(real_t(pc[0]-partition_[k]) );
      index_t p = 0;
      for (index_t j=1;j<pc.size();j++)
      {
       real_t dj = fabs(real_t(pc[j]-partition_[k]) );
        if (dj<d)
        {
          dj = d;
          p = j;
        }
      }
      partition_[k] = pc[p];
    }
  }

}
*/

void
Vertices::toJSON( json& J ) const
{
  J["dimension"] = dim_;
  J["nb_ghost"] = ghost_;
  J["coordinates"] = x_;
  J["parameters"] = u_;

  std::vector<int> geometry(nb(),-1);
  for (index_t k=0;k<nb();k++)
  {
    if (primitive_[k]==NULL) continue;
    geometry[k] = primitive_[k]->identifier();
  }
  J["geometry"] = geometry;
}

/*
void
Vertices::fromJSON( const json& J , const Model* model )
{
  dim_ = J["dimension"];
  std::vector<real_t> X = J.at("coordinates");
  x_ = X;
  std::vector<real_t> U = J.at("parameters");
  u_ = U;
  ghost_ = J["nb_ghost"];

  primitive_.resize( nb() , NULL );
  body_.resize( nb() );

  if (model!=NULL)
  {
    std::vector<geometrics::Primitive*> primitives;
    model->listEntities(primitives);
    std::vector<int> geometry = J.at("geometry");
    ursa_assert(geometry.size()==nb());
    for (index_t k=0;k<nb();k++)
    {
      if (geometry[k]<0) continue;
      for (index_t j=0;j<primitives.size();j++)
      {
        if ((int)primitives[j]->bodyIndex()==geometry[j])
        {
          setPrimitive(k,primitives[j]);
          break;
        }
      }
    }
  }
}
*/

#endif

void
Vertices::clear()
{
  x_.clear();
  u_.clear();
  body_.clear();
  primitive_.clear();
  ghost_ = 0;
}

} // ursa
