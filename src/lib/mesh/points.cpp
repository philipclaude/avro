#include "common/tools.h"

#include "geometry/body.h"
#include "geometry/entity.h"
#include "geometry/model.h"

#include "mesh/topology.h"
#include "mesh/points.h"

#include "numerics/geometry.h"

#include <egads.h>
#include <json/json.hpp>

#include <math.h>
#include <unordered_map>

namespace luna
{

Points::Points() :
	DOF<real_t>(0),
	dim_(0),
	udim_(0),
	u_(0),
	nb_ghost_(0)
{}

Points::Points( const coord_t _dim ) :
	DOF<real_t>(_dim),
	dim_(_dim),
	udim_(dim_-1), // default to assuming parameter space is dim-1
	u_(udim_),
	nb_ghost_(0)
{}

Points::Points( const coord_t _dim , const coord_t _udim ) :
	DOF<real_t>(_dim),
	dim_(_dim),
	udim_(_udim),
	u_(udim_),
	nb_ghost_(0)
{}

Points::~Points()
{
	clear();
}

void
Points::create( const std::vector<real_t>& x )
{
	luna_assert(x.size()==dim_);
	create(x.data());
}

void
Points::create( const real_t* x )
{
	DOF<real_t>::add(x,dim_);

	std::vector<real_t> infty(udim_,INFTY);
	u_.add( infty.data() , udim_ );

	body_.add(0);
	primitive_.add(NULL);
	fixed_.add(false);
}

void
Points::copy( Points& v , const bool ghosts) const
{
  v.set_dim( dim_ );
	v.set_parameter_dim( udim_ );

  // copy the coordinates
  for (index_t k=nb_ghost();k<nb();k++)
    v.create( operator[](k) );

  // copy the ghost vertices
  if (ghosts)
  {
    for (index_t k=0;k<nb_ghost();k++)
      v.create_ghost();
  }

  // copy the meta data
	// TODO account for ghost offset...
  for (index_t k=0;k<nb();k++)
  {
    v.body(k)   = body_[k];
    //v.setEntity( k , primitive_[k] );
    v.set_fixed( k , fixed_[k] );
    //v.setParam( k , u(k) );
  }
}

void
Points::create_ghost()
{
	// add some arbitrary ghost vertices
	std::vector<real_t> x0( dim_ , INFTY );

	DOF<real_t>::insert( nb_ghost_ , x0.data() , x0.size() );

	std::vector<real_t> u0( udim_ , INFTY );

	u_.insert( nb_ghost_ , u0.data() , u0.size() );
	body_.insert( nb_ghost_ , nb_ghost_ );
	primitive_.insert( nb_ghost_ , NULL );
	fixed_.insert( nb_ghost_ , false );

	// increment the ghost counter
	nb_ghost_++;
}

int&
Points::body( const index_t k )
{
	luna_assert_msg( k<nb() , "k = %lu , nb = %lu" , k , nb() );
	return body_[k];
}

#if 1

void
Points::set_entity( const index_t k , Entity* e )
{
	luna_assert( k<nb() );
	primitive_[k] = e;
}

void
Points::set_param( const index_t k , const std::vector<real_t>& u )
{
	luna_assert(u.size()==udim_);
	set_param( k , u.data() );
}

void
Points::set_param( const index_t k , const real_t* u )
{
  luna_assert( k<nb() );
	u_.set( k , u );
}

void
Points::print( bool info ) const
{
	printf("Points (%p):\n",(void*)this);
	for (index_t k=0;k<nb();k++)
	{
		printf("p[%4d]: (",int(k));
		for (index_t d=0;d<dim_;d++)
			printf(" %12.4e ", operator[](k)[d] );
		printf(")");
		if (info)
		{
      std::string geo;
			int num = -1;
  		if (primitive_[k]!=NULL)
  		{
  			geo = "-"+primitive_[k]->name();
				num = primitive_[k]->number();
  		}
  		else geo = "";
  		printf(" : %s , b[ %3d ] , g[ %1d-%p%s ] , u = (",(k<nb_ghost_)? "GHST":"REAL",
  						body_[k],num,(void*)primitive_[k],geo.c_str());
			for (index_t d=0;d<udim_;d++)
				printf(" %12.4e ",u(k)[d]);
			printf(")");
		}
		printf("\n");
	}
}

void
Points::print( index_t k , bool info ) const
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
		printf(" : %s , b[ %3d ] , g[ %p%s ] , u = (",(k<nb_ghost_)? "GHST":"REAL",
						body_[k],(void*)primitive_[k],geo.c_str());
		for (index_t d=0;d<udim_;d++)
			printf("%12.4e ",u(k)[d]);
		printf(")");
	}
	printf("\n");
}

void
Points::remove( const index_t k )
{
	luna_assert_msg( k < nb() ,
		"k = %lu , nb = %lu" , k , nb() );

	DOF<real_t>::remove(k);

	u_.remove(k);
	body_.remove(k);
	primitive_.remove(k);
	fixed_.remove(k);

	if (k<nb_ghost_) nb_ghost_--;
}

void
Points::duplicates( std::vector<index_t>& idx ,real_t tol ) const
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
Points::duplicates( std::vector<index_t>& idx , const Table<int>& F ) const
{
	luna_assert( F.nb() == nb() );

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
			luna_assert( d < 1e-6 ); // pretty loose tolerance
		}
	}
	printf("there are %lu unique vertices\n",symbolic.size());
	luna_implement;
}

void
Points::dump( const std::string& filename ) const
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

void
Points::attach( const Body& body , index_t ibody ,real_t tol )
{
  // useful for when a mesh is read but it does not have a geometry
  // even though we know what it should be
  std::vector<real_t> x(dim_,0.);

  std::vector<Entity*> e_candidates;
  std::vector<std::vector<real_t>> u_candidates;

  bool recheck;

  // get the full list of primitives
  std::vector<Entity*> primitives;
  body.get_tessellatable(primitives);

  std::vector<real_t> distances( primitives.size() , 0. );
  std::vector<real_t> u( udim_*primitives.size() );

  if (dim_ != 4) luna_assert_msg(udim_ > 0, "udim = %u" , udim_ );
  std::vector<real_t> uk(udim_, INFTY);

  for (index_t k=0;k<nb();k++)
  {
    recheck = true;

    const real_t* xk = (*this)[k];

    e_candidates.clear();
    u_candidates.clear();

    // compute the distance to each primitive
    for (index_t j=0;j<primitives.size();j++)
    {
      // re-assign the original coordinates
			x.assign( xk , xk+dim_ );

      Entity* e = primitives[j];

      // project to the primitive and compute the distance
      primitives[j]->inverse( x , uk );

			// compute the distance between the inverse evaluated point and the original one
      distances[j] = numerics::distance2( xk , x.data() , dim_ );

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
      luna_assert(primitive_[k] == NULL);
			//printf("no candidates for vertex %lu\n",k);
      continue;
    }

    // find the candidate with the lowest topological number
    Entity* ek = e_candidates[0];
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
      luna_assert(ek == primitive_[k]);
    }
    else
      set_entity(k,ek);

    set_param(k,uk);
    body_[k] = ibody;

		ek->evaluate( uk , x );

    // get the primitive coordinates
		if (recheck)
		{
    	// check the distance is lower than the tolerance
    	real_t d = numerics::distance2(xk,x.data(),dim_);
    	luna_assert_msg( d < tol , "point %lu , d = %1.16e" , k , d );
		}
  }
}

void
Points::attach( const Model& model , real_t tol )
{
  for (index_t k=0;k<model.nb_bodies();k++)
    attach( model.body(k) , k+1 , tol );
}

/*
void
Points::computePartition( ArrayBase<index_t>& data ,
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
Points::to_json( json& J ) const
{
  J["dimension"] = dim_;
  J["nb_ghost"] = nb_ghost_;
  J["coordinates"] = DOF<real_t>::data();
  J["parameters"] = u_.data();

  std::vector<int> geometry(nb(),-1);
  for (index_t k=0;k<nb();k++)
  {
    if (primitive_[k]==NULL) continue;
    geometry[k] = primitive_[k]->identifier();
  }
  J["geometry"] = geometry;
}

void
Points::from_json( const json& J , const Model* model )
{
  dim_ = J["dimension"];
  std::vector<real_t> X = J.at("coordinates");
  DOF<real_t>::data_ = X;
  //std::vector<real_t> U = J.at("parameters");
  //u_ = U;
  nb_ghost_ = J["nb_ghost"];

  primitive_.resize( nb() , NULL );
  body_.resize( nb() , 0 );

  if (model!=NULL)
  {
    std::vector<Entity*> primitives;
    model->get_entities(primitives);
    std::vector<int> geometry = J.at("geometry");
    luna_assert(geometry.size()==nb());
    for (index_t k=0;k<nb();k++)
    {
      if (geometry[k]<0) continue;
      for (index_t j=0;j<primitives.size();j++)
      {
        if ((int)primitives[j]->identifier()==geometry[j])
        {
          set_entity(k,primitives[j]);
          break;
        }
      }
    }
  }
}

#endif

void
Points::compute_param( index_t k )
{
	if (primitive_[k]==NULL) return;

	std::vector<real_t> x( (*this)[k] , (*this)[k] + dim_ );
	std::vector<real_t> U( u(k) , u(k) + udim_ );
	primitive_[k]->inverse(x,U);

	for (coord_t d=0;d<udim_;d++)
		u(k,d) = U[d];
}

void
Points::compute_params()
{
	for (index_t k=0;k<nb();k++)
		compute_param(k);
}

void
Points::clear()
{
  u_.clear();
  body_.clear();
  primitive_.clear();
  nb_ghost_ = 0;
}

} // luna
