//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
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

namespace avro
{

Points::Points() :
	DOF<real_t>(0),
	dim_(0),
	udim_(0),
	u_(0),
	incidence_(TableLayout_Jagged),
	nb_ghost_(0),
	parameter_space_(false)
{}

Points::Points( const coord_t _dim ) :
	DOF<real_t>(_dim),
	dim_(_dim),
	udim_(dim_-1), // default to assuming parameter space is dim-1
	u_(udim_),
	incidence_(TableLayout_Jagged),
	nb_ghost_(0),
	parameter_space_(false)
{}

Points::Points( const coord_t _dim , const coord_t _udim ) :
	DOF<real_t>(_dim),
	dim_(_dim),
	udim_(_udim),
	u_(udim_),
	incidence_(TableLayout_Jagged),
	nb_ghost_(0),
	parameter_space_(false)
{}

Points::~Points()
{
	clear();
}

void
Points::create( const std::vector<real_t>& x )
{
	avro_assert(x.size()==dim_);
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
	global_.add(nb()); // 1-bias
}

void
Points::copy( Points& v , const bool ghosts) const
{
	v.clear();

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
    v.set_entity( k , primitive_[k] );
    v.set_fixed( k , fixed_[k] );
    v.set_param( k , u(k) );
		v.set_global( k , global_[k] );
  }

	for (index_t k=0;k<incidence_.nb();k++)
	{
		std::vector<int> f = incidence_.get(k);
		v.incidence().add( f.data() , f.size() );
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
	global_.insert( nb_ghost_ , 0 );

	// increment the ghost counter
	nb_ghost_++;
}

int&
Points::body( const index_t k )
{
	avro_assert_msg( k<nb() , "k = %lu , nb = %lu" , k , nb() );
	return body_[k];
}

void
Points::set_entity( const index_t k , Entity* e )
{
	avro_assert( k<nb() );
	primitive_[k] = e;
}

void
Points::set_param( const index_t k , const std::vector<real_t>& u )
{
	avro_assert(u.size()==udim_);
	set_param( k , u.data() );
}

void
Points::set_param( const index_t k , const real_t* u )
{
  avro_assert( k<nb() );
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
  		printf(" : %5s , %5s , g[ %1d-%p%s ] , u = (",(k<nb_ghost_)? "ghost":"real",
  						(fixed_[k])?("fixed"):("free"),num,(void*)primitive_[k],geo.c_str());
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
		int num = -1;
		if (primitive_[k]!=NULL)
		{
			geo = "-"+primitive_[k]->name();
			num = primitive_[k]->number();
		}
		else geo = "";
		printf(" : %5s , %5s , g[ %1d-%p%s ] , u = (",(k<nb_ghost_)? "ghost":"real",
							(fixed_[k])?("fixed"):("free"),num,(void*)primitive_[k],geo.c_str());
		for (index_t d=0;d<udim_;d++)
			printf("%12.4e ",u(k)[d]);
		printf(")");
	}
	printf("\n");
}

void
Points::remove( const index_t k )
{
	avro_assert_msg( k < nb() ,
		"k = %lu , nb = %lu" , k , nb() );

	DOF<real_t>::remove(k);

	u_.remove(k);
	body_.remove(k);
	primitive_.remove(k);
	fixed_.remove(k);
	global_.remove(k);

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
		for (index_t j=k+1;j<nb();j++)
		{
			if (numerics::distance(operator[](k),operator[](j),dim_)<tol)
			{
				idx[j] = k;
				break;
			}
		}
	}
}

void
Points::duplicates( std::vector<index_t>& idx , const Table<int>& F ) const
{
	avro_assert( F.nb() == nb() );
	idx.resize( nb() );

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
			avro_assert( d < 1e-6 ); // pretty loose tolerance
		}
	}
	printf("there are %lu unique vertices\n",symbolic.size());
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

  if (dim_ != 4) avro_assert_msg(udim_ > 0, "udim = %u" , udim_ );
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
      avro_assert(primitive_[k] == NULL);
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
      avro_assert(ek == primitive_[k]);
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
			if (d>=tol)
			{
				print();
				ek->print();
			}
    	avro_assert_msg( d < tol , "point %lu , d = %1.16e" , k , d );
		}
  }
}

void
Points::attach( const Model& model , real_t tol )
{
  for (index_t k=0;k<model.nb_bodies();k++)
    attach( model.body(k) , k+1 , tol );
}

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
    avro_assert(geometry.size()==nb());
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
Points::extract_params( index_t k , Entity* e0 , real_t* param ) const
{

	const real_t* uk = (parameter_space_) ? (*this)[k] : u(k);

	if (e0 == entity(k))
	{
		for (coord_t d=0;d<udim_;d++)
			param[d] = uk[d];
		return;
	}

#if 0
	EGADS::Object* e = (EGADS::Object*) e0;
	EGADS::Object* g = (EGADS::Object*) entity(k);
	if (e->number()==2)
	{
		if (g->number()==1)
			e->extract_params_edge_on_face( uk , g , param );
		else if (g->number()==0)
			e->extract_params_node_on_face( uk , g , param );
		else
			avro_assert_not_reached;
	}
	else if (e->number()==1)
	{
		// extract_params_face_on_edge()
		if (g->number()==0)
			e->extract_params_node_on_edge( uk , g , param );
		else
			avro_assert_not_reached;
	}

	// check the coordinates are close if we have physical coordinates
	if (!parameter_space_)
	{
		/*
		real_t d = 0.0;
		real_t tol = min of tolerance on e and g
		e->evaluate();
		d = numerics::distance( )
		avro_assert( d < tol );
		*/
	}
	#endif
}

void
Points::move_to( index_t k0 , index_t k1 )
{
	std::vector<real_t> x0( (*this)[k0] , (*this)[k0]+dim_ );
	std::vector<real_t> u0( u(k0) , u(k0)+udim_ );
	Entity* e0 = entity(k0);

	avro_assert_msg( k1 <= k0 , "trying to move %lu to %lu"  , k0 , k1 );
	#if 0
	DOF<real_t>::insert( k1*dim_ , (*this)[k0] , dim_ );
	u_.insert( k1*udim_ , u(k0) , udim_ );
	body_.insert( k1 , body(k0) );
	primitive_.insert( k1 , entity(k0) );
	fixed_.insert( k1 , fixed(k0) );
	#else
	DOF<real_t>::insert( k1*dim_ , x0.data() , x0.size() );
	u_.insert( k1*udim_ , u0.data() , u0.size() );
	body_.insert( k1 , body(k0) );
	primitive_.insert( k1 , e0 );
	fixed_.insert( k1 , fixed(k0) );
	#endif
	remove(k0+1); // +1 because we added a vertex in front of k0
}

void
Points::clear()
{
	DOF<real_t>::clear();
  u_.clear();
  body_.clear();
  primitive_.clear();
  nb_ghost_ = 0;
	incidence_.clear();
	fixed_.clear();
}

} // avro
