#include "common/array.h"
#include "common/tools.h"
#include "common/stringify.h"

//#include "geometry/body.h"
//#include "geometry/entity.h"
//#include "geometry/model.h"

//#include "mesh/geometrics.h"
//#include "mesh/mesh.h"
#include "mesh/topology.h"
#include "mesh/vertices.h"

#include <egads.h>
#include <json/json.hpp>

#include <math.h>

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
	mesh_.clear();
	entity_.clear();
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
	mesh_.clear();
	entity_.clear();
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
  mesh_.push_back(NULL);
  entity_.push_back(NULL);
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
    //v.setEntity( k , entity_[k] );
    //v.setMesh( k , mesh_[k] );
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
	mesh_.insert( mesh_.begin() +ghost_ , NULL );
	entity_.insert( entity_.begin() +ghost_ , NULL );
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

#if 0

void
Vertices::setMesh( const index_t k , MeshBase* m )
{
	ursa_assert( k<nb() );
	mesh_[k] = m;
}

void
Vertices::setEntity( const index_t k , Entity* e )
{
	ursa_assert( k<nb() );
	entity_[k] = e;
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
  		if (entity_[k]!=NULL)
  		{
  			geo = "-"+entity_[k]->objectClassName()+"-"+entity_[k]->memberTypeName();
  		}
  		else geo = "";
  		printf(" : %s , b[ %3d ] , g[ %p%s ] , u = (",(k<ghost_)? "GHST":"REAL",
  						body_[k],(void*)entity_[k],geo.c_str());
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
		if (entity_[k]!=NULL)
		{
			geo = "-"+entity_[k]->objectClassName()+"-"+entity_[k]->memberTypeName();
		}
		else geo = "";
		printf(" : %s , b[ %3d ] , g[ %p%s ] , u = (",(k<ghost_)? "GHST":"REAL",
						body_[k],(void*)entity_[k],geo.c_str());
		for (index_t d=0;d<udim_;d++)
			printf("%12.4e ",u(k)[d]);
		printf(")");
	}
	printf("\n");
}

template<typename type>
void
Vertices::print( Topology<type>& topology , index_t v ) const
{
  print(v,true);
  std::vector<index_t> ball;
  topology.allWithSubset( {v} , ball );
  for (index_t k=0;k<ball.size();k++)
  {
    std::string s = (topology.ghost(ball[k])) ? "ghost" : "real";
    printInline( topology.get(ball[k]) , "\telement["+stringify(ball[k])+"] " + s );
  }
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
	mesh_.erase( mesh_.begin() +k );
	entity_.erase( entity_.begin() + k );
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
			if (geometrics::distance(operator[](k),operator[](j),dim_)<tol)
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

			real d = geometrics::distance( (*this)[k] , (*this)[it->second] , dim_ );
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

void
Vertices::intersectGeometry( index_t n0 , index_t n1 , Entity*& e ) const
{
  Entity* e0 = entity(n0);
  Entity* e1 = entity(n1);

  if (e0==NULL || e1==NULL)
  {
    e = NULL;
    return;
  }
  e = e0->intersect(e1);
}

void
Vertices::findGeometry( const Body& body , index_t ibody ,real_t tol )
{
  // useful for when a mesh is read but it does not have a geometry
  // even though we know what it should be
  std::vector<real_t> x(4,0.);
  coord_t dlim = (dim_<=3) ? dim_ : 4;

  Entity *e;
  std::vector<Entity*> e_candidates;
  std::vector<std::vector<real_t>> u_candidates;

  bool recheck;

  // get the full list of entities
  std::vector<Entity*> entities;
  body.listTessellatableEntities(entities);
  std::vector<real_t> distances( entities.size() , 0. );
  std::vector<real_t> u( udim_*entities.size() );

  if (dim_ != 4) ursa_assert_msg(udim_ > 0, "udim = %u" , udim_ );
  std::vector<real_t> uk(udim_, INFTY);

  for (index_t k=0;k<nb();k++)
  {

    recheck = true;

    const real_t* xk = this->operator[](k);

    e_candidates.clear();
    u_candidates.clear();

    // compute the distance to each entity
    for (index_t j=0;j<entities.size();j++)
    {
      // re-assign the original coordinates
      for (coord_t d=0;d<dlim;d++)
        x[d] = xk[d];

      e = entities[j];

      // project to the entity and compute the distance
      entities[j]->projectToGeometry( x.data() , uk.data() );

      distances[j] = geometrics::distance2( xk , x.data() , dim_ );

      #if 0
      if (entities[j]->number()==2 && !entities[j]->tesseractGeometry())
      {
        // better to use inTopology for surface fit
        int icode = EG_inTopology( entities[j]->object() , xk );
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
      ursa_assert(entity_[k] == NULL);
			printf("no candidiates for vertex %lu\n",k);
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
    if (entity_[k] != NULL)
    {
      ursa_assert(ek == entity_[k]);
    }
    else
      setEntity(k,ek);

    setParam(k,uk);
    body_[k] = ibody;

    // get the entity coordinates
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
   real_t d = geometrics::distance2(xk,x.data(),dim_);
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

  std::vector<Entity*> e_candidates;
  std::vector<std::vector<real_t>> u_candidates;

  // get the full list of entities
  std::vector<Entity*> entities;
  body.listTessellatableEntities(entities);
  std::vector<real_t> distances( entities.size() , 0. );
  std::vector<real_t> u( udim_*entities.size() );

  if (dim_ != 4) ursa_assert(udim_ > 0);
  std::vector<real_t> uk(udim_);

  for (index_t k=0;k<nb();k++)
  {
    // find the closest entity with lowest topological dimension
    const real_t* xk = this->operator[](k);

    // compute the distance to each entity
    for (index_t j=0;j<entities.size();j++)
    {
      // re-assign the original coordinates
      for (coord_t d=0;d<dlim;d++)
        x[d] = xk[d];

      // project to the entity and compute the distance
      entities[j]->projectToGeometry(x.data(), u.data() + udim_*j);
      distances[j] = geometrics::distance2( xk , x.data() , dim_ );

    }

   real_t dmin = *std::min_element( distances.begin() , distances.end() );

    e_candidates.clear();
    u_candidates.clear();
    for (index_t j=0;j<distances.size();j++)
    {
      if (distances[j]<dmin+1e-12)
      {
        e_candidates.push_back( entities[j] );
        u_candidates.push_back( std::vector<real_t>(u.begin() + udim_*j, u.begin() + udim_*(j+1)) );
      }
    }

    if (e_candidates.size()==0)
		{
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

    //print(k,true);
    setEntity(k,ek);
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

void
Vertices::reserve( index_t nvert )
{
  x_.reserve( nvert*dim_ );
  u_.reserve( nvert*udim_ );
  body_.reserve( nvert );
  mesh_.reserve( nvert );
  entity_.reserve( nvert );
}

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
    if (entity_[k]==NULL) continue;
    geometry[k] = entity_[k]->bodyIndex();
  }
  J["geometry"] = geometry;
}

void
Vertices::fromJSON( const json& J , const Model* model )
{
  dim_ = J["dimension"];
  std::vector<real_t> X = J.at("coordinates");
  x_ = X;
  std::vector<real_t> U = J.at("parameters");
  u_ = U;
  ghost_ = J["nb_ghost"];

  entity_.resize( nb() , NULL );
  body_.resize( nb() );
  mesh_.resize( nb() );
  tag_.resize( nb() );

  if (model!=NULL)
  {
    std::vector<Entity*> entities;
    model->listEntities(entities);
    std::vector<int> geometry = J.at("geometry");
    ursa_assert(geometry.size()==nb());
    for (index_t k=0;k<nb();k++)
    {
      tag_[k] = 0;
      if (geometry[k]<0) continue;
      tag_[k] = 1;
      for (index_t j=0;j<entities.size();j++)
      {
        if ((int)entities[j]->bodyIndex()==geometry[j])
        {
          setEntity(k,entities[j]);
          break;
        }
      }
    }
  }
  else
  {
    std::vector<int> geometry = J.at("geometry");
    for (index_t k=0;k<nb();k++)
    {
      if (geometry[k]<0)
        tag_[k] = 0;
      else
        tag_[k] = 1;
    }
  }
}

#endif

void
Vertices::clear()
{
  x_.clear();
  u_.clear();
  body_.clear();
  mesh_.clear();
  entity_.clear();
  ghost_ = 0;
}

} // ursa
