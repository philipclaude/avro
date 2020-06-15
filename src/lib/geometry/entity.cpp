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
#include "common/tree.hpp"

#include "geometry/entity.h"

#include "shape/simplex.h"

#include "mesh/topology.h"

#include "numerics/coordinate.h"

#include <algorithm>
#include <set>

namespace avro
{

template <typename T>
int
contains( const T& x , const std::vector<T>& X )
{
	const index_t n = X.size();
	for (index_t k=0;k<n;k++)
		if (X[k]==x) return int(k);
	return -1;
}

template<typename T>
void
intersection( const std::vector<T>& s0 , const std::vector<T>& s1 , std::vector<T>& p )
{
	p.clear();
	const index_t n0 = s0.size();

	p.clear();
	for (index_t i=0;i<n0;i++)
	{
		if ( contains( s0[i] , s1 ) > -1 )
			p.push_back(s0[i]);
	}
	uniquify(p);
}

Entity::Entity( coord_t number ) :
  number_(number),
  name_("unnamed"),
  parent_(nullptr),
  identifier_(0),
	sense_required_(false),
  tessellatable_(false),
  egads_(false),
	sign_(1)
{}

Entity::Entity( coord_t number , const std::string& name ) :
  number_(number),
  name_(name),
  parent_(nullptr),
  identifier_(0),
	sense_required_(false),
  tessellatable_(false),
  egads_(false),
	sign_(1)
{}

Entity*
Entity::intersect( Entity* e1 )
{
  std::vector<Entity*> p0( nb_parents()+1 );
  p0[0] = this;
  for (index_t k=0;k<nb_parents();k++)
    p0[k+1] = parents_[k];

  std::vector<Entity*> p1( e1->nb_parents()+1 );
  p1[0] = e1;
  for (index_t k=0;k<e1->nb_parents();k++)
    p1[k+1] = e1->parents(k);

  // intersect all the parents, including this and e1 itself
  std::vector<Entity*> s;
  intersection( p0 , p1 , s );

  // find the lowest dimension entity
  if (s.size()==0) return NULL;

  Entity* e = NULL;
  for (index_t k=0;k<s.size();k++)
  {
    if (!s[k]->tessellatable()) continue;

    if (e==NULL) e = s[k];
    else if (s[k]->number()<e->number())
    {
      e = s[k];
    }
  }

	if (e==NULL)
	{
		/*printf("parents 0:\n");
		for (index_t k=0;k<p0.size();k++)
			p0[k]->print_header();
		printf("parents 1:\n");
		for (index_t k=0;k<p1.size();k++)
			p1[k]->print_header();*/
	}


  if (e==NULL) return NULL;
  if (!e->tessellatable()) return NULL;

  return e;
}

Entity*
Entity::intersect( Entity* e1 , Entity* e2 , bool only_check )
{
  #if 0
  Entity* e = intersect(e1);
  if (e==NULL) return NULL;
  return e->intersect(e2);
  #else
  // we need to look for a common face between all three entities
  // there should be only 1 (if any)
  std::set<Entity*> p0,p1,p2;
  if (number_==2) p0.insert(this);
  if (e1->number()==2) p1.insert(e1);
  if (e2->number()==2) p2.insert(e2);

  for (index_t k=0;k<nb_parents();k++)
  {
    if (!parents(k)->tessellatable() || parents(k)->number()!=2) continue;
    p0.insert( parents(k) );
  }
  for (index_t k=0;k<e1->nb_parents();k++)
  {
    if (!e1->parents(k)->tessellatable() || e1->parents(k)->number()!=2) continue;
    p1.insert( e1->parents(k) );
  }
  for (index_t k=0;k<e2->nb_parents();k++)
  {
    if (!e2->parents(k)->tessellatable() || e2->parents(k)->number()!=2) continue;
    p2.insert( e2->parents(k) );
  }

  // it is possible to not find a common face, initialize to NULL
  Entity* face = NULL;
  std::set<Entity*>::iterator it,it0;
  for (it=p0.begin();it!=p0.end();it++)
  {
    if (p1.find(*it)!=p1.end() && p2.find(*it)!=p2.end())
    {
      if (face!=NULL && !only_check)
      {
        // we cannot have duplicated faces
        // print some debug information
        face->print();

				print();
        printf("parents 0:\n");
        for (it0=p0.begin();it0!=p0.end();it0++)
          (*it0)->print();
				e1->print();
        printf("parents 1:\n");
        for (it0=p1.begin();it0!=p1.end();it0++)
          (*it0)->print();
				e2->print();
        printf("parents 2:\n");
        for (it0=p2.begin();it0!=p2.end();it0++)
          (*it0)->print();

        avro_assert_not_reached;
      }
      else if (face!=NULL && only_check)
      {
        return NULL;
      }
      face = *it;
    }
  }
  return face;
  #endif
}

Entity*
Entity::intersect( Entity* e1 , Entity* e2 , Entity* e3 )
{
  #if 0
  Entity* e = intersect(e1,e2);
  if (e==NULL) return NULL;
  return e->intersect(e3);
  #else
  // we need to look for a common face between all three entities
  // there should be only 1 (if any)
  std::set<Entity*> p0,p1,p2,p3;
  if (number_==3) p0.insert(this);
  if (e1->number()==3) p1.insert(e1);
  if (e2->number()==3) p2.insert(e2);
  if (e3->number()==3) p3.insert(e3);

  for (index_t k=0;k<nb_parents();k++)
  {
    if (!parents(k)->tessellatable() || parents(k)->number()!=3) continue;
    p0.insert( parents(k) );
  }
  for (index_t k=0;k<e1->nb_parents();k++)
  {
    if (!e1->parents(k)->tessellatable() || e1->parents(k)->number()!=3) continue;
    p1.insert( e1->parents(k) );
  }
  for (index_t k=0;k<e2->nb_parents();k++)
  {
    if (!e2->parents(k)->tessellatable() || e2->parents(k)->number()!=3) continue;
    p2.insert( e2->parents(k) );
  }
  for (index_t k=0;k<e3->nb_parents();k++)
  {
    if (!e3->parents(k)->tessellatable() || e3->parents(k)->number()!=3) continue;
    p3.insert( e3->parents(k) );
  }

  // it is possible to not find a common face, initialize to NULL
  Entity* vol = NULL;
  std::set<Entity*>::iterator it,it0;
  for (it=p0.begin();it!=p0.end();it++)
  {
    if (p1.find(*it)!=p1.end() && p2.find(*it)!=p2.end() && p3.find(*it)!=p3.end())
    {
      if (vol!=NULL)
      {
        // we cannot have duplicated Volumes
        // print some debug information
        vol->print();

        printf("parents 0:\n");
        for (it0=p0.begin();it0!=p0.end();it0++)
          (*it0)->print();
        printf("parents 1:\n");
        for (it0=p1.begin();it0!=p1.end();it0++)
          (*it0)->print();
        printf("parents 2:\n");
        for (it0=p2.begin();it0!=p2.end();it0++)
          (*it0)->print();
        for (it0=p3.begin();it0!=p3.end();it0++)
          (*it0)->print();
        avro_assert_not_reached;
      }
      vol = *it;
    }
  }
  return vol;
  #endif
}

bool
Entity::above( const Entity* e ) const
{
	for (index_t k=0;k<nb_children();k++)
	{
		if (child_ptr(k)==e) return true;
		if (child(k).above(e))
			return true;
	}
	return false;
}

template<typename type>
void
Entity::construct( std::shared_ptr<Topology<type>>& node , Topology<type>& root ) const
{
	node = std::make_shared<Topology<type>>(root.points(),number_,root.shape().order());
}

void
Entity::print_header() const
{
	printf("entity %p: number = %u: %s\n",(void*)this,number_,name_.c_str());
}

template class Tree<Entity>;
template void Entity::construct( std::shared_ptr<Topology<Simplex>>& node , Topology<Simplex>& root ) const;

} // avro
