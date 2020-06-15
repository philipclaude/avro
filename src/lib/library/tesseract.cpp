//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "common/table.h"
#include "common/tools.h"

#include "geometry/psc/facet.h"
#include "geometry/psc/node.h"

#include "library/tesseract.h"

#include "shape/polytope.h"

namespace avro
{

namespace library
{

static void
getNodeIndices( const Entity* e , std::vector<index_t>& N )
{
  for (index_t k=0;k<e->nb_children();k++)
    getNodeIndices( e->child_ptr(k) , N );
  if (e->number()==0)
    N.push_back(e->identifier());
}

static index_t
identify( std::vector<index_t>& nodes )
{
  std::sort(nodes.begin(),nodes.end());
  index_t f = 0;
  for (index_t j=0;j<nodes.size();j++)
    f += nodes[j]*pow(16,j);
  return f;
}

static int
distancei( int* x , int* y , coord_t dim )
{
  int d = 0;
  for (coord_t i=0;i<dim;i++)
    d += (x[i] -y[i])*(x[i] -y[i]);
  return d;
}

static std::string
label( std::vector<index_t>& f )
{
  avro_assert( f.size()>0 );
  std::sort( f.begin() , f.end() );
  std::string s = stringify(f[0]);
  for (index_t i=1;i<f.size();i++)
    s += "|"+stringify(f[i]);
  return s;
}

index_t
findEdgeIndex( index_t p0 , index_t p1 , std::vector<index_t> edges )
{
  for (index_t k=0;k<edges.size()/2;k++)
  {
    if (edges[2*k]==p0 && edges[2*k+1]==p1) return k;
    if (edges[2*k]==p1 && edges[2*k+1]==p0) return k;
  }
  avro_assert_not_reached;
  return 0;
}

void
Tesseract::build()
{
  typedef std::shared_ptr<Entity> Entity_ptr;
  std::vector<Entity_ptr> cube_entities(8);
  std::vector<Entity_ptr> square_entities(24);
  std::vector<Entity_ptr> edge_entities(32);
  std::vector<Entity_ptr> node_entities(16);

  // construct the 8 bisectors with the corresponding coordinate they are on
  std::vector<int> bisector = {-1,1,-2,2,-3,3,-4,4};

  // create the signs of all vertex coordinates
  Table<int> P(TableLayout_Rectangular,4);
  int X0[4] =  { 1, 1, 1, 1}; P.add(X0,4);
  int X1[4] =  {-1, 1, 1, 1}; P.add(X1,4);
  int X2[4] =  {-1,-1, 1, 1}; P.add(X2,4);
  int X3[4] =  {-1,-1,-1, 1}; P.add(X3,4);
  int X4[4] =  {-1, 1,-1, 1}; P.add(X4,4);
  int X5[4] =  {-1, 1, 1,-1}; P.add(X5,4);
  int X6[4] =  {-1,-1, 1,-1}; P.add(X6,4);
  int X7[4] =  {-1, 1,-1,-1}; P.add(X7,4);
  int X8[4] =  { 1,-1, 1, 1}; P.add(X8,4);
  int X9[4] =  { 1,-1,-1, 1}; P.add(X9,4);
  int X10[4] = { 1, 1,-1,-1}; P.add(X10,4);
  int X11[4] = { 1,-1, 1,-1}; P.add(X11,4);
  int X12[4] = { 1,-1,-1,-1}; P.add(X12,4);
  int X13[4] = { 1, 1, 1,-1}; P.add(X13,4);
  int X14[4] = { 1, 1,-1, 1}; P.add(X14,4);
  int X15[4] = {-1,-1,-1,-1}; P.add(X15,4);

  // sanity check each vertex is unique
  avro_assert( P.nb()==16 );
  for (index_t k=0;k<P.nb();k++)
  {
    for (index_t j=k+1;j<P.nb();j++)
    {
      int dist = distancei( P(k) , P(j) , 4 );
      avro_assert_msg( dist>0 , "k = %lu , j = %lu , dist = %d\n",k,j,dist );
    }
  }

  // now create the tesseract node entities
  real_t x[4];
  for (index_t k=0;k<P.nb();k++)
  {
    // map the coordinates using the appropriate sign
    for (index_t d=0;d<4;d++)
      x[d] = x0_[d] +P(k,d)*.5*length_[d];

    node_entities[k] = std::make_shared<PSC::Node>(this,x);

    node_entities[k]->set_identifier(k);
  }

  // add the facets to each vertex
  Table<int> vfm(TableLayout_Jagged);
  for (index_t k=0;k<P.nb();k++)
  {
    std::vector<int> facets;
    for (index_t i=0;i<bisector.size();i++)
    {
      // decode the bisector
      int b = bisector[i];
      coord_t dim = std::abs(b) -1;
      int x = (b<0) ? -1 : 1;

      // check if the coordinate matches
      if (P(k,dim)==x)
        facets.push_back( b );
    }

    // the number of facets attached to each vertex should be 4
    avro_assert(facets.size()==4);
    vfm.add(facets.data(),4);
  }

  // the tesseract "mesh" (1 polytopal element)
  std::vector<index_t> T = linspace(16);

  // retrieve the bounding cubes
  Polytope refTess(4,1,vfm);
  std::vector<int> facets3;
  refTess.hrep(T.data(),16,facets3);
  avro_assert_msg( facets3.size()==8 ,
    "nb_cubes = %lu but should = %d",facets3.size(),8 );

  // get the edges
  std::vector<index_t> edges;
  refTess.edges( T.data() , T.size() , edges );
  avro_assert( edges.size()==64 ); // 32 edges

  // create the edges from the nodes
  for (index_t k=0;k<edges.size()/2;k++)
  {
    std::vector<Entity_ptr> edge = {node_entities[edges[2*k]],node_entities[edges[2*k+1]] };
    edge_entities[k] = std::make_shared<PSC::Facet>( this, edge );
  }

  // loop through all the bisectors
  std::map<std::string,index_t> square_labels;
  Table<index_t> cubes(TableLayout_Jagged);
  Table<index_t> squares(TableLayout_Jagged);
  for (index_t k=0;k<bisector.size();k++)
  {
    Polytope refCube(3,1,vfm);

    // first construct the cube
    std::vector<index_t> cube;
    refCube.vrep( T.data() , T.size() , bisector[k] , cube );

    cubes.add( cube.data() , cube.size() );

    avro_assert( cube.size()==8 );

    // get the squares bounding each cube
    std::vector<int> facets;
    refCube.hrep( cube.data() , cube.size() , facets );

    // we should get 7 facets (the expected 6 plus the 3-facet we are on)
    avro_assert( facets.size()==7 );
    std::vector<index_t> square_index;
    for (index_t j=0;j<facets.size();j++)
    {
      // skip the redundant bisector
      if (facets[j]==bisector[k]) continue;

      // ask the ref square for the vrep along this facet
      Polytope refSquare(2,1,vfm);
      std::vector<index_t> square;
      refSquare.vrep( cube.data() , cube.size() , facets[j] , square );

      avro_assert( square.size()==4 );

      // generate a unique label for the square
      std::string s = label(square);
      std::map<std::string,index_t>::iterator it = square_labels.find(s);
      if (it==square_labels.end())
      {
        index_t idx = square_labels.size();

        // the square index is the current size of the set
        square_index.push_back( idx );
        square_labels.insert( std::pair<std::string,index_t>(s,idx) );
        squares.add( square.data() , square.size() );

        // get the facets of this square (edges)
        std::vector<int> squareFacets;
        refSquare.hrep( square.data() , square.size() , squareFacets );

        // we should get 4 edges + two redundant ones from the cube and square
        avro_assert( squareFacets.size()==6 );

        std::vector<Entity_ptr> square_edges;
        for (index_t i=0;i<squareFacets.size();i++)
        {
          // skip the redundant bisectors
          if (squareFacets[i]==facets[j] || squareFacets[i]==bisector[k])
            continue;

          // retrieve the edge of the square
          std::vector<index_t> ve;
          refSquare.vrep( square.data() , square.size() , squareFacets[i] , ve );
          avro_assert( ve.size()==2 );


          // look for the edge
          index_t edgeIndex = findEdgeIndex( ve[0] , ve[1] , edges );

          square_edges.push_back( edge_entities[edgeIndex] );

        }

        avro_assert(square_edges.size()==4);

        // create the square
        square_entities[idx] = std::make_shared<PSC::Facet>( this , square_edges );
      }
      else
      {
        // get the square index, its edges were referenced upon construction
        index_t idx = it->second;
        square_index.push_back( idx );
      }

    }

    // add the box which references the square indices
    avro_assert( square_index.size()==6 );
    std::vector<Entity_ptr> cube_squares( 6 );
    for (index_t ii=0;ii<6;ii++)
      cube_squares[ii] = square_entities[square_index[ii]];
    cube_entities[k] = std::make_shared<PSC::Facet>( this, cube_squares );
  }
  avro_assert( cubes.nb()==8 );
  avro_assert( squares.nb()==24 );
  avro_assert( square_labels.size()==24 );
  avro_assert( cube_entities.size()==8 );

  for (index_t k=0;k<cube_entities.size();k++)
    add(cube_entities[k]);

  build_parents();

  std::map<index_t,Entity*> id2entity;
  std::vector<index_t> identifiers;

  // set the levels (for printing)
  for (index_t k=0;k<cube_entities.size();k++)
  {
    //TODOcube_entities[k]->setLevel(1);

    std::vector<index_t> node_indices;
    getNodeIndices(cube_entities[k].get(),node_indices);
    uniquify(node_indices);
    avro_assert_msg( node_indices.size()==8 , "|node_indices| = %lu" , node_indices.size() );

    index_t id = identify(node_indices);
    id2entity[id] = cube_entities[k].get();
    identifiers.push_back(id);

    cube_entities[k]->set_identifier(k+1);
    //TODOcube_entities[k]->setGraphLabel("V "+stringify(k));
  }

  for (index_t k=0;k<square_entities.size();k++)
  {
    //TODOsquare_entities[k]->setLevel(2);

    std::vector<index_t> node_indices;
    getNodeIndices(square_entities[k].get(),node_indices);
    uniquify(node_indices);
    avro_assert( node_indices.size()==4 );

    index_t id = identify(node_indices);
    id2entity[id] = square_entities[k].get();
    identifiers.push_back(id);

    square_entities[k]->set_identifier(k+1);
    //TODOsquare_entities[k]->setGraphLabel("F "+stringify(k));
  }

  for (index_t k=0;k<edge_entities.size();k++)
  {
    //TODOedge_entities[k]->setLevel(3);

    std::vector<index_t> node_indices;
    getNodeIndices(edge_entities[k].get(),node_indices);
    uniquify(node_indices);
    avro_assert( node_indices.size()==2 );

    index_t id = identify(node_indices);
    id2entity[id] = edge_entities[k].get();
    identifiers.push_back(id);

    edge_entities[k]->set_identifier(k+1);
    //TODOedge_entities[k]->setGraphLabel("E "+stringify(k));
  }

  for (index_t k=0;k<node_entities.size();k++)
  {
    //TODOnode_entities[k]->setLevel(4);
    //node_entities[k]->set_identifier(k+1);
    index_t id = k;
    id2entity[id] = node_entities[k].get();
    identifiers.push_back(id);
    //TODOnode_entities[k]->setGraphLabel("N "+stringify(k));
  }

  // adjust the body indices to account for the unique labels
  std::sort( identifiers.begin() , identifiers.end() );
  for (index_t k=0;k<identifiers.size();k++)
  {
    index_t pid = id2entity[identifiers[k]]->identifier();
    id2entity[identifiers[k]]->set_identifier(k+1);
    if (id2entity[identifiers[k]]->number()==3)
      printf("entity %lu (%lu) maps to %lu\n",pid,identifiers[k],id2entity.size()-k );
  }

  // I really want the body indices of the cubes to come first
  // specifically because this makes the boundary group association easier
  // for codes that will call avro...in order to know which boundary group
  // is associated with the particular cube, just call identifier -1
  // so flip everything around (since cubes have the highest sorting number)
  index_t nb_entities = id2entity.size();
  for (index_t k=0;k<nb_entities;k++)
    id2entity[identifiers[k]]->set_identifier( nb_entities-k );

  //name_ = "TESSERACT";

  // check all cubes have no parents
  for (index_t k=0;k<cube_entities.size();k++)
  {
    avro_assert_msg( cube_entities[k]->nb_parents()==0 ,
      "nb_parents = %lu" , cube_entities[k]->nb_parents() );
  }

  // check the parents of all squares
  for (index_t k=0;k<square_entities.size();k++)
  {
    std::vector<index_t> count(4);
    for (index_t j=0;j<square_entities[k]->nb_parents();j++)
      count[ square_entities[k]->parents(j)->number() ]++;

    avro_assert( count[0]==0 );
    avro_assert( count[1]==0 );
    avro_assert( count[2]==0 );
    avro_assert( count[3]==2 );
  }

  // check the parents of all edges
  for (index_t k=0;k<edge_entities.size();k++)
  {
    std::vector<index_t> count(4);
    for (index_t j=0;j<edge_entities[k]->nb_parents();j++)
      count[ edge_entities[k]->parents(j)->number() ]++;

    avro_assert( count[0]==0 );
    avro_assert( count[1]==0 );
    avro_assert( count[2]==3 ); // three squares touch this edge
    avro_assert( count[3]==3 ); // three cubes touch this edge
  }

  // check the parents of all nodes
  for (index_t k=0;k<node_entities.size();k++)
  {
    std::vector<index_t> count(4);
    for (index_t j=0;j<node_entities[k]->nb_parents();j++)
      count[ node_entities[k]->parents(j)->number() ]++;

    avro_assert( count[0]==0 );
    avro_assert( count[1]==4 ); // four edges touch this node
    avro_assert( count[2]==6 ); // six squares touch this node
    avro_assert( count[3]==4 ); // four cubes touch this node
  }

}

} // library

} // avro
