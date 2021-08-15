#ifndef AVRO_SANDBOX_SDOT_VISUALIZE_H_
#define AVRO_SANDBOX_SDOT_VISUALIZE_H_

#include "common/parallel_for.h"
#include "common/set.h"

#include "library/meshb.h"
#include "numerics/linear_algebra.h"

#include "voronoi/optimal_transport.h"

typedef avro::real_t REAL;
typedef void VOID;
#define TETLIBRARY
#include <tetgen1.5.0/tetgen.h>

#include <iomanip>
#include <fstream>

namespace avro
{

using namespace voronoi;

bool
isclipped_kernel( const index_t* v , index_t nv , coord_t dim , const Points& points , const std::vector<real_t>& center , const std::vector<real_t>& normal )
{
  // determine if all points are on the appropriate side of the hyperplane
  int s0 = 0;
  for (index_t j = 0; j < nv; j++)
  {
    // compute the value and sign of the dot product
    real_t dp = 0.0;
    for (coord_t d = 0; d < dim; d++)
      dp += ( points[ v[j] ][d] - center[d] ) * normal[d];
    int s = ( dp < 0.0 ) ? -1 : 1;

    if (j == 0)
    {
      s0 = s;
      continue;
    }

    // if the sign is different for any two vertices, then we have an intersection
    if (s*s0 < 0)
    {
      return true;
    }
  }
  return false;
}

template<typename type>
class ClippedCell : public Topology<Simplex>
{
public:
  ClippedCell( const LaguerreDiagram<type>& diagram , coord_t dim , index_t cell ,
      const matd<real_t>& B , const matd<real_t>& At ) :
    Topology<Simplex>(points_,3),
    diagram_(diagram),
    dim_(dim),
    points_(3),
    cell_(cell),
    B_(B), At_(At)
  {}

  void clip( const std::vector<real_t>& center , const std::vector<real_t>& normal )
  {
    // retrieve the laguerre cell
    const LaguerreCell<type>& polytope = diagram_.cell(cell_);

    for (index_t k = 0; k < polytope.nb(); k++)
    {
      // determine if this polytope is clipped by the plane
      bool clipped = isclipped_kernel( polytope(k) , polytope.nv(k) , dim_ , polytope.points() , center , normal );
      if (!clipped) continue;

      std::vector<real_t> polyhedron;

      // compute the edges of this polytope
      std::vector<index_t> edges;
      for (index_t i = 0; i < polytope.nv(k); i++)
      for (index_t j = i+1; j < polytope.nv(k); j++)
      {
        if (polytope.element().is_edge( polytope(k,i) , polytope(k,j) ) )
        {
          edges.push_back( polytope(k,i) );
          edges.push_back( polytope(k,j) );
        }
      }

      Table<int> incidence(TableLayout_Jagged);

      // compute the intersection point of every edge
      for (index_t i = 0; i < edges.size()/2; i++)
      {
        index_t e0 = edges[2*i];
        index_t e1 = edges[2*i+1];

        // compute the intersection point
        real_t num = 0.0, den = 0.0;
        for (coord_t d = 0; d < dim_; d++)
        {
          num += ( center[d] - polytope.points()[e0][d] )*normal[d];
          den += ( polytope.points()[e1][d] - polytope.points()[e0][d] )*normal[d];
        }
        if (den == 0.0) continue;
        real_t s = num / den;

        // check if there is an intersection
        if (s < 0.0 || s > 1.0) continue;

        // compute the geomeetric intersection
        vecd<real_t> xs(dim_);
        for (coord_t d = 0; d < dim_; d++)
          xs(d) = polytope.points()[e0][d] + s*( polytope.points()[e1][d] - polytope.points()[e0][d] );

        // compute the uvw coordinates
        vecd<real_t> xu(3);
        xu = B_* (At_*xs );

        for (coord_t d = 0; d < xu.m(); d++)
          polyhedron.push_back( xu(d) );

        // add the incidence relations for this point (to extract edges)
        std::vector<int> b0 = polytope.points().incidence().get(e0);
        std::vector<int> b1 = polytope.points().incidence().get(e1);

        std::vector<int> bu;
        Set::intersection( b0 , b1 , bu );
        incidence.add( bu.data() , bu.size() );
      }

      // compute the tetrahedralization of the polyhedron
      tetgenio input,output;

      // fill the vertices
      input.numberofpoints = polyhedron.size() / 3;
      input.pointlist = new REAL[input.numberofpoints*3];
      for (index_t i = 0; i < polyhedron.size();i++)
        input.pointlist[i] = polyhedron[i];

      std::string switches = "Q";
      try
      {
        tetrahedralize( (char*)switches.c_str() , &input , &output );
      }
      catch(...)
      {
        printf("there was an error with tetgen :(\n");
      }

      Polytope element(3,1,incidence);

      // add the tetrahedra
      index_t nb_points = points_.nb();
      for (index_t j = 0; j < (index_t)output.numberoftetrahedra; j++)
      {
        std::vector<index_t> tet( output.tetrahedronlist+4*j , output.tetrahedronlist+4*(j+1) );
        for (index_t d = 0; d < 4; d++)
          tet[d] += nb_points;
        add(tet.data(),tet.size());
        tet2site_.push_back( polytope.cell2site(k) );
      }

      // add the points
      for (index_t j = 0; j < polyhedron.size() / 3; j++)
        points_.create( &polyhedron[3*j] );

      avro_assert( polyhedron.size() / 3 == incidence.nb() );
      for (index_t j = 0; j < incidence.nb(); j++)
      for (index_t i = j+1; i < incidence.nb(); i++)
      {
        if (element.is_edge( i , j ) )
        {
          // add this edge
          edges_.push_back(i);
          edges_.push_back(j);
        }
      }
    }
  }

  index_t tet2site( index_t j ) const { return tet2site_[j]; }

  const std::vector<index_t>& edges() const { return edges_; }

private:
  const LaguerreDiagram<type>& diagram_;
  coord_t dim_;
  Points points_;
  index_t cell_; // which laguerre cell is this from?
  index_t elem_; // which domain element?
  index_t site_; // which voronoi site does this correspond to?

  const matd<real_t>& B_;
  const matd<real_t>& At_;

  std::vector<index_t> tet2site_;

  std::vector<index_t> edges_;
};

template<typename type>
class HyperSlice
{
  typedef HyperSlice thisclass;

public:
  HyperSlice( const LaguerreDiagram<type>& diagram ) :
    diagram_(diagram),
    dim_(diagram.number()), // not dim because they could be embedded in dim+1 for sdot
    uvw_(3),
    tetrahedra_(uvw_,3),
    edges_(uvw_,1),
    normal_(dim_,0)
  {}

  void isclipped( index_t k )
  {
    clipped_[k] = isclipped_kernel( diagram_(k) , diagram_.nv(k) , dim_ , diagram_.points() , center_ , normal_ );
  }

  void clip( index_t k )
  {
    // now actually do the clipping
    cells_[k]->clip(center_,normal_);
  }

  void compute( const std::vector<real_t>& c , const coord_t dir )
  {
    center_.assign( c.begin() , c.end() );

    std::vector<real_t> e0(dim_,0),e1(dim_,0),e2(dim_,0),e3(dim_,0);
    e0[0] = 1;
    e1[1] = 1;
    e2[2] = 1;
    e3[3] = 1;
    for (coord_t d = 0; d < dim_; d++)
      normal_[d] = 0;
    normal_[dir] = 1;

    real_t* u = nullptr;
    real_t* v = nullptr;
    real_t *w = nullptr;
    if (dir == 0)
    {
      u = e1.data();
      v = e2.data();
      w = e3.data();
    }
    else if (dir == 1)
    {
      u = e0.data();
      v = e2.data();
      w = e3.data();
    }
    else if (dir == 2)
    {
      u = e0.data();
      v = e1.data();
      w = e3.data();
    }
    else if (dir == 3)
    {
      u = e0.data();
      v = e1.data();
      w = e2.data();
    }

    // determine all cells which are clipped
    clipped_.clear();
    clipped_.resize( diagram_.nb() , false );
    ProcessCPU::parallel_for(
      parallel_for_member_callback( this , &thisclass::isclipped ),
      0,clipped_.size() );


    // precompute the transformation matrix
    matd<real_t> A(4,3);
    matd<real_t> At(3,4);
    matd<real_t> B(3,3);

    for (coord_t d = 0; d < 4; d++)
    {
      A(d,0) = u[d];
      A(d,1) = v[d];
      A(d,2) = w[d];
      At(0,d) = u[d];
      At(1,d) = v[d];
      At(2,d) = w[d];
    }
    B = At*A;
    B = numerics::inverse( B );
    //std::cout << A << std::endl;
    //std::cout << B << std::endl;

    // create the clipped cells
    cells_.clear();
    cells_.reserve( clipped_.size() );
    for (index_t k = 0; k < clipped_.size(); k++)
    {
      if (!clipped_[k]) continue;
      cells_.push_back( std::make_shared<ClippedCell<type>>(diagram_,dim_,k,B,At) );
    }

    printf("need to clip %lu cells out of %lu\n",cells_.size(),diagram_.nb());

    // clip the cells
    ProcessCPU::parallel_for(
      parallel_for_member_callback( this , &thisclass::clip ),
      0,cells_.size() );

    for (index_t k = 0; k < cells_.size(); k++)
    {
      const ClippedCell<type>& cell = *cells_[k].get();

      index_t nb_points = uvw_.nb();
      for (index_t j = 0; j < cell.nb(); j++)
      {
        std::vector<index_t> tet = cell.get(j);
        for (coord_t d = 0; d < 4; d++)
          tet[d] += nb_points;
        tetrahedra_.add( tet.data() , tet.size() );
        tet2site_.push_back( cell.tet2site(j) );
      }

      const std::vector<index_t>& edges = cell.edges();
      for (index_t j = 0; j < edges.size() / 2; j++)
      {
        index_t e0 = edges[2*j] + nb_points;
        index_t e1 = edges[2*j+1] + nb_points;
        index_t e[2] = {e0,e1};
        edges_.add( e , 2 );
      }

      for (index_t j = 0; j < cell.points().nb(); j++)
        uvw_.create( cell.points()[j] );
    }
  }

  void save( const std::string& prefix )
  {
    // export the mesh
    library::meshb writer;
    writer.open( 3 , prefix + "_tet.mesh" );
    writer.write( uvw_ );
    std::vector<index_t> refs( tetrahedra_.nb() , 0 );
    writer.write( tetrahedra_ , refs );
    writer.close();

    // export the sites
    json J;
    J["field"] = tet2site_;
    std::ofstream output(prefix+"_sites.json");
    output << std::setw(4) << J << std::endl;
  }

  const Topology<Simplex>& tetrahedra() const { return tetrahedra_; }
  Topology<Simplex>& tetrahedra() { return tetrahedra_; }
  const std::vector<index_t>& tet2site() const { return tet2site_; }
  const Topology<Simplex>& edges() const { return edges_; }

private:
  const LaguerreDiagram<type>& diagram_;
  coord_t dim_;
  Points uvw_;
  Topology<Simplex> tetrahedra_;
  Topology<Simplex> edges_;

  std::vector<bool> clipped_;
  std::vector<std::shared_ptr<ClippedCell<type>>> cells_;

  std::vector<real_t> normal_;
  std::vector<real_t> center_;

  std::vector<index_t> tet2site_;
};

class SliceSites : public Field<Simplex,real_t>
{
public:
  SliceSites( Topology<Simplex>& slice , const std::vector<index_t>& sites ) :
    Field<Simplex,real_t>(slice,0,DISCONTINUOUS)
  {
    this->build();
    this->element().set_basis( BasisFunctionCategory_Lagrange );
    for (index_t k=0;k<slice.nb();k++)
    {
      this->value(k) = sites[k];
    }
  }

  index_t nb_rank() const { return 1; }

  std::vector<std::string> ranknames() const
   {std::vector<std::string> result; result.push_back("sites"); return result;}
};

}

#endif
