#include "element/trace2cell.h"

namespace avro
{

void
mapTetTraceCoordsToCanonical( const real_t& s0 , const real_t& t0 , const real_t& u0,
                              const int orientation , real_t& s, real_t& t, real_t& u )
{
  // see accompanying matlab script which runs through all the
  // permutations of 1-2-3-4 (there are 24)
  // these are written out in the order they appear in the script
  // TODO use std::next_permutation along with the math
  switch (orientation)
  {
    case 1: // 0-1-2-3
    {
      s = s0;
      t = t0;
      u = u0;
      break;
    }
    case -1: // 0-1-3-2
    {
      s = s0;
      t = u0;
      u = t0;
      break;
    }
    case -2: // 0-2-1-3
    {
      s = t0;
      t = s0;
      u = u0;
      break;
    }
    case 2: // 0-2-3-1
    {
      s = t0;
      t = u0;
      u = s0;
      break;
    }
    case 3: // 0-3-1-2
    {
      s = u0;
      t = s0;
      u = t0;
      break;
    }
    case -3: // 0-3-2-1
    {
      s = u0;
      t = t0;
      u = s0;
      break;
    }
    case -4: // 1-0-2-3
    {
      s = 1 -s0 -t0 -u0;
      t = t0;
      u = u0;
      break;
    }
    case 4: // 1-0-3-2
    {
      s = 1 -s0 -t0 -u0;
      t = u0;
      u = t0;
      break;
    }
    case 5: // 1-2-0-3
    {
      s = t0;
      t = 1 -s0 -t0 -u0;
      u = u0;
      break;
    }
    case -5: // 1-2-3-0
    {
      s = t0;
      t = u0;
      u = 1 -s0 -t0 -u0;
      break;
    }
    case -6: // 1-3-0-2
    {
      s = u0;
      t = 1 -s0 -t0 -u0;
      u = t0;
      break;
    }
    case 6: // 1-3-2-0
    {
      s = u0;
      t = t0;
      u = 1 -s0 -t0 -u0;
      break;
    }
    case 7: // 2-0-1-3
    {
      s = 1 -s0 -t0 -u0;
      t = s0;
      u = u0;
      break;
    }
    case -7: // 2-0-3-1
    {
      s = 1 -s0 -t0 -u0;
      t = u0;
      u = s0;
      break;
    }
    case -8: // 2-1-0-3
    {
      s = s0;
      t = 1 -s0 -t0 -u0;
      u = u0;
      break;
    }
    case 8: // 2-1-3-0
    {
      s = s0;
      t = u0;
      u = 1 -s0 -t0 -u0;
      break;
    }
    case 9: // 2-3-0-1
    {
      s = u0;
      t = 1 -s0 -t0 -u0;
      u = s0;
      break;
    }
    case -9: // 2-3-1-0
    {
      s = u0;
      t = s0;
      u = 1 -s0 -t0 -u0;
      break;
    }
    case -10: // 3-0-1-2
    {
      s = 1 -s0 -t0 -u0;
      t = s0;
      u = t0;
      break;
    }
    case 10: // 3-0-2-1
    {
      s = 1 -s0 -t0 -u0;
      t = t0;
      u = s0;
      break;
    }
    case 11: // 3-1-0-2
    {
      s = s0;
      t = 1 -s0 -t0 -u0;
      u = t0;
      break;
    }
    case -11: // 3-1-2-0
    {
      s = s0;
      t = t0;
      u = 1 -s0 -t0 -u0;
      break;
    }
    case -12: // 3-2-0-1
    {
      s = t0;
      t = 1 -s0 -t0 -u0;
      u = s0;
      break;
    }
    case 12: // 3-2-1-0
    {
      s = t0;
      t = s0;
      u = 1 -s0 -t0 -u0;
      break;
    }
    default:
      avro_assert_not_reached;
  }
}

void
mapTetTraceCoordsToPentatope( const real_t* ut , int trace , real_t* uc )
{
  // the faces are properly oriented (see test in avro TraceToCellRefCoord_ut)
  // the intuition behind this is that the coordinate axes are defined starting
  // from the first vertex, then connecting to the second, third, etc.
  // i.e. facet (tet) 0-2-4-3 has edges (0,2), (0,4) and (0,3) defining
  // the coordinate axes as they are in the canonical trace (tet)
  // therefore, (0,2) means that t gets the first reference coordinate (s)
  // then (0,4) means the v gets the second (t) and (0,3) means u gets the
  // third reference coordinate (u).
  const real_t s = ut[0];
  const real_t t = ut[1];
  const real_t u = ut[2];
  switch (trace)
  {
  case 0:      // 1-2-4-3, 1+s+t+u = 0 plane
  {
    uc[0] = 1 - s - t - u;
    uc[1] = s;
    uc[2] = u;
    uc[3] = t;
    break;
  }
  case 1:      // 0-2-3-4, s = 0 plane
  {
    uc[0] = 0;
    uc[1] = s;
    uc[2] = t;
    uc[3] = u;
    break;
  }
  case 2:      // 0-1-4-3, t = 0 plane
  {
    uc[0] = s;
    uc[1] = 0;
    uc[2] = u;
    uc[3] = t;
    break;
  }
  case 3:      // 0-1-2-4, u = 0 plane
  {
    uc[0] = s;
    uc[1] = t;
    uc[2] = 0;
    uc[3] = u;
    break;
  }
  case 4:     // 0-1-3-2, v = 0 plane
  {
    uc[0] = s;
    uc[1] = u;
    uc[2] = t;
    uc[3] = 0;
    break;
  }
  default:
    avro_assert_not_reached;
  }
}

void
mapLineTraceCoordsToTriangle( const real_t* ut , int trace , real_t* uc )
{
  const real_t s = ut[0];
  switch (trace)
  {
  case 0:      // 1-2
  {
    uc[0] = 1 - s;
    uc[1] =     s;
    break;
  }
  case 1:      // 2-0
  {
    uc[0] = 0;
    uc[1] = 1 - s;
    break;
  }
  case 2:      // 0-1
  {
    uc[0] = s;
    uc[1] = 0;
    break;
  }
  default:
    avro_assert_not_reached;
  }
}

void
mapTriangleTraceCoordsToTet( const real_t* ut , int trace , real_t* uc )
{
  const real_t s = ut[0];
  const real_t t = ut[1];
  switch ( trace )
  {
    case 0:      // 1-2-3
    {
      uc[0] = 1 - s - t;
      uc[1] = s;
      uc[2] = t;
      break;
    }
    case 1:      // 0-3-2
    {
      uc[0] = 0;
      uc[1] = t;
      uc[2] = s;
      break;
    }
    case 2:      // 0-1-3
    {
      uc[0] = s;
      uc[1] = 0;
      uc[2] = t;
      break;
    }
    case 3:      // 0-2-1
    {
      uc[0] = t;
      uc[1] = s;
      uc[2] = 0;
      break;
    }
  default:
    avro_assert_not_reached;
  }
}

void
mapTraceCoordsToCanonical( coord_t number , const real_t* ut0 , int orientation , real_t* ut )
{
  switch (number)
  {
    case 1:
    {
      if (orientation < 0) ut[0] = 1 - ut0[0];
      else ut[0] = ut0[0];
      break;
    }
    case 2:
    {
      real_t s0 = ut0[0];
      real_t t0 = ut0[1];
      if (orientation < 0)
      {
        // swap s and t to correct the normal if needed
        std::swap(s0,t0);
      }
      switch ( abs(orientation) )
      {
      case 1:
      {
        ut[0] = s0;
        ut[1] = t0;
        break;
      }
      case 2:
      {
        ut[0] = 1.0 - s0 - t0;
        ut[1] = s0;
        break;
      }
      case 3:
      {
        ut[0] = t0;
        ut[1] = 1.0 - s0 - t0;
        break;
      }
      default:
        printf("unknown orientation %d\n",orientation);
        avro_assert_not_reached
      }
      break;
    }
    case 3:
    {
      mapTetTraceCoordsToCanonical( ut0[0] , ut0[1] , ut0[2] , orientation , ut[0] , ut[1] , ut[2] );
    }
  }
}

void
mapTraceCoordsToCell( coord_t number , const real_t* ut , int trace , real_t* uc )
{
  switch (number)
  {
    case 1:
      mapLineTraceCoordsToTriangle( ut , trace , uc );
      break;

    case 2:
      mapTriangleTraceCoordsToTet( ut , trace , uc );
      break;

    case 3:
      mapTetTraceCoordsToPentatope( ut , trace , uc );
  }
}

void
TraceToCellRefCoord::eval( const CanonicalTraceToCell& canonicalFace,
                           const real_t* ut0,
                           real_t* uc ) const
{
  std::vector<real_t> ut( trace_.number() );

  // map the volume (tet) coordinates to their canonical representation
  // based on the permutation decided by the orientation (sign and value)
  mapTraceCoordsToCanonical( trace_.number() , ut0 , canonicalFace.orient, ut.data() );
  mapTraceCoordsToCell( trace_.number() , ut.data() , canonicalFace.trace , uc );
}

} // avro
