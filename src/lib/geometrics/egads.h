#ifndef URSA_LIB_LIBRARY_EGADS_H_
#define URSA_LIB_LIBRARY_EGADS_H_

struct egObject;
typedef egObject* ego;

namespace ursa
{

namespace numerics
{
  class Coordinate;
}

namespace geometrics
{

namespace EGADS
{

class Object
{
public:
  Object();
  Object( ego* );

  void project( numerics::Coordinate& x ) const;

  ego* object();
  const ego* object() const;

private:
  ego* ego_;
};

} // EGADS

} // geometrics

} // ursa

#endif
