namespace ursa
{

namespace geometrics
{

template<typename T>
void
Primitive<T>::project( numerics::Coordinate& x ) const
{
  prim_.project(x);
}

} // geometrics

} // ursa
