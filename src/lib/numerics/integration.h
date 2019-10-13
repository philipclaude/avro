#ifndef URSA_LIB_NUMERICS_INTEGRATION_H_
#define URSA_LIB_NUMERICS_INTEGRATION_H_

namespace ursa
{

template<typename> class Topology;

template<typename Function_t>
struct Integrand
{
  template<typename T>
  void operator() ( const real_t* X , T& f ) const
  {
    cast()(X,f);
  }

  Function_t& cast() { return static_cast<Function_t>(*this); }
  const Function_t& cast() const { return static_cast<const Function_t>(*this); }
};

template<typename Basis_t,typename Function_t>
class Integral
{
  Integral( const Basis_t& basis );

  template<typename T>
  void integrate( const Topology<Basis_t>& topology , index_t k , T& f ) const
  {
    T df;
    f = 0;
    for (index_t j=0;j<master_.nb_quad();j++)
    {
      integrand_( master_.x(k) , df );
      f += master_.w(k)*df*master_.jacobian( topology.vertices() , topology(k) , topology.nv(k) );
    }
  }

private:
  const Basis_t& master_;

};

} // ursa

#endif
