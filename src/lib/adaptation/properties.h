#ifndef AVRO_MESH_LOCAL_PROPERTIES_H_
#define AVRO_MESH_LOCAL_PROPERTIES_H_

#include <string>
#include <vector>

namespace avro
{

template<typename type> class Topology;

namespace numerics
{
template<typename type> class MetricField;
}
class Distribution
{
public:
  Distribution( const std::vector<real_t>& data ) :
    data_(data)
  {}

  index_t nb() const { return data_.size(); }

  real_t min() const { return *std::min_element(data_.begin(),data_.end()); }
  real_t max() const { return *std::max_element(data_.begin(),data_.end()); }
  real_t avg() const
  {
    real_t A = 0.0;
    for (index_t k=0;k<nb();k++)
      A += data_[k];
    return A/nb();
  }

  real_t stdev() const
  {
    real_t S = 0.0;
    real_t A = avg();
    for (index_t k=0;k<nb();k++)
      S += (data_[k] -A)*(data_[k] -A);
    return sqrt(S/nb());
  }

  index_t count( real_t dlo , real_t dhi ) const
  {
    index_t n = 0;
    for (index_t k=0;k<nb();k++)
    {
      if (data_[k]>dlo && data_[k]<dhi)
        n++;
    }
    return n;
  }

private:
  const std::vector<real_t>& data_;
};

class Properties
{

public:
  template<typename type>
  Properties( const Topology<type>& topology ,
              MetricField<type>& metric );

  template<typename type>
  void compute( const Topology<type>& topology ,
                MetricField<type>& metric );

  void bin( const std::vector<real_t>& llims=std::vector<real_t>() ,
            const std::vector<real_t>& qlims=std::vector<real_t>() );
  void print( const std::string& title="metric properties" ,
              const index_t nt=0 ) const;
  void dump( const std::string& filename ) const;

  real_t qavg() const { return qstats_.avg(); }

  void conformity( real_t& lunit , real_t& qunit , index_t& nb_elem ) const;

private:
  Distribution lstats_;
  std::vector<real_t> length_;
  std::vector<index_t> lbin_;
  std::vector<real_t> llims_;

  Distribution qstats_;
  std::vector<real_t> quality_;
  std::vector<index_t> qbin_;
  std::vector<real_t> qlims_;
};

} // avro

#endif
