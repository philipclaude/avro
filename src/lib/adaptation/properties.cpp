#include "adaptation/metric.h"
#include "adaptation/properties.h"

#include "master/simplex.h"

#include "mesh/topology.h"

#include <json/json.hpp>

#include <fstream>

namespace avro
{

template<typename type>
Properties::Properties( const Topology<type>& topology ,
                        MetricField<type>& metric ) :
  lstats_(length_),
  qstats_(quality_)
{
  // compute the length and quality values, then bin with the defaults
  compute(topology,metric);
  bin();
}

template<typename type>
void
Properties::compute( const Topology<type>& topology ,
                     MetricField<type>& metric )
{
  // evaluate the edge lengths
  metric.lengths( topology , length_ );

  // evaluate the cell quality
  quality_.clear();
  for (index_t k=0;k<topology.nb();k++)
  {
    if (topology.ghost(k)) continue;
    quality_.push_back( metric.quality( topology , k ) );
  }

  std::sort( length_.begin() , length_.end() );
  std::sort( quality_.begin() , quality_.end() );
  bin();
}

void
Properties::bin( const std::vector<real_t>& llims ,
                 const std::vector<real_t>& qlims )
{

  if (llims.size()==0)
  {
    real_t lmin = lstats_.min();
    real_t lmax = lstats_.max();
    llims_ = { lmin , 0.5, sqrt(0.5) , 1. , sqrt(2.) , 2. ,lmax };
  }
  else
    llims_ = llims;

  if (qlims.size()==0)
  {
    real_t qmin = qstats_.min();
    real_t qmax = qstats_.max();
    qlims_ = { qmin , 0.4 , 0.8 , qmax };
  }
  else
    qlims_ = qlims;

  lbin_.resize( llims_.size() -1 );
  qbin_.resize( qlims_.size() -1 );

  // compute the bins
  for (index_t j=1;j<llims_.size();j++)
    lbin_[j-1] = lstats_.count( llims_[j-1] , llims_[j] );
  for (index_t j=1;j<qlims_.size();j++)
    qbin_[j-1] = qstats_.count( qlims_[j-1] , qlims_[j] );
}

void
tab( index_t nt )
{
  for (index_t j=0;j<nt;j++)
    printf("\t");
}

void
Properties::print( const std::string& title , index_t nt ) const
{
  index_t ltotal = 0;
  for (index_t j=0;j<lbin_.size();j++)
    ltotal += lbin_[j];
  index_t qtotal = 0;
  for (index_t j=0;j<qbin_.size();j++)
    qtotal += qbin_[j];

  if (ltotal==0) ltotal = 1;
  if (qtotal==0) qtotal = 1;

  tab(nt); printf("%s:\n",title.c_str());
  tab(nt); printf("\tlengths: [%3.2f,%3.2f], avg = %3.2f for %lu edges\n",
                    lstats_.min(),lstats_.max(),lstats_.avg(),lstats_.nb());
  for (index_t j=0;j<lbin_.size();j++)
    tab(nt),printf("\t\t%3.4f < l < %3.4f = %6lu (%3.2f %%)\n",
            llims_[j],llims_[j+1],lbin_[j],lbin_[j]*100./ltotal);

  tab(nt); printf("\tquality: [%3.2f,%3.2f], avg = %3.2f for %lu elems\n",
                    qstats_.min(),qstats_.max(),qstats_.avg(),qstats_.nb());
  for (index_t j=0;j<qbin_.size();j++)
    printf("\t\t%3.4f < q < %3.4f = %6lu (%3.2f %%)\n",
            qlims_[j],qlims_[j+1],qbin_[j],qbin_[j]*100./qtotal);
}

void
Properties::dump( const std::string& filename ) const
{
  nlohmann::json jfile;
  nlohmann::json jlength;
  nlohmann::json jquality;

  jlength["values"] = length_;
  jlength["min"] = lstats_.min();
  jlength["max"] = lstats_.max();
  jlength["avg"] = lstats_.avg();
  jlength["stdev"] = lstats_.stdev();
  jlength["nb"] = lstats_.nb();

  jquality["values"] = quality_;
  jquality["min"] = qstats_.min();
  jquality["max"] = qstats_.max();
  jquality["avg"] = qstats_.avg();
  jquality["stdev"] = qstats_.stdev();
  jquality["nb"] = qstats_.nb();

  jfile["length"] = jlength;
  jfile["quality"] = jquality;

  std::ofstream file(filename);
  file << std::setw(2) << jfile << std::endl;
}


// instantiate the properties for simplices
template Properties::Properties( const Topology<Simplex>& ,
                                 MetricField<Simplex>& );
template void Properties::compute( const Topology<Simplex>& , MetricField<Simplex>& );

} // avro
