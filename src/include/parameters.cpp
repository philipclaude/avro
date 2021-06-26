#include "parameters.h"

namespace avro
{

ParameterSet::ParameterSet() {
  register_parameter( "output_redirect" , std::string() , "defines where avro output will be redirected (if supported)" );
  register_parameter( "directory" , std::string("./") , "defines where output files will be written" );
  register_parameter( "prefix" , std::string("avro_") , "output prefix" );
  register_parameter( "write_mesh" , true , "controls whether the mesh is written upon completion of the algortihm" );
  register_parameter( "limit_metric" , false , "controls whether the metric is limited from the current mesh" );
  register_parameter( "metric limiting factor" , 2.0 , "factor by which entries of the step matrix are limited when limiting the metric" );
  register_parameter( "nb_iter_cvt" , 10 , "controls the number of iterations performed when computing a Centroidal Voronoi Tesselation" );
  register_parameter( "nb_iter_sdot" , 10 , "controls the number of iterations performed when optimizing the weights of a power diagram" );
  register_parameter( "max_passes" , 3 , "controls the maximum number of passes used in a parallel adaptation run" );
  register_parameter( "is_curved" , true , "controls whether an input geometry is curved (which will make for a faster adaptation if false)" );
}

} // avro
