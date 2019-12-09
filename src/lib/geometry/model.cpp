#include "geometry/model.h"

namespace luna
{

void
Model::determine_number()
{
  number_ = 0;
  for (index_t k=0;k<nb_bodies();k++)
  {
    if (body_[k]->number()>number_)
      number_ = body_[k]->number();
  }
}

void
Model::print() const
{
  for (index_t k=0;k<nb_bodies();k++)
    body_[k]->print();
}

} // luna
