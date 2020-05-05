#ifdef COMPUTE_CLASS

ComputeStyle(temp/autopas,ComputeTempAutoPas)

#else

#ifndef LMP_COMPUTE_TEMP_AUTOPAS_H
#define LMP_COMPUTE_TEMP_AUTOPAS_H

#include "compute_temp.h"

namespace LAMMPS_NS {

class ComputeTempAutoPas : public ComputeTemp {
public:
  using ComputeTemp::ComputeTemp;

  double compute_scalar() override;

  void compute_vector() override;

};

}
#endif
#endif
