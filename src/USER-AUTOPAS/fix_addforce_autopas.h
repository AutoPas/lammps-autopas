#ifdef FIX_CLASS

FixStyle(addforce/autopas,FixAddForceAutoPas)

#else

#ifndef LMP_FIX_ADDFORCE_AUTOPAS_H
#define LMP_FIX_ADDFORCE_AUTOPAS_H

#include "fix_addforce.h"

namespace LAMMPS_NS {

class FixAddForceAutoPas : public FixAddForce {
public:
  using FixAddForce::FixAddForce;

  void post_force(int) override;
};
};

#endif
#endif
