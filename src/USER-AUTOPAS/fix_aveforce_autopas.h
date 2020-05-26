#ifdef FIX_CLASS

FixStyle(aveforce/autopas,FixAveForceAutoPas)

#else

#ifndef LMP_FIX_AVEFORCE_AUTOPAS_H
#define LMP_FIX_AVEFORCE_AUTOPAS_H

#include "fix_aveforce.h"

namespace LAMMPS_NS {

class FixAveForceAutoPas : public FixAveForce {
public:
  using FixAveForce::FixAveForce;

  void post_force(int) override;

  void post_force_respa(int, int, int) override;

};

}

#endif
#endif
