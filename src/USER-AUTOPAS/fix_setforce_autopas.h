#ifdef FIX_CLASS

FixStyle(setforce/autopas,FixSetForceAutoPas)

#else

#ifndef LMP_FIX_SET_FORCE_AUTOPAS_H
#define LMP_FIX_SET_FORCE_AUTOPAS_H

#include "fix_setforce.h"

namespace LAMMPS_NS {

class FixSetForceAutoPas : public FixSetForce {
public:
  using FixSetForce::FixSetForce;

  void post_force(int) override;

  void post_force_respa(int, int, int) override;

};

}

#endif
#endif
