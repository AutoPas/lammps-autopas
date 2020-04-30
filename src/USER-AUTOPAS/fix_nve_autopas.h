#ifdef FIX_CLASS

FixStyle(nve/autopas,FixNVEAutoPas)

#else

#ifndef LMP_FIX_NVE_AUTOPAS_H
#define LMP_FIX_NVE_AUTOPAS_H

#include "fix_nve.h"

namespace LAMMPS_NS {

class FixNVEAutoPas : public FixNVE {

public:
  FixNVEAutoPas(class LAMMPS *, int, char **);

  void initial_integrate(int) override;

  void final_integrate() override;

protected:
  template<bool initial>
  void do_integrate();
};

}
#endif
#endif
