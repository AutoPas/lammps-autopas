#ifdef FIX_CLASS

FixStyle(enforce2d/autopas,FixEnforce2DAutoPas)

#else

#ifndef LMP_FIX_ENFORCE2D_AUTOPAS_H
#define LMP_FIX_ENFORCE2D_AUTOPAS_H

#include "fix_enforce2d.h"

namespace LAMMPS_NS {

class FixEnforce2DAutoPas : public FixEnforce2D {
public:
  using FixEnforce2D::FixEnforce2D;

  void post_force(int) override;
};

}

#endif
#endif
