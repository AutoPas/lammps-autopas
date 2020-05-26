#ifdef FIX_CLASS

FixStyle(temp/rescale/autopas,FixTempRescaleAutoPas)

#else

#ifndef LMP_FIX_TEMP_RESCALE_AUTOPAS_H
#define LMP_FIX_TEMP_RESCALE_AUTOPAS_H

#include "fix_temp_rescale.h"

namespace LAMMPS_NS {

class FixTempRescaleAutoPas : public FixTempRescale {
public:
  using FixTempRescale::FixTempRescale;

  void end_of_step() override;

};
};

#endif
#endif
