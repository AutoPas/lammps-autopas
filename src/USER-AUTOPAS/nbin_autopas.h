#ifdef NBIN_CLASS

NBinStyle(autopas,NBinAutoPas,0)

#else

#ifndef LMP_NBIN_AUTOPAS_H
#define LMP_NBIN_AUTOPAS_H

#include "nbin_standard.h"

namespace LAMMPS_NS {

class NBinAutoPas : public NBinStandard {
public:

  explicit NBinAutoPas(class LAMMPS *);
  ~NBinAutoPas() {}

  void bin_atoms() override;
};

}

#endif
#endif
