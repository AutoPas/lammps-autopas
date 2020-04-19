
#ifdef INTEGRATE_CLASS

IntegrateStyle(verlet/autopas,VerletAutoPas)

#else

#ifndef LMP_VERLET_AUTOPAS_H
#define LMP_VERLET_AUTOPAS_H

#include "verlet.h"

namespace LAMMPS_NS {

class VerletAutoPas : public Verlet {
public:
  explicit VerletAutoPas(class LAMMPS *, int, char **);

  ~VerletAutoPas() override = default;

  void setup(int flag) override;

  void setup_minimal(int i) override;

  void run(int i) override;

protected:
  void force_clear() override;

};


}
#endif
#endif