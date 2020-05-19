#ifndef LMP_DOMAIN_AUTOPAS_H
#define LMP_DOMAIN_AUTOPAS_H

#include "domain.h"

#include "autopas.h"

namespace LAMMPS_NS {

class DomainAutoPas : public Domain {
public:
  using Domain::Domain;

  void pbc() override;

  void image_check() override;

  void box_too_small_check() override;

  int closest_image(int, int) override;

  int closest_image(const double *, int) override;

  void lamda2x(int) override;

  void x2lamda(int) override;

  void reset_box() override;

protected:
  virtual bool
  pbc(AutoPasLMP::ParticleType &, double *lo, double *hi, double *period);

};

}

#endif
