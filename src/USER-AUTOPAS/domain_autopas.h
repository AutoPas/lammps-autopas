#pragma once

#include "domain.h"
#include "autopas.h"

namespace LAMMPS_NS {

class DomainAutoPas : public Domain {
public:
  explicit DomainAutoPas(class LAMMPS *);

  ~DomainAutoPas() override = default;

  void pbc() override;

  void image_check() override;

  void box_too_small_check() override;

  int closest_image(int, int) override;

  int closest_image(const double *const, int) override;

  void lamda2x(int) override;

  void x2lamda(int) override;

  /*
   * AutoPas only functions
   */

protected:
  virtual bool pbc(AutoPasLMP::ParticleType &, double *lo, double *hi, double *period);

};

}