#pragma once

#include "domain.h"

namespace LAMMPS_NS {

class DomainAutoPas : public Domain {
public:
  explicit DomainAutoPas(class LAMMPS *);

  ~DomainAutoPas() override = default;

  void pbc() override;

  void image_check() override;

  void box_too_small_check() override;

  int closest_image(int i, int i1) override;

  int closest_image(const double *const pDouble, int i) override;

  void lamda2x(int i) override;

  void x2lamda(int i) override;

};

}