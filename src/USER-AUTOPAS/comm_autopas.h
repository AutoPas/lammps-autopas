#pragma once

#include "comm_brick.h"

namespace LAMMPS_NS {

class CommAutoPas : public CommBrick {
public:
  explicit CommAutoPas(class LAMMPS *);
  ~CommAutoPas() override = default;

  void forward_comm(int) override;

  void reverse_comm() override;

  void exchange() override;

  void borders() override;

private:

  template<bool haloOnly=false>
  void border_impl(int idxfirst, int idxlast, double lo, double hi, int dim, std::vector<AutoPasLMP::ParticleType *> &sendparticles) const;
  template<bool haloOnly=false>
  void border_impl(int idxfirst, int idxlast, double *mlo, double *mhi, int dim, std::vector<AutoPasLMP::ParticleType *> &sendparticles) const;

};
}
