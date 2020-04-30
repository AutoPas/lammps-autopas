#ifndef LMP_COMM_BRICK_AUTOPAS_H
#define LMP_COMM_BRICK_AUTOPAS_H

#include "comm_brick.h"

#include <vector>

#include "autopas.h"

namespace LAMMPS_NS {

class CommAutoPas : public CommBrick {
public:
  explicit CommAutoPas(class LAMMPS *);

  void forward_comm(int) override;

  void reverse_comm() override;

  void exchange() override;

  void borders() override;

  void setup() override;

private:

  template<bool haloOnly = false>
  void border_impl(int idxfirst, int idxlast, double lo, double hi, int dim,
                   std::vector<AutoPasLMP::ParticleType *> &sendparticles) const;

  template<bool haloOnly = false>
  void border_impl(int idxfirst, int idxlast, double *mlo, double *mhi, int dim,
                   std::vector<AutoPasLMP::ParticleType *> &sendparticles) const;


  void reverse_comm_impl_self(int iswap) const;

  void reverse_comm_impl_other(int iswap) const;

  std::vector<std::vector<AutoPasLMP::ParticleType *>> _sendlist_particles;
};
}

#endif
