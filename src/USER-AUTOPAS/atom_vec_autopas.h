#ifndef LMP_ATOM_VEC_AUTOPAS_H
#define LMP_ATOM_VEC_AUTOPAS_H

#include "atom_vec.h"
#include "autopas.h"

#include <vector>

namespace LAMMPS_NS {

class AtomVecAutopas : public AtomVec {

public:
  explicit AtomVecAutopas(class LAMMPS *);

  int pack_comm(int, int *, double *, int, int *) override;

  int pack_comm_vel(int, int *, double *, int, int *) override;

  void unpack_comm(int, int, double *) override;

  void unpack_comm_vel(int, int, double *) override;

  int pack_reverse(int, int, double *) override;

  void unpack_reverse(int, int *, double *) override;

  virtual int pack_exchange(const AutoPasLMP::ParticleType &, double *) = 0;

  /*
   * AutoPas versions
   */

  virtual int
  pack_border_autopas(const std::vector<AutoPasLMP::ParticleType *> &, double *,
                      int, const int *) = 0;

  virtual int
  pack_border_vel_autopas(const std::vector<AutoPasLMP::ParticleType *> &,
                          double *,
                          int, const int *) = 0;

  void unpack_reverse_autopas(
      int sendnum,
      const int* sendlist,
      const std::vector<std::array<double, 3>>& force_buf, int firstrecv);

protected:
  int *mask{};

};
}

#endif
