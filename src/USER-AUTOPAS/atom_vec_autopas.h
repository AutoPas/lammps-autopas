#pragma once

#include "atom_vec.h"
#include "autopas.h"

namespace LAMMPS_NS {

class AtomVecAutopas : public AtomVec {

public:
  explicit AtomVecAutopas(class LAMMPS *);
  ~AtomVecAutopas() override = default;

  int pack_comm(int, int *, double *, int, int *) override;
  int pack_comm_vel(int, int *, double *, int, int *) override;
  void unpack_comm(int, int, double *) override;
  void unpack_comm_vel(int, int, double *) override;
  int pack_reverse(int, int, double *) override;
  void unpack_reverse(int, int *, double *) override;

  virtual int pack_exchange(const AutoPasLMP::ParticleType &, double *) = 0;

protected:
  int *mask;

};
}