#ifdef ATOM_CLASS

AtomStyle(atomic/autopas,AtomVecAtomicAutopas)

#else

#ifndef LMP_ATOM_VEC_ATOMIC_AUTOPAS_H
#define LMP_ATOM_VEC_ATOMIC_AUTOPAS_H

#include "atom_vec_autopas.h"

namespace LAMMPS_NS {

class AtomVecAtomicAutopas : public AtomVecAutopas {

public:
  explicit AtomVecAtomicAutopas(class LAMMPS *);

  void grow(int i) override;

  void grow_reset() override;

  void copy(int, int, int) override;

  int pack_border(int, int *, double *, int, int *) override;

  int pack_border_vel(int, int *, double *, int, int *) override;

  void unpack_border(int, int, double *) override;

  void unpack_border_vel(int, int, double *) override;

  int pack_exchange(int, double *) override;

  int unpack_exchange(double *) override;

  int size_restart() override;

  int pack_restart(int, double *) override;

  int unpack_restart(double *) override;

  void create_atom(int, double *) override;

  void data_atom(double *, imageint, char **) override;

  void pack_data(double **) override;

  void write_data(FILE *, int, double **) override;

  bigint memory_usage() override;

  /*
   * AutoPas variants
   */

  int pack_exchange(const AutoPasLMP::ParticleType &, double *) override;


  int
  pack_border_autopas(const std::vector<AutoPasLMP::ParticleType *> &, double *,
                      int, const int *) override;

  int pack_border_vel_autopas(const std::vector<AutoPasLMP::ParticleType *> &,
                              double *,
                              int, const int *) override;

protected:
  tagint *tag;
  int *type;
  imageint *image;
  double **x, **v;//,**f;
};
}

#endif
#endif
