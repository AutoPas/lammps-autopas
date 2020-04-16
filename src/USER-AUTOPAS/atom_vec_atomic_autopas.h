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
  ~AtomVecAtomicAutopas() override = default;

  void grow(int i) override;

  void grow_reset() override;

  void copy(int i, int i1, int i2) override;

  int
  pack_border(int i, int *pInt, double *pDouble, int i1, int *pInt1) override;

  int pack_border_vel(int i, int *pInt, double *pDouble, int i1,
                      int *pInt1) override;

  void unpack_border(int i, int i1, double *pDouble) override;

  void unpack_border_vel(int i, int i1, double *pDouble) override;

  int pack_exchange(int i, double *pDouble) override;

  int unpack_exchange(double *pDouble) override;

  int size_restart() override;

  int pack_restart(int i, double *pDouble) override;

  int unpack_restart(double *pDouble) override;

  void create_atom(int i, double *pDouble) override;

  void data_atom(double *pDouble, imageint imageint1, char **pString) override;

  void pack_data(double **pDouble) override;

  void write_data(FILE *file, int i, double **pDouble) override;

  bigint memory_usage() override;

protected:
  tagint *tag;
  int *type,*mask;
  imageint *image;
};
}

#endif
#endif