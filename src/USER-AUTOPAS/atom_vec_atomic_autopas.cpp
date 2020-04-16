#include "atom_vec_atomic_autopas.h"

LAMMPS_NS::AtomVecAtomicAutopas::AtomVecAtomicAutopas(LAMMPS_NS::LAMMPS *lmp)
    : AtomVecAutopas(lmp) {

}

void LAMMPS_NS::AtomVecAtomicAutopas::grow(int i) {
  throw "Not implemented";

}
void LAMMPS_NS::AtomVecAtomicAutopas::grow_reset() {
  // This is only called from special.cpp line 718
  throw "Not implemented";
}


void LAMMPS_NS::AtomVecAtomicAutopas::copy(int i, int i1, int i2) {
  throw "Not implemented";

}

int
LAMMPS_NS::AtomVecAtomicAutopas::pack_border(int i, int *pInt, double *pDouble,
                                             int i1, int *pInt1) {
  throw "Not implemented";
  return 0;
}

int LAMMPS_NS::AtomVecAtomicAutopas::pack_border_vel(int i, int *pInt,
                                                     double *pDouble, int i1,
                                                     int *pInt1) {
  throw "Not implemented";
  return 0;
}

void
LAMMPS_NS::AtomVecAtomicAutopas::unpack_border(int i, int i1, double *pDouble) {
  throw "Not implemented";

}

void LAMMPS_NS::AtomVecAtomicAutopas::unpack_border_vel(int i, int i1,
                                                        double *pDouble) {
  throw "Not implemented";

}

int LAMMPS_NS::AtomVecAtomicAutopas::pack_exchange(int i, double *pDouble) {
  throw "Not implemented";
  return 0;
}

int LAMMPS_NS::AtomVecAtomicAutopas::unpack_exchange(double *pDouble) {
  throw "Not implemented";
  return 0;
}

int LAMMPS_NS::AtomVecAtomicAutopas::size_restart() {
  throw "Not implemented";
  return 0;
}

int LAMMPS_NS::AtomVecAtomicAutopas::pack_restart(int i, double *pDouble) {
  throw "Not implemented";
  return 0;
}

int LAMMPS_NS::AtomVecAtomicAutopas::unpack_restart(double *pDouble) {
  throw "Not implemented";
  return 0;
}

void LAMMPS_NS::AtomVecAtomicAutopas::create_atom(int i, double *pDouble) {
  throw "Not implemented";

}

void LAMMPS_NS::AtomVecAtomicAutopas::data_atom(double *pDouble,
                                                LAMMPS_NS::imageint imageint1,
                                                char **pString) {
  throw "Not implemented";

}

void LAMMPS_NS::AtomVecAtomicAutopas::pack_data(double **pDouble) {
  throw "Not implemented";

}

void LAMMPS_NS::AtomVecAtomicAutopas::write_data(FILE *file, int i,
                                                 double **pDouble) {
  throw "Not implemented";

}

LAMMPS_NS::bigint LAMMPS_NS::AtomVecAtomicAutopas::memory_usage() {
  throw "Not implemented";
  return 0;
}
