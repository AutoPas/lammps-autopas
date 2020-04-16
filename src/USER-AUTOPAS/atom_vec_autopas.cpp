#include "atom_vec_autopas.h"
// #include <iterator>     // std::next
#include "domain.h"

using namespace LAMMPS_NS;

AtomVecAutopas::AtomVecAutopas(LAMMPS_NS::LAMMPS *lmp)
    : AtomVec(lmp) {

}

AutoPasLMP::FloatVecType LAMMPS_NS::AtomVecAutopas::h_x(int i) {
  auto it = this->lmp->autopas->_autopas->begin();
  //std::next(it, 3)
  return {0,0,0};
}



int
LAMMPS_NS::AtomVecAutopas::pack_comm(int n, int *list, double *buf,
                                           int pbc_flag, int *pbc) {
  throw "Not implemented";
return 0;
}

int LAMMPS_NS::AtomVecAutopas::pack_comm_vel(int, int *, double *, int,
                                                   int *) {
  throw "Not implemented";
  return 0;
}

void LAMMPS_NS::AtomVecAutopas::unpack_comm(int, int, double *) {

  throw "Not implemented";
}

void LAMMPS_NS::AtomVecAutopas::unpack_comm_vel(int, int, double *) {
  throw "Not implemented";

}

int LAMMPS_NS::AtomVecAutopas::pack_reverse(int, int, double *) {
  throw "Not implemented";
  return 0;
}

void LAMMPS_NS::AtomVecAutopas::unpack_reverse(int, int *, double *) {
  throw "Not implemented";

}
