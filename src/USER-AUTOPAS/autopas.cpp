#include "autopas.h"

using namespace LAMMPS_NS;

AutoPasLMP::AutoPasLMP(class LAMMPS *lmp, int narg, char **args) : Pointers(lmp) {
  autopas_exists = 1;
  lmp->autopas = this;
}

AutoPasLMP::~AutoPasLMP() {

}

int AutoPasLMP::neigh_count(int) {
  return 0;
}

void AutoPasLMP::accelerator(int, char **) {

}