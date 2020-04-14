#ifndef AUTOPAS_LMP_H
#define AUTOPAS_LMP_H

#include "pointers.h"
// #include "autopas_type.h" // TODO_AP: Same types as LAMMPS
// #include "pair_autopas.h"  // TODO_AP: No default pair style for AP

namespace LAMMPS_NS {

class AutoPasLMP : protected Pointers {
public:
  int autopas_exists;

  AutoPasLMP(class LAMMPS *, int, char **);
  ~AutoPasLMP();
  void accelerator(int, char **);
  int neigh_count(int);

};

}

#endif