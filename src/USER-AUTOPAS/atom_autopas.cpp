#include "atom_autopas.h"

#include "error.h"

using namespace LAMMPS_NS;

void AtomAutoPas::sort() {
  sortfreq = 0;
  error->warning(FLERR, "Atom sorting is not supported using AutoPas "
                        "and will be disabled.");
}
