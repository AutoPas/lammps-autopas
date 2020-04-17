#include "atom_autopas.h"
#include "error.h"

using namespace LAMMPS_NS;

AtomAutoPas::AtomAutoPas(class LAMMPS *lmp) : Atom(lmp) {

}

AtomAutoPas::~AtomAutoPas() {

}

void AtomAutoPas::sort() {
  sortfreq = 0;
  error->warning(FLERR,"Atom sorting is not supported using AutoPas "
                       "and will thus be disabled.");
}
