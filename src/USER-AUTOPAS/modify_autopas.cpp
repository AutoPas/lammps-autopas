#include "modify_autopas.h"
#include "atom_autopas.h"
#include "update.h"
#include "fix.h"
#include "compute.h"
#include "autopas.h"

using namespace LAMMPS_NS;

ModifyAutoPas::ModifyAutoPas(LAMMPS *lmp) : Modify(lmp)
{
  // atomKK = (AtomKokkos *) atom;
}