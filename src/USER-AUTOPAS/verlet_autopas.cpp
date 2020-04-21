#include "verlet_autopas.h"

#include "neighbor.h"
#include "domain.h"
#include "autopas.h"
#include "comm.h"
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "output.h"
#include "update.h"
#include "modify.h"
#include "timer.h"
#include "error.h"

using namespace LAMMPS_NS;

VerletAutoPas::VerletAutoPas(LAMMPS *lmp, int narg, char **arg) :
    Verlet(lmp, narg, arg) {

}

void VerletAutoPas::setup(int flag) {
  Verlet::setup(flag);
}

void VerletAutoPas::setup_minimal(int i) {
  Verlet::setup_minimal(i);
}

void VerletAutoPas::run(int i) {
  Verlet::run(i);
}

void VerletAutoPas::force_clear() {

#pragma omp parallel default(none)
  for (auto iter = lmp->autopas->_autopas->begin(
      autopas::IteratorBehavior::ownedOnly); iter.isValid(); ++iter) {
    iter->setF({0, 0, 0});
  }
}
