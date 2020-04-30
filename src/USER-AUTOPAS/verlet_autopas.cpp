#include "verlet_autopas.h"

#include "angle.h"
#include "atom.h"
#include "atom_vec.h"
#include "autopas.h"
#include "bond.h"
#include "comm.h"
#include "dihedral.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "improper.h"
#include "kspace.h"
#include "modify.h"
#include "neighbor.h"
#include "output.h"
#include "pair.h"
#include "timer.h"
#include "update.h"

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
  for (auto iter = lmp->autopas->iterate<autopas::IteratorBehavior::ownedOnly>(); iter.isValid(); ++iter) {
    iter->setF({0, 0, 0});
  }
}
