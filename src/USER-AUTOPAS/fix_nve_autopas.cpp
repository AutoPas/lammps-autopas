#include "fix_nve_autopas.h"
#include "atom.h"
#include "autopas.h"
#include <autopas/utils/ArrayMath.h>


LAMMPS_NS::FixNVEAutoPas::FixNVEAutoPas(LAMMPS *lmp, int narg, char **arg) :
    FixNVE(lmp, narg, arg) {

}

void LAMMPS_NS::FixNVEAutoPas::initial_integrate(int /*vflag*/) {
  do_integrate</*initial*/ true>();
}

void LAMMPS_NS::FixNVEAutoPas::final_integrate() {
  do_integrate</*initial*/ false>();
}

template<bool initial>
void LAMMPS_NS::FixNVEAutoPas::do_integrate() {
  using autopas::utils::ArrayMath::mulScalar;

  double dtfm;

  // update v and x of atoms in group

  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  auto &autopas = *lmp->autopas->_autopas;
  // int nlocal = atom->nlocal;
  // if (igroup == atom->firstgroup) nlocal = atom->nfirst;

#pragma omp parallel default(none) shared(autopas, rmass, mass, type, mask) private(dtfm)
  for (auto iter = autopas.begin(
      autopas::IteratorBehavior::ownedOnly); iter.isValid(); ++iter) {
    auto &particle = *iter;
    int idx = particle.getID(); // TODO Handle global index to local mapping?
    if (mask[idx] & groupbit) {

      if (rmass) {
        dtfm = dtf / rmass[idx];
      } else {
        dtfm = dtf / mass[type[idx]];
      }

      particle.addV(mulScalar(particle.getF(), dtfm));

      if constexpr (initial) {
        particle.addR(mulScalar(particle.getV(), dtv));
      }
    }
  }
}

