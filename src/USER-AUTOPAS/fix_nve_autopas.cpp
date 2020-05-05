#include "fix_nve_autopas.h"

#include "atom.h"
#include "autopas.h"

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

  // update v and x of atoms in group

  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;

#pragma omp parallel default(none) shared(rmass, mass, type, mask)
  for (auto iter = lmp->autopas->iterate<autopas::ownedOnly>(); iter.isValid(); ++iter) {
    auto &particle = *iter;
    auto idx{lmp->autopas->particle_to_index(particle)};
    if (mask[idx] & groupbit) {

      const double dtfm = dtf / (rmass ? rmass[idx] : mass[type[idx]]);

      particle.addV(mulScalar(particle.getF(), dtfm));

      if constexpr (initial) {
        particle.addR(mulScalar(particle.getV(), dtv));
      }
    }
  }
}
