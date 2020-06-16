#include "pair_lj_cut_autopas.h"

#include "atom.h"
#include "autopas.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "suffix.h"
#include "timer.h"

#include <domain.h>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairLJCutAutoPas::PairLJCutAutoPas(LAMMPS *lmp) :
    PairLJCut(lmp) {
  suffix_flag |= Suffix::AUTOPAS;
  respa_enable = 0;
  cut_respa = nullptr;
}

/* ---------------------------------------------------------------------- */

void PairLJCutAutoPas::compute(int eflag, int vflag) {

  ev_init(eflag, vflag);

  // Force calculation
  //PairFunctorType functor{_autopas.getCutoff(), *_particlePropertiesLibrary};
  AutoPasLMP::PairFunctorType functor{lmp->autopas->get_cutoff(), lmp->autopas->get_interaction_map()};
  functor.setParticleProperties(24 * epsilon[1][1], sigma[1][1] * sigma[1][1]);

  lmp->autopas->iterate_pairwise(&functor);

  timer->stamp(Timer::AUTOPAS);

  // Copy virial
  std::copy_n(functor.getVirial()->begin(), 6, virial);
  auto upot = functor.getUpot();
  eng_vdwl = upot;

}

/* ---------------------------------------------------------------------- */

double PairLJCutAutoPas::memory_usage() {
  double bytes = 0; //TODO Any additional memory usage for pair calculation?
  bytes += PairLJCut::memory_usage();

  return bytes;
}

void PairLJCutAutoPas::init_style() {
  PairLJCut::init_style();
  lmp->autopas->init_autopas(cut_global, epsilon, sigma);
}
