#include <cmath>
#include <domain.h>
#include "pair_lj_cut_autopas.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "timer.h"
#include "autopas.h"

#include "suffix.h"

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
  auto &autopas = lmp->autopas->_autopas;

  //printf("Simulating computation in AutoPas\n");

  ev_init(eflag, vflag);



  /*
  //printf("LAMMPS -> AutoPas\n");
  // Copy to AutoPas
#ifdef AUTOPAS_OPENMP
#pragma omp parallel for default(none)
#endif
  for (int i = 0; i < atom->nlocal + atom->nghost; ++i) {

    floatVecType pos{atom->x[i][0], atom->x[i][1], atom->x[i][2]};
    floatVecType vel{atom->v[i][0], atom->v[i][1], atom->v[i][2]};
    unsigned long moleculeId = i;
    unsigned long typeId = atom->type[i];

    auto particle = ParticleType(pos, vel, moleculeId, typeId);

    if (i < atom->nlocal) {
      _autopas.addParticle(particle);
    } else {
      _autopas.addOrUpdateHaloParticle(particle);
    }


    // TODO_AUTOPAS Needs rvalue ref support in AutoPas
    // _autopas.addParticle(ParticleType(pos, vel, moleculeId, typeId))
  }*/


  //printf("AutoPas computation\n");


  // Force calculation
  //PairFunctorType functor{_autopas.getCutoff(), *_particlePropertiesLibrary};
  AutoPasLMP::PairFunctorType functor{autopas->getCutoff()};
  functor.setParticleProperties(24, 1);

  autopas->iteratePairwise(&functor);

  timer->stamp(Timer::AUTOPAS);

  /*
  //printf("AutoPas -> LAMMPS\n");
  // Copy from AutoPas
#ifdef AUTOPAS_OPENMP
#pragma omp parallel default(none)
#endif
  for (auto iter = _autopas.begin(
      autopas::IteratorBehavior::ownedOnly); iter.isValid(); ++iter) {
    auto force = iter->getF();
    int moleculeId = iter->getID();

    atom->f[moleculeId][0] += force[0];
    atom->f[moleculeId][1] += force[1];
    atom->f[moleculeId][2] += force[2];
  }

   */

  std::copy_n(functor.getVirial()->begin(), 6, virial);
  auto upot = functor.getUpot();
  eng_vdwl = upot;

  //printf("AutoPas complete\n");
 // _autopas.deleteAllParticles();

}

/* ---------------------------------------------------------------------- */

double PairLJCutAutoPas::memory_usage() {
  double bytes = 0; // memory_usage_thr(); //TODO_AUTOPAS Get memory usage from AutoPas
  bytes += PairLJCut::memory_usage();

  return bytes;
}

void PairLJCutAutoPas::init_style() {
  PairLJCut::init_style();
  lmp->autopas->init_autopas(cut_global, epsilon, sigma);
}
