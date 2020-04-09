/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include <cmath>
#include <domain.h>
#include "pair_lj_cut_autopas.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "timer.h"

#include "suffix.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairLJCutAutoPas::PairLJCutAutoPas(LAMMPS *lmp) :
    PairLJCut(lmp) {
  suffix_flag |= Suffix::AUTOPAS;
  respa_enable = 0;
  cut_respa = NULL;
}

/* ---------------------------------------------------------------------- */

void PairLJCutAutoPas::compute(int eflag, int vflag) {
  if (!_isInitialized) {
    init_autopas();
    _isInitialized = true;
  }

  //printf("Simulating computation in AutoPas\n");

  ev_init(eflag, vflag);

  const int nall = atom->nlocal + atom->nghost;
  const int inum = list->inum;

  auto[invalidParticles, updated] = _autopas.updateContainer();

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
  }

  timer->stamp(Timer::PAIR);

  //printf("AutoPas computation\n");


  // Force calculation
  // PairFunctorType functor{_autopas.getCutoff(), *_particlePropertiesLibrary};
  PairFunctorType functor{_autopas.getCutoff()};
  functor.setParticleProperties(24,1);

  _autopas.iteratePairwise(&functor);

  timer->stamp(Timer::AUTOPAS);

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

  std::copy_n(functor.getVirial()->begin(), 6, virial);
  auto upot = functor.getUpot();
  eng_vdwl = upot;

  //printf("AutoPas complete\n");
  _autopas.deleteAllParticles();

}
/* ---------------------------------------------------------------------- */

double PairLJCutAutoPas::memory_usage() {
  double bytes = 0; // memory_usage_thr(); //TODO_AUTOPAS Get memory usage from AutoPas
  bytes += PairLJCut::memory_usage();

  return bytes;
}

void PairLJCutAutoPas::init_autopas() {

  // Initialize particle properties
  _particlePropertiesLibrary = std::make_unique<ParticlePropertiesLibraryType>(
      cut_global);

  for (int i = 1; i <= atom->ntypes; ++i) {
    std::cout << "Type, Eps, Sig: " << i << " " << epsilon[i][i] << " " << sigma[i][i] << "\n";
    _particlePropertiesLibrary->addType(
        i, epsilon[i][i], sigma[i][i], atom->mass[i]
    );
  }

  // _autopas.setAllowedCellSizeFactors(*cellSizeFactors);
  _autopas.setAllowedContainers({autopas::ContainerOption::verletLists, autopas::ContainerOption::linkedCells});
  _autopas.setAllowedDataLayouts({autopas::DataLayoutOption::aos, autopas::DataLayoutOption::soa});
  // _autopas.setAllowedNewton3Options({autopas::Newton3Option::disabled, autopas::Newton3Option::enabled});
  // _autopas.setAllowedTraversals({autopas::TraversalOption::c08,autopas::TraversalOption::c04,autopas::TraversalOption::c04SoA, autopas::TraversalOption::sliced, autopas::TraversalOption::slicedVerlet});

  floatVecType boxMax{}, boxMin{};
  std::copy(std::begin(domain->boxhi), std::end(domain->boxhi), boxMax.begin());
  std::copy(std::begin(domain->boxlo), std::end(domain->boxlo), boxMin.begin());
  _autopas.setBoxMax(boxMax);
  _autopas.setBoxMin(boxMin);

  _autopas.setCutoff(
      cut_global); // TODO_AUTOPAS Test: cut_global (PairLJCut) or cutforce (Pair)
  //_autopas.setNumSamples(tuningSamples);
  //_autopas.setSelectorStrategy(selectorStrategy);
  //_autopas.setTuningInterval(tuningInterval);
  //_autopas.setTuningStrategyOption(tuningStrategy);

  //_autopas.setVerletClusterSize(_config->verletClusterSize);
  _autopas.setVerletRebuildFrequency(neighbor->every);
  _autopas.setVerletSkin(neighbor->skin);

  autopas::Logger::create();
  autopas::Logger::get()->set_level(autopas::Logger::LogLevel::warn);

  _autopas.init();

  neighbor->every = 1; //TODO AutoPas cant handle adding particles that are out of bounds
}
