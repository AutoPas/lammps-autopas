#include "autopas.h"
#include "atom.h"
#include "domain.h"
#include "neighbor.h"

using namespace LAMMPS_NS;

AutoPasLMP::AutoPasLMP(class LAMMPS *lmp, int narg, char **args) : Pointers(
    lmp) {
  autopas_exists = 1;
  lmp->autopas = this;
}

AutoPasLMP::~AutoPasLMP() = default;

void AutoPasLMP::init_autopas(double cutoff, double **epsilon, double **sigma) {

  _autopas = std::make_unique<AutoPasType>();

  // Initialize particle properties
  _particlePropertiesLibrary = std::make_unique<ParticlePropertiesLibraryType>(
      cutoff);

  for (int i = 1; i <= atom->ntypes; ++i) {
    std::cout << "Type, Eps, Sig: " << i << " " << epsilon[i][i] << " "
              << sigma[i][i] << "\n";
    _particlePropertiesLibrary->addType(
        i, epsilon[i][i], sigma[i][i], atom->mass[i]
    );
  }

  // _autopas.setAllowedCellSizeFactors(*cellSizeFactors);
  auto sensibleContainerOptions = autopas::ContainerOption::getAllOptions();
  sensibleContainerOptions.erase(
      autopas::ContainerOption::directSum); // Never good choice
  _autopas->setAllowedContainers(sensibleContainerOptions);

  auto sensibleTraversalOptions = autopas::TraversalOption::getAllOptions();
  sensibleTraversalOptions.erase(
      autopas::TraversalOption::verletClusters); //  Segfault
  sensibleTraversalOptions.erase(
      autopas::TraversalOption::verletClustersColoring); // Segfault
  sensibleTraversalOptions.erase(
      autopas::TraversalOption::verletClustersStatic); // Segfault
  sensibleTraversalOptions.erase(
      autopas::TraversalOption::verletClusterCells); // Segfault
  _autopas->setAllowedTraversals(sensibleTraversalOptions);

  _autopas->setAllowedContainers({autopas::ContainerOption::linkedCells});
  _autopas->setAllowedDataLayouts({autopas::DataLayoutOption::soa});
  _autopas->setAllowedNewton3Options({autopas::Newton3Option::enabled});
  _autopas->setAllowedTraversals({autopas::TraversalOption::c04});

  FloatVecType boxMax{}, boxMin{};
  std::copy(std::begin(domain->boxhi), std::end(domain->boxhi), boxMax.begin());
  std::copy(std::begin(domain->boxlo), std::end(domain->boxlo), boxMin.begin());
  _autopas->setBoxMax(boxMax);
  _autopas->setBoxMin(boxMin);

  _autopas->setCutoff(
      cutoff); // TODO_AUTOPAS Test: cut_global (PairLJCut) or cutforce (Pair)
  //_autopas.setNumSamples(tuningSamples);
  //_autopas.setSelectorStrategy(selectorStrategy);
  _autopas->setTuningInterval(1000);
  _autopas->setTuningStrategyOption(autopas::TuningStrategyOption::fullSearch);

  neighbor->every = 1; //TODO AutoPas cant handle adding particles that are out of bounds
  //_autopas.setVerletClusterSize(_config->verletClusterSize);
  _autopas->setVerletRebuildFrequency(neighbor->every);
  std::cout << neighbor->skin << "\n";
  // _autopas.setVerletSkin(neighbor->skin);
  _autopas->setVerletSkin(0); // No skin needed when rebuilding every step

  autopas::Logger::create();
  autopas::Logger::get()->set_level(autopas::Logger::LogLevel::warn);

  _autopas->init();

  // Handle particles that got added before AutoPas was initialized
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {

    FloatVecType pos{atom->x[i][0], atom->x[i][1], atom->x[i][2]};
    FloatVecType vel{atom->v[i][0], atom->v[i][1], atom->v[i][2]};
    unsigned long moleculeId = i;
    unsigned long typeId = atom->type[i];

    _autopas->addParticle(ParticleType(pos, vel, moleculeId, typeId));
  }

  delete atom->x;
  delete atom->v;
  // TODO When to copy back?
}

AutoPasLMP::ParticleType *AutoPasLMP::particle_by_index(int idx) {

  ParticleType *particle = nullptr;

  // Assuming the particle ids are unique, no reduction is necessary
#pragma omp parallel default(none) shared(idx, particle)
  for (auto iter = _autopas->begin(
      autopas::IteratorBehavior::haloAndOwned); iter.isValid(); ++iter) {
    auto &p = *iter;
    if (p.getID() == idx) {
      particle = &p;
    }
  }

  return particle;

}

