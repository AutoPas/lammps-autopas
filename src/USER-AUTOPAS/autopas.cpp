#include "autopas.h"
#include "atom.h"
#include "domain.h"
#include "neighbor.h"

using namespace LAMMPS_NS;

AutoPasLMP::AutoPasLMP(class LAMMPS *lmp, int narg, char **args) : Pointers(lmp) {
  autopas_exists = 1;
  lmp->autopas = this;
}

AutoPasLMP::~AutoPasLMP() = default;

void AutoPasLMP::init_autopas(double cutoff, double** epsilon, double** sigma) {

  _autopas = std::make_unique<AutoPasType>();

  // Initialize particle properties
  _particlePropertiesLibrary = std::make_unique<ParticlePropertiesLibraryType>(cutoff);

  for (int i = 1; i <= atom->ntypes; ++i) {
    std::cout << "Type, Eps, Sig: " << i << " " << epsilon[i][i] << " "
              << sigma[i][i] << "\n";
    _particlePropertiesLibrary->addType(
        i, epsilon[i][i], sigma[i][i], atom->mass[i]
    );
  }

  // _autopas.setAllowedCellSizeFactors(*cellSizeFactors);
  auto sensibleContainerOptions = autopas::ContainerOption::getAllOptions();
  sensibleContainerOptions.erase(autopas::ContainerOption::directSum); // Never good choice
  _autopas->setAllowedContainers(sensibleContainerOptions);

  auto sensibleTraversalOptions = autopas::TraversalOption::getAllOptions();
  sensibleTraversalOptions.erase(autopas::TraversalOption::verletClusters); //  Segfault
  sensibleTraversalOptions.erase(autopas::TraversalOption::verletClustersColoring); // Segfault
  sensibleTraversalOptions.erase(autopas::TraversalOption::verletClustersStatic); // Segfault
  sensibleTraversalOptions.erase(autopas::TraversalOption::verletClusterCells); // Segfault
  _autopas->setAllowedTraversals(sensibleTraversalOptions);

  _autopas->setAllowedContainers({autopas::ContainerOption::linkedCells});
  _autopas->setAllowedDataLayouts({autopas::DataLayoutOption::soa});
  _autopas->setAllowedNewton3Options({autopas::Newton3Option::enabled});
  _autopas->setAllowedTraversals({autopas::TraversalOption::c04});

  floatVecType boxMax{}, boxMin{};
  std::copy(std::begin(domain->boxhi), std::end(domain->boxhi), boxMax.begin());
  std::copy(std::begin(domain->boxlo), std::end(domain->boxlo), boxMin.begin());
  _autopas->setBoxMax(boxMax);
  _autopas->setBoxMin(boxMin);

  _autopas->setCutoff(cutoff); // TODO_AUTOPAS Test: cut_global (PairLJCut) or cutforce (Pair)
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
}
