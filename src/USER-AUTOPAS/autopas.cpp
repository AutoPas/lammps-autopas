#include "autopas.h"

#include <algorithm>
#include <limits>

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neighbor.h"


using namespace LAMMPS_NS;

AutoPasLMP::AutoPasLMP(class LAMMPS *lmp, int narg, char **args) : Pointers(
    lmp) {
  lmp->autopas = this;
}

void AutoPasLMP::init_autopas(double cutoff, double **epsilon, double **sigma) {

  // TODO SUPPORT FOR: pair_modify shift yes

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

  // TODO Reenable mixings
  if (atom->ntypes > 1) {
    error->warning(FLERR, "Mixings are currently disabled with AutoPas");
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

  {
  //  _autopas->setAllowedContainers({autopas::ContainerOption::linkedCells});
  //  _autopas->setAllowedDataLayouts({autopas::DataLayoutOption::soa});
  //  _autopas->setAllowedTraversals({autopas::TraversalOption::c04});
  }

  // TODO AutoPas always calculates FullShell.
  //  Turn LAMMPS newton setting off to disable force exchange?
  //  Or just leave reverse_comm empty? What about other forces?
  force->newton = force->newton_pair = force->newton_bond = false;

  _autopas->setAllowedNewton3Options(
      {autopas::Newton3Option::enabled, autopas::Newton3Option::disabled});

  _autopas->setBoxMax({domain->boxhi[0], domain->boxhi[1], domain->boxhi[2]});
  _autopas->setBoxMin({domain->boxlo[0], domain->boxlo[1], domain->boxlo[2]});

  _autopas->setCutoff(
      cutoff); // TODO Test: cut_global (PairLJCut) or cutforce (Pair)
  //_autopas.setNumSamples(tuningSamples);
  //_autopas.setSelectorStrategy(selectorStrategy);
  _autopas->setTuningInterval(1000);
  _autopas->setTuningStrategyOption(autopas::TuningStrategyOption::fullSearch);

  //_autopas.setVerletClusterSize(_config->verletClusterSize);

  _autopas->setVerletRebuildFrequency(
      std::max(neighbor->every, neighbor->delay));
  _autopas->setVerletSkin(neighbor->skin);

  { // TODO Necessary until verlet skin works correctly
    neighbor->every = 1;
    _autopas->setVerletRebuildFrequency(1);
    _autopas->setVerletSkin(0); // No skin needed when rebuilding every step
  }

  autopas::Logger::create();
  autopas::Logger::get()->set_level(autopas::Logger::LogLevel::trace);

  _autopas->init();

  // Handle particles that got added before AutoPas was initialized
  move_into();

  // TODO When to move_back()? After run command finished? Destructor?
}

AutoPasLMP::ParticleType *AutoPasLMP::particle_by_index(int idx) {

  if (!_index_structure_valid) {
    update_index_structure();
  }

  if (_use_index_map) {
    return _index_map.at(idx);
  } else {
    return _index_vector.at(idx);
  }

}

void AutoPasLMP::update_autopas() {
  auto&&[invalidParticles, updated] = _autopas->updateContainer();
  _leavingParticles = std::move(invalidParticles);
  _index_structure_valid = false;
}

void AutoPasLMP::move_into() {

#pragma omp parallel for default(none)
  for (int i = 0; i < atom->nlocal; i++) {
    FloatVecType pos{atom->x[i][0], atom->x[i][1], atom->x[i][2]};
    FloatVecType vel{atom->v[i][0], atom->v[i][1], atom->v[i][2]};
    unsigned long moleculeId = i;
    unsigned long typeId = atom->type[i];

    add_particle(ParticleType(pos, vel, moleculeId, typeId));
  }

  // Destroy memory for debugging purposes so we segfault instead of accidentally using the old particles outside of AutoPas
  memory->destroy(atom->x);
  memory->destroy(atom->v);

  _initialized = true;
}

void AutoPasLMP::move_back() {
  copy_back();
  _autopas->deleteAllParticles();
  _initialized = false;
}

void AutoPasLMP::copy_back() const {
  auto nmax = atom->nlocal;
  atom->x = memory->grow(atom->x, nmax, 3, "atom:x");
  atom->v = memory->grow(atom->v, nmax, 3, "atom:v");
  atom->f = memory->grow(atom->f, nmax, 3, "atom:f");

#pragma omp parallel default(none)
  for (auto iter = const_iterate<autopas::ownedOnly>(); iter.isValid(); ++iter) {
    auto &p = *iter;
    auto idx = particle_to_index(p);
    std::copy_n(p.getR().begin(), 3, atom->x[idx]);
    std::copy_n(p.getV().begin(), 3, atom->v[idx]);
    std::copy_n(p.getF().begin(), 3, atom->f[idx]);
  }

  for (auto &p: _leavingParticles) {
    auto idx = particle_to_index(p);
    std::copy_n(p.getR().begin(), 3, atom->x[idx]);
    std::copy_n(p.getV().begin(), 3, atom->v[idx]);
    std::copy_n(p.getF().begin(), 3, atom->f[idx]);
  }
}

void AutoPasLMP::update_index_structure() {
  auto n = _autopas->getNumberOfParticles(autopas::haloAndOwned);
  if (_use_index_map) {
    _index_map.clear();
    _index_map.reserve(n);
  } else {
    _index_vector.clear();
    _index_vector.resize(n);
  }

#pragma omp parallel default(none)
  for (auto &p: *lmp->autopas->_autopas) { // owned and halo
    auto idx = particle_to_index(p);
    if (_use_index_map) {
      _index_map.emplace(idx, &p);
      // TODO Test if no race when emplacing in map
    } else {
      _index_vector[idx] = &p;
    }
  }
  _index_structure_valid = true;
}

#pragma clang diagnostic push
#pragma ide diagnostic ignored "UnreachableCode"

template<bool halo>
void AutoPasLMP::add_particle(const ParticleType &p) {

  if constexpr (halo) {
    _autopas->addOrUpdateHaloParticle(p);
  } else {
    _autopas->addParticle(p);
  }

  if (_index_structure_valid) {
    _index_structure_valid = false;
    //std::cout << "Invalidate map\n";
  }
}

#pragma clang diagnostic pop

template<bool haloOnly>
autopas::ParticleIteratorWrapper<AutoPasLMP::ParticleType, true>
AutoPasLMP::particles_by_slab(int dim, double lo, double hi) const {

  std::array<double, 3> low{};
  std::array<double, 3> high{};

  low.fill(-std::numeric_limits<double>::max());
  high.fill(std::numeric_limits<double>::max());

  low[dim] = lo;
  high[dim] = hi;

  if constexpr (haloOnly) {
    return _autopas->getRegionIterator(low, high, autopas::haloOnly);
  } else {
    return _autopas->getRegionIterator(low, high, autopas::haloAndOwned);
  }
}

std::vector<AutoPasLMP::ParticleType> &
AutoPasLMP::get_leaving_particles() {
  return _leavingParticles;
}

bool AutoPasLMP::is_initialized() {
  return _initialized;
}


std::vector<int>
AutoPasLMP::particle_to_index(
    const std::vector<AutoPasLMP::ParticleType *> &particles) {
  std::vector<int> list(particles.size());
  std::transform(particles.begin(),
                 particles.end(), list.begin(),
                 [](auto *p) { return particle_to_index(*p); });
  return list;
}

AutoPasLMP::FloatType AutoPasLMP::get_cutoff() {
  return _autopas->getCutoff();
}

bool AutoPasLMP::iterate_pairwise(AutoPasLMP::PairFunctorType *functor) {
  return _autopas->iteratePairwise(functor);
}

AutoPasLMP::AutoPasType::const_iterator_t
AutoPasLMP::const_iterate_auto(int first, int last) {
  auto nlocal{atom->nlocal};
  if (last < nlocal) {
    return const_iterate<autopas::ownedOnly>();
  } else if (first >= nlocal) {
    return const_iterate<autopas::haloOnly>();
  } else {
    return const_iterate<autopas::haloAndOwned>();
  }

}

AutoPasLMP::AutoPasType::iterator_t
AutoPasLMP::iterate_auto(int first, int last) {
  auto nlocal{atom->nlocal};
  if (last < nlocal) {
    return iterate<autopas::ownedOnly>();
  } else if (first >= nlocal) {
    return iterate<autopas::haloOnly>();
  } else {
    return iterate<autopas::haloAndOwned>();
  }
}

void AutoPasLMP::update_domain_size() {
  // For changing the domain it is necessary to reinitialize AutoPas
  // All particles will be deleted
  if (_initialized) {
    error->one(FLERR,
               "Loosing all particles, use move_back() before modifying particles outside of AutoPas");
  }

  _autopas->setBoxMax({domain->boxhi[0], domain->boxhi[1], domain->boxhi[2]});
  _autopas->setBoxMin({domain->boxlo[0], domain->boxlo[1], domain->boxlo[2]});

  _autopas->init();
}

template void AutoPasLMP::add_particle<false>(const ParticleType &p);

template void AutoPasLMP::add_particle<true>(const ParticleType &p);

template
autopas::ParticleIteratorWrapper<AutoPasLMP::ParticleType, true>
AutoPasLMP::particles_by_slab<true>(int dim, double lo, double hi) const;

template
autopas::ParticleIteratorWrapper<AutoPasLMP::ParticleType, true>
AutoPasLMP::particles_by_slab<false>(int dim, double lo, double hi) const;
