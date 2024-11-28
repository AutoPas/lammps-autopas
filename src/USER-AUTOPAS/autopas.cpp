#include "autopas.h"

#include <algorithm>
#include <limits>
#include <string>

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neighbor.h"

#include <autopas/utils/MemoryProfiler.h>

using namespace LAMMPS_NS;

AutoPasLMP::AutoPasLMP(class LAMMPS *lmp, int narg, char **argc)
    : Pointers(lmp) {
  lmp->autopas = this;

  // Remove additional traversals from the defaults provided by AutoPas
  _opt.allowed_traversals.erase(
      autopas::TraversalOption::vcl_cluster_iteration); //  Wrong results
  _opt.allowed_traversals.erase(
      autopas::TraversalOption::vcl_c06); // Wrong results
  _opt.allowed_traversals.erase(
      autopas::TraversalOption::vcl_c01_balanced); // Wrong results

  // Command line parsing
  std::vector<std::string> args(argc, argc + narg);

  int iarg = 0;
  while (iarg < narg) {
    if (args[iarg] == "log") {
      if (iarg + 2 > narg)
        error->all(FLERR, "Invalid AutoPas command-line arg: log");
      std::string level = args[iarg + 1];
      if (level == "trace")
        _opt.log_level = autopas::Logger::LogLevel::trace;
      else if (level == "debug")
        _opt.log_level = autopas::Logger::LogLevel::debug;
      else if (level == "info")
        _opt.log_level = autopas::Logger::LogLevel::info;
      else if (level == "warn")
        _opt.log_level = autopas::Logger::LogLevel::warn;
      else if (level == "err")
        _opt.log_level = autopas::Logger::LogLevel::err;
      else if (level == "critical")
        _opt.log_level = autopas::Logger::LogLevel::critical;
      else if (level == "off")
        _opt.log_level = autopas::Logger::LogLevel::off;
      else
        error->all(FLERR, "Invalid AutoPas command-line arg: log");
      iarg += 2;
    } else if (args[iarg] == "n" || args[iarg] == "newton") {
      if (iarg + 2 > narg)
        error->all(FLERR, "Invalid AutoPas command-line arg: newton");
      _opt.allowed_newton3 =
          autopas::Newton3Option::parseOptions(args[iarg + 1]);
      iarg += 2;
    } else if (args[iarg] == "t" || args[iarg] == "traversals") {
      if (iarg + 2 > narg)
        error->all(FLERR, "Invalid AutoPas command-line arg: traversal");
      _opt.allowed_traversals =
          autopas::TraversalOption::parseOptions(args[iarg + 1]);
      iarg += 2;
    } else if (args[iarg] == "c" || args[iarg] == "containers") {
      if (iarg + 2 > narg)
        error->all(FLERR, "Invalid AutoPas command-line arg: containers");
      _opt.allowed_containers =
          autopas::ContainerOption::parseOptions(args[iarg + 1]);
      iarg += 2;
    } else if (args[iarg] == "d" || args[iarg] == "data") {
      if (iarg + 2 > narg)
        error->all(FLERR, "Invalid AutoPas command-line arg: data");
      _opt.allowed_data_layouts =
          autopas::DataLayoutOption::parseOptions(args[iarg + 1]);
      iarg += 2;
    } else if (args[iarg] == "i" || args[iarg] == "interval") {
      if (iarg + 2 > narg)
        error->all(FLERR, "Invalid AutoPas command-line arg: interval");
      _opt.tuning_interval = std::stoi(args[iarg + 1]);
      iarg += 2;
    } else if (args[iarg] == "s" || args[iarg] == "strategies") {
      if (iarg + 2 > narg)
        error->all(FLERR, "Invalid AutoPas command-line arg: strategies");
      _opt.tuning_strategies = autopas::TuningStrategyOption::parseOptions<std::vector<autopas::TuningStrategyOption>>(args[iarg + 1]);
      iarg += 2;
    } else if (args[iarg] == "f" || args[iarg] == "rule_file") {
      if (iarg + 2 > narg)
        error->all(FLERR, "Invalid AutoPas command-line arg: rule_file");
      _opt.rule_filename = args[iarg + 1];
      iarg += 2;
    } else if (args[iarg] == "selector") {
      if (iarg + 2 > narg)
        error->all(FLERR, "Invalid AutoPas command-line arg: selector");
      auto opts = autopas::SelectorStrategyOption::parseOptions(args[iarg + 1]);
      if (opts.size() != 1)
        error->all(FLERR, "Invalid AutoPas command-line arg: selector");
      _opt.selector_strategy = *opts.begin();
      iarg += 2;
    } else if (args[iarg] == "vc_size") {
      if (iarg + 2 > narg)
        error->all(FLERR, "Invalid AutoPas command-line arg: vc_size");
      _opt.verlet_cluster_size = std::stoi(args[iarg + 1]);
      iarg += 2;
    } else if (args[iarg] == "samples") {
      if (iarg + 2 > narg)
        error->all(FLERR, "Invalid AutoPas command-line arg: samples");
      _opt.num_samples = std::stoi(args[iarg + 1]);
      iarg += 2;
    } else if (args[iarg] == "evidence") {
      if (iarg + 2 > narg)
        error->all(FLERR, "Invalid AutoPas command-line arg: evidence");
      _opt.max_evidence = std::stoi(args[iarg + 1]);
      iarg += 2;
    } else if (args[iarg] == "pred_ror") {
      if (iarg + 2 > narg)
        error->all(FLERR, "Invalid AutoPas command-line arg: pred_ror");
      _opt.predictive_tuning.relative_optimum_range = std::stod(args[iarg + 1]);
      iarg += 2;
    } else if (args[iarg] == "pred_mtpwt") {
      if (iarg + 2 > narg)
        error->all(FLERR, "Invalid AutoPas command-line arg: pred_mtpwt");
      _opt.predictive_tuning.max_tuning_phases_without_test =
          std::stoi(args[iarg + 1]);
      iarg += 2;
    } else if (args[iarg] == "bayesian_af") {
      if (iarg + 2 > narg)
        error->all(FLERR, "Invalid AutoPas command-line arg: bayesian_af");
      auto opts =
          autopas::AcquisitionFunctionOption::parseOptions(args[iarg + 1]);
      if (opts.size() != 1)
        error->all(FLERR, "Invalid AutoPas command-line arg: bayesian_af");
      _opt.acquisition_function = *opts.begin();
      iarg += 2;
    } else if (args[iarg] == "csf") {
      if (iarg + 2 > narg)
        error->all(FLERR, "Invalid AutoPas command-line arg: csf");
      auto tokens = autopas::utils::StringUtils::tokenize(
          args[iarg + 1], autopas::utils::StringUtils::delimiters);
      std::set<double> csf;
      for (auto &t : tokens)
        csf.insert(std::stod(t));
      _opt.allowed_cell_size_factors =
          std::make_unique<autopas::NumberSetFinite<double>>(csf);
      iarg += 2;
    } else if (args[iarg] == "estimator") {
      if (iarg + 2 > narg)
        error->all(FLERR, "Invalid AutoPas command-line arg: estimator");
      _opt.allowed_load_estimators =
          autopas::LoadEstimatorOption::parseOptions(args[iarg + 1]);
      iarg += 2;
    } else if (args[iarg] == "notune") {
      _opt.tuning = false;
      iarg += 1;
    } else
      error->all(FLERR, "Invalid AutoPas command-line args");
  }
}

void AutoPasLMP::init_autopas(double cutoff, double **epsilon, double **sigma) {

  if (!lmp->atom->tag_enable) {
    error->one(FLERR, "AutoPas requires particle IDs");
  }

  if (!lmp->atom->map_user) {
    error->one(FLERR, "AutoPas requires global id mapping");
  }

  // TODO SUPPORT FOR: pair_modify shift yes

  _autopas = std::make_unique<AutoPasType>();

  // Initialize particle properties
  _particlePropertiesLibrary =
      std::make_unique<ParticlePropertiesLibraryType>(cutoff);

  print_config(epsilon, sigma);

  // TODO Reenable mixings
  if (atom->ntypes > 1) {
    auto msg = "Mixings are currently disabled with AutoPas! Using epsilon " +
               std::to_string(epsilon[1][1]) + " and sigma " +
               std::to_string(sigma[1][1]) + " for all interactions!";
    error->warning(FLERR, msg.c_str());
  }

  _autopas->setBoxMax({domain->subhi[0], domain->subhi[1], domain->subhi[2]});
  _autopas->setBoxMin({domain->sublo[0], domain->sublo[1], domain->sublo[2]});

  // TODO Test: cut_global (PairLJCut) or cutforce (Pair)
  _autopas->setCutoff(cutoff);

  _autopas->setAllowedCellSizeFactors(*_opt.allowed_cell_size_factors);

  _autopas->setVerletRebuildFrequency(
      std::max(neighbor->every, neighbor->delay));
  _autopas->setVerletSkinPerTimestep(neighbor->skin / _autopas->getVerletRebuildFrequency());
  _autopas->setVerletClusterSize(_opt.verlet_cluster_size);

  _autopas->setTuningInterval(_opt.tuning_interval);
  _autopas->setNumSamples(_opt.num_samples);
  _autopas->setMaxEvidence(_opt.max_evidence);

  _autopas->setRelativeOptimumRange(
      _opt.predictive_tuning.relative_optimum_range);
  _autopas->setMaxTuningPhasesWithoutTest(
      _opt.predictive_tuning.max_tuning_phases_without_test);

  _autopas->setAcquisitionFunction(_opt.acquisition_function);
  _autopas->setSelectorStrategy(_opt.selector_strategy);
  _autopas->setAllowedLoadEstimators(_opt.allowed_load_estimators);

  _autopas->setAllowedContainers(_opt.allowed_containers);
  _autopas->setAllowedTraversals(_opt.allowed_traversals);
  _autopas->setAllowedDataLayouts(_opt.allowed_data_layouts);

  _autopas->setRuleFileName(_opt.rule_filename);
  // TODO AutoPas always calculates FullShell.
  //  Turn LAMMPS newton setting off to disable force exchange?
  //  Or just leave reverse_comm empty? What about other forces?
  force->newton = force->newton_pair = force->newton_bond = false;
  _autopas->setAllowedNewton3Options(_opt.allowed_newton3);

  _autopas->setTuningStrategyOption(_opt.tuning_strategies);

  autopas::Logger::create();
  autopas::Logger::get()->set_level(_opt.log_level);

  _autopas->init();

  // Handle particles that got added before AutoPas was initialized
  move_into();

  // TODO When to move_back()? After run command finished? Destructor?
}

void AutoPasLMP::print_config(double *const *epsilon,
                              double *const *sigma) const {
  // Helper function for printing selected options
  auto printOpt = [](auto name, auto set, auto nameMap) {
    std::cout << "  " << name << " - ";
    for (auto it = set.begin(); it != set.end(); ++it) {
      std::cout << nameMap[*it];
      if (it != std::prev(set.end()))
        std::cout << ", "; // fence post
    }
    std::cout << std::endl;
  };

  // Helper function for printing selected values
  auto printVal = [](auto name, auto set) {
    std::cout << "  " << name << " - ";
    for (auto it = set.begin(); it != set.end(); ++it) {
      std::cout << *it;
      if (it != std::prev(set.end()))
        std::cout << ", "; // fence post
    }
    std::cout << std::endl;
  };

  std::cout << "AutoPas setup:\n";
  printOpt("Traversals", _opt.allowed_traversals,
           autopas::TraversalOption::getOptionNames());
  printOpt("Containers", _opt.allowed_containers,
           autopas::ContainerOption::getOptionNames());

  printOpt("Data Layout", _opt.allowed_data_layouts,
           autopas::DataLayoutOption::getOptionNames());
  printOpt("Newton3", _opt.allowed_newton3,
           autopas::Newton3Option::getOptionNames());

  if (_opt.allowed_traversals.count(
          autopas::TraversalOption::vcl_cluster_iteration) ||
      _opt.allowed_traversals.count(autopas::TraversalOption::vcl_c06) ||
      _opt.allowed_traversals.count(
          autopas::TraversalOption::vcl_c01_balanced)) {
    printVal("Verlet Cluster Size", std::set{_opt.verlet_cluster_size});
  }

  if (_opt.tuning) {
    printVal("Tuning Interval", std::set{_opt.tuning_interval});
    printVal("Num Samples", std::set{_opt.num_samples});
    printVal("Max Evidence", std::set{_opt.max_evidence});

    printOpt("Tuning Strategies", _opt.tuning_strategies,
             autopas::TuningStrategyOption::getOptionNames());
    printOpt("Selector Strategy", std::set{_opt.selector_strategy},
             autopas::SelectorStrategyOption::getOptionNames());

    for (const auto &strategy : _opt.tuning_strategies) {
      if (strategy == autopas::TuningStrategyOption::bayesianSearch ||
          strategy == autopas::TuningStrategyOption::bayesianClusterSearch) {
        printOpt("Acquisition Function", std::set{_opt.acquisition_function},
                 autopas::AcquisitionFunctionOption::getOptionNames());
      } else if (strategy ==
                 autopas::TuningStrategyOption::predictiveTuning) {
        printVal("Relative Optimum Range",
                 std::set{_opt.predictive_tuning.relative_optimum_range});
        printVal("Max Tuning Phases Without Test",
                 std::set{_opt.predictive_tuning.max_tuning_phases_without_test});
      }
    }

    if (_opt.allowed_traversals.count(
            autopas::TraversalOption::lc_sliced_balanced) ||
        _opt.allowed_traversals.count(
            autopas::TraversalOption::vlc_sliced_balanced)) {
      printOpt("Load Estimator", _opt.allowed_load_estimators,
               autopas::LoadEstimatorOption::getOptionNames());
    }
  } else {
    // Double check that AutoPas is not tuning when we are asked to not tune:
    if (_opt.allowed_traversals.size() != 1)
      error->all(FLERR, "AutoPas will tune: Multiple traversal options");
    if (_opt.allowed_containers.size() != 1)
      error->all(FLERR, "AutoPas will tune: Multiple container options");
    if (_opt.allowed_data_layouts.size() != 1)
      error->all(FLERR, "AutoPas will tune: Multiple data layout options");
    if (_opt.allowed_newton3.size() != 1)
      error->all(FLERR, "AutoPas will tune: Multiple newton3 options");
    if (_opt.allowed_cell_size_factors->getAll().size() != 1)
      error->all(FLERR, "AutoPas will tune: Multiple cell size factors");
    if (_opt.allowed_load_estimators.size() != 1) {
      if (_opt.allowed_traversals.count(
              autopas::TraversalOption::lc_sliced_balanced) ||
          _opt.allowed_traversals.count(
              autopas::TraversalOption::vlc_sliced_balanced)) {
        error->all(FLERR, "AutoPas will tune: Multiple load estimators");
      }
    }
  }

  // dummy because LAMMPS starts with typeID 1 and autopas with 0
  _particlePropertiesLibrary->addSiteType(0, 1, 1, 1);
  std::cout << "  Particle Properties\n";
  for (int i = 1; i <= atom->ntypes; ++i) {
    std::cout << "    Type, Eps, Sig: " << i << " " << epsilon[i][i] << " "
              << sigma[i][i] << "\n";
    _particlePropertiesLibrary->addSiteType(i, epsilon[i][i], sigma[i][i],
                                        atom->mass[i]);
  }
  _particlePropertiesLibrary->calculateMixingCoefficients();
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

bool AutoPasLMP::update_autopas() {
  auto invalidParticles = _autopas->updateContainer();
  if (!invalidParticles.empty()) {
    _leavingParticles = std::move(invalidParticles);
    _index_structure_valid = false;
    return true;
  }

  return false;
}

void AutoPasLMP::move_into() {

#pragma omp parallel for default(none)
  for (int i = 0; i < atom->nlocal; i++) {
    FloatVecType pos{atom->x[i][0], atom->x[i][1], atom->x[i][2]};
    FloatVecType vel{atom->v[i][0], atom->v[i][1], atom->v[i][2]};
    unsigned long moleculeId = atom->tag[i];
    unsigned long typeId = atom->type[i];

    add_particle(ParticleType(pos, vel, moleculeId, i, typeId));
  }

  // Destroy memory for debugging purposes so we segfault instead of
  // accidentally using the old particles outside of AutoPas
  memory->destroy(atom->x);
  memory->destroy(atom->v);
  memory->destroy(atom->f);

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
  for (auto iter = const_iterate<autopas::IteratorBehavior::owned>(); iter.isValid();
       ++iter) {
    auto &p = *iter;
    auto idx{iter->getLocalID()};
    std::copy_n(p.getR().begin(), 3, atom->x[idx]);
    std::copy_n(p.getV().begin(), 3, atom->v[idx]);
    std::copy_n(p.getF().begin(), 3, atom->f[idx]);
  }

  for (auto &p : _leavingParticles) {
    auto idx{p.getLocalID()};
    std::copy_n(p.getR().begin(), 3, atom->x[idx]);
    std::copy_n(p.getV().begin(), 3, atom->v[idx]);
    std::copy_n(p.getF().begin(), 3, atom->f[idx]);
  }
}

void AutoPasLMP::update_index_structure() {
  auto n = _autopas->getNumberOfParticles(autopas::IteratorBehavior::ownedOrHalo);
  if (_use_index_map) {
    _index_map.clear();
    _index_map.reserve(n);
  } else {
    _index_vector.clear();
    _index_vector.resize(n);
  }

#pragma omp parallel default(none)
  for (auto &p : *lmp->autopas->_autopas) { // owned and halo
    auto idx{p.getLocalID()};
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

template <bool halo> void AutoPasLMP::add_particle(const ParticleType &p) {
  assert(p.getTypeId() > 0);
  assert(p.getTypeId() <= atom->ntypes);

  if constexpr (halo) {
    _autopas->addHaloParticle(p);
  } else {
    _autopas->addParticle(p);
  }

  if (_index_structure_valid) {
    _index_structure_valid = false;
    // std::cout << "Invalidate map\n";
  }
}

#pragma clang diagnostic pop

template <bool haloOnly>
AutoPasLMP::AutoPasType::RegionIteratorT
AutoPasLMP::particles_by_slab(int dim, double lo, double hi) const {

  std::array<double, 3> low{};
  std::array<double, 3> high{};

  low.fill(std::numeric_limits<double>::lowest());
  high.fill(std::numeric_limits<double>::max());

  low[dim] = lo;
  high[dim] = hi;

  if constexpr (haloOnly) {
    return _autopas->getRegionIterator(low, high, autopas::IteratorBehavior::halo);
  } else {
    return _autopas->getRegionIterator(low, high, autopas::IteratorBehavior::ownedOrHalo);
  }
}

std::vector<AutoPasLMP::ParticleType> &AutoPasLMP::get_leaving_particles() {
  return _leavingParticles;
}

bool AutoPasLMP::is_initialized() { return _initialized; }

std::vector<int> AutoPasLMP::particle_to_index(
    const std::vector<AutoPasLMP::ParticleType *> &particles) {
  std::vector<int> list(particles.size());
  std::transform(particles.begin(), particles.end(), list.begin(),
                 [](auto *p) { return p->getLocalID(); });
  return list;
}

AutoPasLMP::FloatType AutoPasLMP::get_cutoff() { return _autopas->getCutoff(); }

bool AutoPasLMP::iterate_pairwise(AutoPasLMP::PairFunctorType *functor) {
  return _autopas->iteratePairwise(functor);
}

AutoPasLMP::AutoPasType::ConstIteratorT
AutoPasLMP::const_iterate_auto(int first, int last) {
  auto nlocal{atom->nlocal};
  if (last < nlocal) {
    return const_iterate<autopas::IteratorBehavior::owned>();
  } else if (first >= nlocal) {
    return const_iterate<autopas::IteratorBehavior::halo>();
  } else {
    return const_iterate<autopas::IteratorBehavior::ownedOrHalo>();
  }
}

AutoPasLMP::AutoPasType::IteratorT AutoPasLMP::iterate_auto(int first,
                                                             int last) {
  auto nlocal{atom->nlocal};
  if (last < nlocal) {
    return iterate<autopas::IteratorBehavior::owned>();
  } else if (first >= nlocal) {
    return iterate<autopas::IteratorBehavior::halo>();
  } else {
    return iterate<autopas::IteratorBehavior::ownedOrHalo>();
  }
}

void AutoPasLMP::update_domain_size() {
  // For changing the domain it is necessary to reinitialize AutoPas
  // All particles will be deleted
  if (_initialized) {
    error->one(FLERR, "Loosing all particles, use move_back() before modifying "
                      "particles outside of AutoPas");
  }

  _autopas->setBoxMax({domain->subhi[0], domain->subhi[1], domain->subhi[2]});
  _autopas->setBoxMin({domain->sublo[0], domain->sublo[1], domain->sublo[2]});

  _autopas->init();
}

std::vector<std::vector<int>> AutoPasLMP::get_interaction_map() {
  if (lmp->neighbor->nex_group > 0 || lmp->neighbor->nex_mol > 0) {
    error->one(
        FLERR,
        "Excluding interactions based on the group or the molecule is "
        "currently not supported with AutoPas. Use exclusion by type instead.");
  }

  // Initialize n x n map with 1 for every type
  std::vector<std::vector<int>> map(lmp->atom->ntypes,
                                    std::vector<int>(lmp->atom->ntypes, 1));

  // Update value in map for every excluded type pairing
  for (int i = 0; i < lmp->neighbor->nex_type; ++i) {
    map[lmp->neighbor->ex1_type[i] - 1][lmp->neighbor->ex2_type[i] - 1] = 0;
    map[lmp->neighbor->ex2_type[i] - 1][lmp->neighbor->ex1_type[i] - 1] = 0;
  }

  return map;
}

double AutoPasLMP::get_box_grow_factor() {
  // Fraction of total domain size, LAMMPS default is 1.0e-4
  return 1.0e-2;
}

size_t AutoPasLMP::get_memory_usage() {
  return autopas::memoryProfiler::currentMemoryUsage();
}

template void AutoPasLMP::add_particle<false>(const ParticleType &p);

template void AutoPasLMP::add_particle<true>(const ParticleType &p);

template AutoPasLMP::AutoPasType::RegionIteratorT
AutoPasLMP::particles_by_slab<true>(int dim, double lo, double hi) const;

template AutoPasLMP::AutoPasType::RegionIteratorT
AutoPasLMP::particles_by_slab<false>(int dim, double lo, double hi) const;
