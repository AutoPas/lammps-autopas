#ifndef LMP_AUTOPAS_H
#define LMP_AUTOPAS_H

#include "pointers.h"
#include "atom.h"

#include <cassert>

#include <autopas/AutoPas.h>
#include <autopas/cells/FullParticleCell.h>
#include <autopas/molecularDynamics/ParticlePropertiesLibrary.h>

#include "autopas_lj_functor.h"
#include "autopas_lj_functor_avx.h"
#include "autopas_particle.h"

namespace LAMMPS_NS {

class AutoPasLMP : protected Pointers {
public:

  // AutoPas Types
  using FloatType = double;
  using FloatVecType = std::array<FloatType, 3>;
  using ParticleType = MoleculeLJLammps<FloatType>;
  using ParticleCellType = autopas::FullParticleCell<ParticleType>;
  using AutoPasType = autopas::AutoPas<ParticleType, ParticleCellType>;
  using ParticlePropertiesLibraryType = ParticlePropertiesLibrary<FloatType, size_t>;
  using PairFunctorType = LJFunctorAVXLammps<ParticleType, ParticleCellType, /*applyShift*/ false, /*useMixing*/ false, /*useNewton3*/ autopas::FunctorN3Modes::Both, /*calculateGlobals*/ true>;

  /*
   * Flag used to differentiate when LAMMPS is build with and without
   * AutoPas support and when using AutoPas enabled styles is ok
   */
  int autopas_exists = 1;

  /**
   * Create the AutoPas management object
   */
  AutoPasLMP(class LAMMPS *, int, char **);

  /**
   * Initialize the AutoPas container
   * @param cutoff LJ Cutoff
   * @param epsilon LJ Epsilon
   * @param sigma LJ Sigma
   */
  void init_autopas(double cutoff, double **epsilon, double **sigma);

  /**
   * Rebuild the AutoPas container and update the leaving particles with the ones returned by AutoPas.
   * @param must_rebuild When false, let AutoPas decide if rebuild is necessary
   * @return Was container rebuild?
   */
  bool update_autopas(bool must_rebuild);

  /**
   * Get a particle by its global index / particle ID.
   * Modifications to that particle are persistent but might have unexpected results if done the wrong way.
   * E.g. never use this method to update particle positions!
   * @param idx Particle ID
   * @return The pointer to the particle inside of AutoPas.
   * It is only valid when used immediately and must not be stored for later reference.
   */
  [[nodiscard]] ParticleType *
  particle_by_index(int idx); // TODO: Remove if possible

  /**
   * Copies all particles from the AutoPas container into the LAMMPS arrays.
   */
  void copy_back() const;

  /**
   * Moves all particles from the AutoPas container into the LAMMPS arrays and sets AutoPas to uninitialized..
   */
  void move_back();

  /**
   * Moves all particles from the LAMMPS arrays into the AutoPas container and sets AutoPas to initialized.
   */
  void move_into();

  /**
   * Adds a particle to the AutoPas container.
   * @tparam halo True: New particle is a ghost atom; False: New particle is a local atom
   * @param p Particle to add
   */
  template<bool halo = false>
  void add_particle(const ParticleType &p);

  /**
   * Iterate over particles in the specified slab of the domain.
   * @tparam haloOnly Iterate only over ghost atoms
   * @param dim Dimension
   * @param lo Low corner
   * @param hi High corner
   * @return
   */
  template<bool haloOnly = false>
  [[nodiscard]] autopas::ParticleIteratorWrapper<ParticleType, true>
  particles_by_slab(int dim, double lo, double hi) const;

  /**
   * Iterate over all particles in the AutoPas container.
   * @tparam iterateBehavior Local, ghost or both
   * @return Particle iterator
   */
  template<autopas::IteratorBehavior iterateBehavior>
  [[nodiscard]] AutoPasLMP::AutoPasType::iterator_t iterate() {
    return _autopas->begin(iterateBehavior);
  }

  /**
   * Const iterate over all particles in the AutoPas container.
   * @tparam iterateBehavior Local, ghost or both
   * @return Const particle iterator
   */
  template<autopas::IteratorBehavior iterateBehavior>
  [[nodiscard]] AutoPasLMP::AutoPasType::const_iterator_t const_iterate() const {
    return _autopas->cbegin(iterateBehavior);
  }

  /**
   * Iterate over all particles in the AutoPas container.
   * Some particles might be dropped from the iterator based on the given local index bounds.
   * The user MUST still performing an additional check if the index of the particle is inside the given bounds.
   * @param first Lower local index
   * @param last  Upper local index
   * @return Particle iterator
   */
  [[nodiscard]] AutoPasLMP::AutoPasType::const_iterator_t
  const_iterate_auto(int first, int last);

  /**
   * Const iterate over all particles in the AutoPas container.
   * Some particles might be dropped from the iterator based on the given local index bounds.
   * The user MUST still performing an additional check if the index of the particle is inside the given bounds.
   * @param first Lower local index
   * @param last  Upper local index
   * @return Const particle iterator
   */
  [[nodiscard]] AutoPasLMP::AutoPasType::iterator_t
  iterate_auto(int first, int last);


  /**
   * Iterate all particle pairs (considering the set cutoff radius)
   * @param functor Function to execute for every pair
   * @return True if AutoPas tuned in this iteration
   */
  bool iterate_pairwise(PairFunctorType *functor);

  /**
   * Get the local indices of given particles.
   * @param particles Vector of particles
   * @return Vector of local indices
   */
  static std::vector<int>
  particle_to_index(const std::vector<ParticleType *> &particles);

  /**
   * Returns the leaving particles of the last container rebuild.
   * Particles can be deleted from this when they were handled.
   * @return Leaving particles.
   */
  std::vector<ParticleType> &get_leaving_particles();

  /**
   * Flag determining if the particles are currently stored in the AutoPas container or in the LAMMPS arrays.
   * This should always be checked before performing any operations using AutoPas.
   * In most cases a fallback to an LAMMPS only implementation has to be provided if AutoPas is currently not initialized.
   * @return True if particles are stored in the AutoPas container
   */
  bool is_initialized();

  /**
   * Gets the cutoff radius used by the AutoPas container.
   * @return Cutoff radius
   */
  FloatType get_cutoff();

  /**
   * Grow or shrink the AutoPas container to the current domain size.
   * This operation can only be performed if AutoPas is not initialized!
   */
  void update_domain_size();

  /**
   * Compute the interaction mapping which particle type should interact with which.
   * The entry for type i can be found at index (i-1)
   * @return Type x Type => bool mapping
   */
  [[nodiscard]] std::vector<std::vector<int>> get_interaction_map();

  /**
   * Gets the factor the domain should grow every time it has to be resized.
   * @return Regrow factor
   */
  double get_box_grow_factor();

  /**
   * Get the memory usage of AutoPas.
   * @return Memory usage in kB.
   */
  static size_t get_memory_usage();

private:
  std::unique_ptr<AutoPasType> _autopas;
  std::unique_ptr<ParticlePropertiesLibraryType> _particlePropertiesLibrary;
  std::vector<ParticleType> _leavingParticles;

  std::unordered_map<int, ParticleType *> _index_map;
  std::vector<ParticleType *> _index_vector;
  bool _use_index_map = false;
  bool _index_structure_valid = false;

  void update_index_structure();
  void print_config(double *const *epsilon, double *const *sigma) const;

  bool _initialized = false;

  /**
   * AutoPas options
   * defaults from autopas/AutoPas.h
   */
  struct {
    bool tuning{true};
    autopas::Logger::LogLevel log_level{autopas::Logger::LogLevel::warn};

    unsigned int verlet_cluster_size{4};
    unsigned int tuning_interval{5000};
    unsigned int num_samples{3};
    unsigned int max_evidence{10};

    // (only relevant for predictive tuning)
    struct {
      double relative_optimum_range{1.2};
      unsigned int max_tuning_phases_without_test{5};
    } predictive_tuning;

    // / only relevant for Bayesian search)
    autopas::AcquisitionFunctionOption acquisition_function = {
        autopas::AcquisitionFunctionOption::upperConfidenceBound};

    autopas::TuningStrategyOption tuning_strategy{
        autopas::TuningStrategyOption::fullSearch};
    autopas::SelectorStrategyOption selector_strategy{
        autopas::SelectorStrategyOption::fastestAbs};

    std::set<autopas::ContainerOption> allowed_containers{
        autopas::ContainerOption::getMostOptions()};
    std::set<autopas::TraversalOption> allowed_traversals{
        autopas::TraversalOption::getMostOptions()};
    std::set<autopas::DataLayoutOption> allowed_data_layouts{
        autopas::DataLayoutOption::getMostOptions()};
    std::set<autopas::Newton3Option> allowed_newton3{
        autopas::Newton3Option::getMostOptions()};

    // (only relevant for LinkedCells, VerletLists and VerletListsCells)
    std::unique_ptr<autopas::NumberSet<double>> allowed_cell_size_factors{
        std::make_unique<autopas::NumberSetFinite<double>>(
            std::set<double>({1.}))
    };

    // (only relevant for BalancedSlicedTraversal and BalancedSlicedTraversalVerlet).
    std::set<autopas::LoadEstimatorOption> allowed_load_estimators{
        autopas::LoadEstimatorOption::getMostOptions()};
  } _opt;

};
}

#endif
