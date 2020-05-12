#ifndef LMP_AUTOPAS_H
#define LMP_AUTOPAS_H

#include "pointers.h"
#include "atom.h"

#include <cassert>

#include <autopas/AutoPas.h>
#include <autopas/cells/FullParticleCell.h>
#include <autopas/molecularDynamics/ParticlePropertiesLibrary.h>

#include "autopas_lj_functor.h"
#include "autopas_particle.h"

namespace LAMMPS_NS {

class AutoPasLMP : protected Pointers {
public:
  using FloatType = double;
  using FloatVecType = std::array<FloatType, 3>;
  using ParticleType = MoleculeLJLammps<FloatType>;
  using ParticleCellType = autopas::FullParticleCell<ParticleType>;
  using AutoPasType = autopas::AutoPas<ParticleType, ParticleCellType>;
  using ParticlePropertiesLibraryType = ParticlePropertiesLibrary<FloatType, size_t>;
  using PairFunctorType = LJFunctorLammps<ParticleType, ParticleCellType, /*applyShift*/ false, /*useMixing*/ false, /*useNewton3*/ autopas::FunctorN3Modes::Both, /*calculateGlobals*/ true>;

  /*
   * Flag used to differentiate when LAMMPS is build with and without
   * AutoPas support and using AutoPas enabled styles is ok
   */
  int autopas_exists = 1;

  AutoPasLMP(class LAMMPS *, int, char **);

  /**
   * Initialize the AutoPas container
   * @param cutoff LJ Cutoff
   * @param epsilon LJ Epsilon
   * @param sigma LJ Sigma
   */
  void init_autopas(double cutoff, double **epsilon, double **sigma);


  void update_autopas();

  // [[deprecated]]
  [[nodiscard]] ParticleType *particle_by_index(int idx);

  void copy_back();

  template<bool halo = false>
  void add_particle(const ParticleType &p);

  template<bool haloOnly = false>
  autopas::ParticleIteratorWrapper<ParticleType, true>
  particles_by_slab(int i, double d, double d1) const;

  std::vector<ParticleType> &get_leaving_particles();

  template<autopas::IteratorBehavior iterateBehavior>
  AutoPasLMP::AutoPasType::iterator_t iterate() {
    return _autopas->begin(iterateBehavior);
  }

  template<autopas::IteratorBehavior iterateBehavior>
  AutoPasLMP::AutoPasType::const_iterator_t const_iterate() {
    return _autopas->cbegin(iterateBehavior);
  }

  AutoPasLMP::AutoPasType::const_iterator_t
  const_iterate_auto(int first, int last);

  AutoPasLMP::AutoPasType::iterator_t
  iterate_auto(int first, int last);


  bool iterate_pairwise(PairFunctorType *functor);


  bool is_initialized();

  FloatType get_cutoff();

  std::vector<int>
  particle_to_index(const std::vector<ParticleType *> &particles);

  inline int particle_to_index(const ParticleType &particle) {
    //auto idx {atom->map(particle.getID())};
    auto idx{particle.getLocalID()};
    assert(idx != -1);
    return idx;
  }

private:
  std::unique_ptr<AutoPasType> _autopas;
  std::unique_ptr<ParticlePropertiesLibraryType> _particlePropertiesLibrary;
  std::vector<ParticleType> _leavingParticles;

  std::unordered_map<int, ParticleType *> _index_map;
  std::vector<ParticleType *> _index_vector;
  bool _use_index_map = false;
  bool _index_structure_valid = false;

  void update_index_structure();

};
}

#endif
