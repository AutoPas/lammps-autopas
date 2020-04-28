#pragma once

#include <autopas/molecularDynamics/MoleculeLJ.h>
#include <autopas/cells/FullParticleCell.h>
#include <autopas/AutoPas.h>
#include <autopas/molecularDynamics/ParticlePropertiesLibrary.h>

#include "pointers.h"
#include "autopas_lj_functor.h"

namespace LAMMPS_NS {

class AutoPasLMP : protected Pointers {
public:
  using FloatType = double;
  using FloatVecType = std::array<FloatType, 3>;
  using ParticleType = autopas::MoleculeLJ<>;
  using ParticleCellType = autopas::FullParticleCell<ParticleType>;
  using AutoPasType = autopas::AutoPas<ParticleType, ParticleCellType>;
  using ParticlePropertiesLibraryType = ParticlePropertiesLibrary<FloatType, size_t>;
  using PairFunctorType = autopas::LJFunctorLammps<ParticleType, ParticleCellType, /*applyShift*/ false, /*useMixing*/ false, /*useNewton3*/ autopas::FunctorN3Modes::Both, /*calculateGlobals*/ true>;

  int autopas_exists;

  AutoPasLMP(class LAMMPS *, int, char **);

  ~AutoPasLMP() override;

  std::unique_ptr<AutoPasType> _autopas;
  std::unique_ptr<ParticlePropertiesLibraryType> _particlePropertiesLibrary;
  std::vector<ParticleType> _leavingParticles;

  void init_autopas(double cutoff, double** epsilon, double** sigma);
  void update_autopas();

  [[nodiscard]] ParticleType* particle_by_index(int idx);

  unsigned long idx(const ParticleType &p);

  void copy_back();

  template<bool halo=false>
  void add_particle(const ParticleType &p);

private:
  std::unordered_map<int,AutoPasLMP::ParticleType *> index_map;
  std::vector<AutoPasLMP::ParticleType *> index_vector;
  bool use_index_map = false;
  bool index_structure_valid = false;
  void update_index_structure();

};

}