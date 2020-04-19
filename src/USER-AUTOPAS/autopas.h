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

  void init_autopas(double cutoff, double** epsilon, double** sigma);

  template<bool halo=false>
  [[nodiscard]] ParticleType* particle_by_index(int idx);

  unsigned long idx(const ParticleType &p);

};

}