#pragma once

#include <autopas/molecularDynamics/MoleculeLJ.h>
#include <autopas/cells/FullParticleCell.h>
#include <autopas/AutoPas.h>
#include <autopas/molecularDynamics/ParticlePropertiesLibrary.h>

#include "pointers.h"
#include "autopas_lj_functor.h"

// #include "autopas_type.h" // TODO_AP: Same types as LAMMPS
// #include "pair_autopas.h"  // TODO_AP: No default pair style for AP

namespace LAMMPS_NS {

class AutoPasLMP : protected Pointers {
public:
  int autopas_exists;

  AutoPasLMP(class LAMMPS *, int, char **);

  ~AutoPasLMP() override;

  //void accelerator(int, char **);

  //int neigh_count(int);

  using FloatType = double;
  using FloatVecType = std::array<FloatType, 3>;
  using ParticleType = autopas::MoleculeLJ<>;
  using ParticleCellType = autopas::FullParticleCell<ParticleType>;
  using AutoPasType = autopas::AutoPas<ParticleType, ParticleCellType>;
  using ParticlePropertiesLibraryType = ParticlePropertiesLibrary<FloatType, size_t>;
  using PairFunctorType = autopas::LJFunctorLammps<ParticleType, ParticleCellType, /*applyShift*/ false, /*useMixing*/ false, /*useNewton3*/ autopas::FunctorN3Modes::Both, /*calculateGlobals*/ true>;

  std::unique_ptr<AutoPasType> _autopas;
  std::unique_ptr<ParticlePropertiesLibraryType> _particlePropertiesLibrary;

  void init_autopas(double cutoff, double** epsilon, double** sigma);

};

}