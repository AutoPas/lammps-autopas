/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(lj/cut/autopas,PairLJCutAutoPas)

#else

#ifndef LMP_PAIR_LJ_CUT_AUTOPAS_H
#define LMP_PAIR_LJ_CUT_AUTOPAS_H

#include "pair_lj_cut.h"

#include "autopas/AutoPas.h"
#include "autopas_lj_functor.h"
#include "autopas/molecularDynamics/MoleculeLJ.h"

namespace LAMMPS_NS {

class PairLJCutAutoPas : public PairLJCut {

public:
  explicit PairLJCutAutoPas(class LAMMPS *);

  void compute(int, int) override;

  double memory_usage() override;

private:
  using floatType = double;
  using floatVecType = std::array<floatType, 3>;
  using ParticleType = autopas::MoleculeLJ<>;
  using ParticleCellType = autopas::FullParticleCell<ParticleType>;
  using AutoPasType = autopas::AutoPas<ParticleType, ParticleCellType>;
  using ParticlePropertiesLibraryType = ParticlePropertiesLibrary<floatType, size_t>;
  using PairFunctorType = autopas::LJFunctorLammps<ParticleType, ParticleCellType, /*applyShift*/ false, /*useMixing*/ true, /*useNewton3*/ autopas::FunctorN3Modes::Both, /*calculateGlobals*/ true>;

  bool _isInitialized = false;
  AutoPasType _autopas;
  std::unique_ptr<ParticlePropertiesLibraryType> _particlePropertiesLibrary;


//  template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
// void eval(int ifrom, int ito, ThrData * const thr);
  void init_autopas();
};

}

#endif
#endif
