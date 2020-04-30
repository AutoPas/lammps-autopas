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

namespace LAMMPS_NS {

class PairLJCutAutoPas : public PairLJCut {
public:
  explicit PairLJCutAutoPas(class LAMMPS *);

  void init_style() override;

  void compute(int, int) override;

  double memory_usage() override;
};

}

#endif
#endif
