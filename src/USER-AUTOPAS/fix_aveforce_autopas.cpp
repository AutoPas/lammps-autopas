#include "fix_aveforce_autopas.h"

#include <mpi.h>

#include "atom.h"
#include "autopas.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "region.h"
#include "respa.h"
#include "input.h"
#include "variable.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum {
  NONE, CONSTANT, EQUAL
};

void FixAveForceAutoPas::post_force(int /*vflag*/) {

  // update region if necessary

  Region *region = nullptr;
  if (iregion >= 0) {
    region = domain->regions[iregion];
    region->prematch();
  }

  // sum forces on participating atoms
  double foriginal[4];
  foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;

#pragma omp parallel default(none) shared(region) reduction(+:foriginal[:4])
  for (auto iter = lmp->autopas->iterate<autopas::ownedOnly>(); iter.isValid(); ++iter) {
    auto &x = iter->getR();
    auto &f = iter->getF();
    int idx = AutoPasLMP::particle_to_index(*iter);
    if (atom->mask[idx] & groupbit) {
      if (region && !region->match(x[0], x[1], x[2])) continue;
      foriginal[0] += f[0];
      foriginal[1] += f[1];
      foriginal[2] += f[2];
      foriginal[3] += 1.0;
    }
  }

  // average the force on participating atoms
  // add in requested amount, computed via variable evaluation if necessary
  // wrap variable evaluation with clear/add

  MPI_Allreduce(foriginal, foriginal_all, 4, MPI_DOUBLE, MPI_SUM, world);

  int ncount = static_cast<int> (foriginal_all[3]);
  if (ncount == 0) return;

  if (varflag == EQUAL) {
    modify->clearstep_compute();
    if (xstyle == EQUAL) xvalue = input->variable->compute_equal(xvar);
    if (ystyle == EQUAL) yvalue = input->variable->compute_equal(yvar);
    if (zstyle == EQUAL) zvalue = input->variable->compute_equal(zvar);
    modify->addstep_compute(update->ntimestep + 1);
  }

  double fave[3];
  fave[0] = foriginal_all[0] / ncount + xvalue;
  fave[1] = foriginal_all[1] / ncount + yvalue;
  fave[2] = foriginal_all[2] / ncount + zvalue;

  // set force of all participating atoms to same value
  // only for active dimensions

#pragma omp parallel default(none) shared(region, fave)
  for (auto iter = lmp->autopas->iterate<autopas::ownedOnly>(); iter.isValid(); ++iter) {
    auto &x = iter->getR();
    auto &f = iter->getF();
    int idx = AutoPasLMP::particle_to_index(*iter);
    if (atom->mask[idx] & groupbit) {
      if (region && !region->match(x[0], x[1], x[2])) continue;
      AutoPasLMP::FloatVecType fave_set{f};
      if (xstyle) fave_set[0] = fave[0];
      if (ystyle) fave_set[1] = fave[1];
      if (zstyle) fave_set[2] = fave[2];
      iter->setF(fave_set);
    }
  }
}

void
FixAveForceAutoPas::post_force_respa(int vflag, int ilevel, int /*iloop*/) {

  // ave + extra force on selected RESPA level
  // just ave on all other levels

  if (ilevel == ilevel_respa) post_force(vflag);
  else {
    Region *region = nullptr;
    if (iregion >= 0) {
      region = domain->regions[iregion];
      region->prematch();
    }

    double foriginal[4];
    foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;

#pragma omp parallel default(none) shared(region) reduction(+:foriginal[:4])
    for (auto iter = lmp->autopas->iterate<autopas::ownedOnly>(); iter.isValid(); ++iter) {
      auto &x = iter->getR();
      auto &f = iter->getF();
      int idx = AutoPasLMP::particle_to_index(*iter);
      if (atom->mask[idx] & groupbit) {
        if (region && !region->match(x[0], x[1], x[2])) continue;
        foriginal[0] += f[0];
        foriginal[1] += f[1];
        foriginal[2] += f[2];
        foriginal[3] += 1.0;
      }
    }

    MPI_Allreduce(foriginal, foriginal_all, 4, MPI_DOUBLE, MPI_SUM, world);

    int ncount = static_cast<int> (foriginal_all[3]);
    if (ncount == 0) return;

    double fave[3];
    fave[0] = foriginal_all[0] / ncount;
    fave[1] = foriginal_all[1] / ncount;
    fave[2] = foriginal_all[2] / ncount;

#pragma omp parallel default(none) shared(region, fave)
    for (auto iter = lmp->autopas->iterate<autopas::ownedOnly>(); iter.isValid(); ++iter) {
      auto &x = iter->getR();
      auto &f = iter->getF();
      int idx = AutoPasLMP::particle_to_index(*iter);
      if (atom->mask[idx] & groupbit) {
        if (region && !region->match(x[0], x[1], x[2])) continue;
        AutoPasLMP::FloatVecType fave_set{f};
        if (xstyle) fave_set[0] = fave[0];
        if (ystyle) fave_set[1] = fave[1];
        if (zstyle) fave_set[2] = fave[2];
        iter->setF(fave_set);
      }
    }
  }
}
