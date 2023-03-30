#include "fix_setforce_autopas.h"

#include "atom.h"
#include "autopas.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "region.h"
#include "respa.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum {
  NONE, CONSTANT, EQUAL, ATOM
};

void FixSetForceAutoPas::post_force(int vflag) {

  // update region if necessary

  Region *region = nullptr;
  if (iregion >= 0) {
    region = domain->regions[iregion];
    region->prematch();
  }

  // reallocate sforce array if necessary

  if (varflag == ATOM && atom->nmax > maxatom) {
    maxatom = atom->nmax;
    memory->destroy(sforce);
    memory->create(sforce, maxatom, 3, "setforce:sforce");
  }

  foriginal[0] = foriginal[1] = foriginal[2] = 0.0;
  force_flag = 0;

  if (varflag == CONSTANT) {
#pragma omp parallel default(none) shared(region) reduction(+:foriginal[:3])
    for (auto iter = lmp->autopas->iterate<autopas::IteratorBehavior::owned>(); iter.isValid(); ++iter) {
      auto &x = iter->getR();
      auto &f = iter->getF();
      auto idx{iter->getLocalID()};
      if (atom->mask[idx] & groupbit) {
        if (region && !region->match(x[0], x[1], x[2])) continue;
        foriginal[0] += f[0];
        foriginal[1] += f[1];
        foriginal[2] += f[2];
        AutoPasLMP::FloatVecType f_set{f};
        if (xstyle) f_set[0] = xvalue;
        if (ystyle) f_set[1] = yvalue;
        if (zstyle) f_set[2] = zvalue;
        iter->setF(f_set);
      }
    }

    // variable force, wrap with clear/add

  } else {

    modify->clearstep_compute();

    if (xstyle == EQUAL) xvalue = input->variable->compute_equal(xvar);
    else if (xstyle == ATOM)
      input->variable->compute_atom(xvar, igroup, &sforce[0][0], 3, 0);
    if (ystyle == EQUAL) yvalue = input->variable->compute_equal(yvar);
    else if (ystyle == ATOM)
      input->variable->compute_atom(yvar, igroup, &sforce[0][1], 3, 0);
    if (zstyle == EQUAL) zvalue = input->variable->compute_equal(zvar);
    else if (zstyle == ATOM)
      input->variable->compute_atom(zvar, igroup, &sforce[0][2], 3, 0);

    modify->addstep_compute(update->ntimestep + 1);

#pragma omp parallel default(none) shared(region) reduction(+:foriginal[:3])
    for (auto iter = lmp->autopas->iterate<autopas::IteratorBehavior::owned>(); iter.isValid(); ++iter) {
      auto &x = iter->getR();
      auto &f = iter->getF();
      auto idx{iter->getLocalID()};
      if (atom->mask[idx] & groupbit) {
        if (region && !region->match(x[0], x[1], x[2])) continue;
        foriginal[0] += f[0];
        foriginal[1] += f[1];
        foriginal[2] += f[2];
        AutoPasLMP::FloatVecType f_set{f};
        if (xstyle == ATOM) f_set[0] = sforce[idx][0];
        else if (xstyle) f_set[0] = xvalue;
        if (ystyle == ATOM) f_set[1] = sforce[idx][1];
        else if (ystyle) f_set[1] = yvalue;
        if (zstyle == ATOM) f_set[2] = sforce[idx][2];
        else if (zstyle) f_set[2] = zvalue;
        iter->setF(f_set);
      }
    }
  }
}

void FixSetForceAutoPas::post_force_respa(int vflag, int ilevel, int iloop) {

  // set force to desired value on requested level, 0.0 on other levels

  if (ilevel == ilevel_respa) post_force(vflag);
  else {
    Region *region = nullptr;
    if (iregion >= 0) {
      region = domain->regions[iregion];
      region->prematch();
    }

#pragma omp parallel default(none) shared(region)
    for (auto iter = lmp->autopas->iterate<autopas::IteratorBehavior::owned>(); iter.isValid(); ++iter) {
      auto &x = iter->getR();
      auto &f = iter->getF();
      auto idx{iter->getLocalID()};
      if (atom->mask[idx] & groupbit) {
        if (region && !region->match(x[0], x[1], x[2])) continue;
        AutoPasLMP::FloatVecType f_set{f};
        if (xstyle) f_set[0] = 0.0;
        if (ystyle) f_set[1] = 0.0;
        if (zstyle) f_set[2] = 0.0;
        iter->setF(f_set);
      }
    }
  }
}
