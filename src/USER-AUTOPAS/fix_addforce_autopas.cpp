#include "fix_addforce_autopas.h"

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
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum {
  NONE, CONSTANT, EQUAL, ATOM
};

void FixAddForceAutoPas::post_force(int vflag) {

  if (update->ntimestep % nevery) return;

  // energy and virial setup

  if (vflag) v_setup(vflag);
  else evflag = 0;

  // update region if necessary

  Region *region = NULL;
  if (iregion >= 0) {
    region = domain->regions[iregion];
    region->prematch();
  }

  // reallocate sforce array if necessary

  if ((varflag == ATOM || estyle == ATOM) && atom->nmax > maxatom) {
    maxatom = atom->nmax;
    memory->destroy(sforce);
    memory->create(sforce, maxatom, 4, "addforce:sforce");
  }

  // foriginal[0] = "potential energy" for added force
  // foriginal[123] = force on atoms before extra force added

  foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;
  force_flag = 0;

  // constant force
  // potential energy = - x dot f in unwrapped coords

  if (varflag == CONSTANT) {
#pragma omp parallel default(none) shared(region)
    for (auto iter = lmp->autopas->iterate<autopas::ownedOnly>(); iter.isValid(); ++iter) {
      auto &x = iter->getR();
      auto &f = iter->getF();
      int idx = AutoPasLMP::particle_to_index(*iter);
      if (atom->mask[idx] & groupbit) {
        if (region && !region->match(x[0], x[1], x[2])) continue;
        double unwrap[3];
        domain->unmap(x.data(), atom->image[idx], unwrap);
        foriginal[0] -=
            xvalue * unwrap[0] + yvalue * unwrap[1] + zvalue * unwrap[2];
        foriginal[1] += f[0];
        foriginal[2] += f[1];
        foriginal[3] += f[2];
        iter->addF({xvalue, yvalue, zvalue});
        if (evflag) {
          double v[6];
          v[0] = xvalue * unwrap[0];
          v[1] = yvalue * unwrap[1];
          v[2] = zvalue * unwrap[2];
          v[3] = xvalue * unwrap[1];
          v[4] = xvalue * unwrap[2];
          v[5] = yvalue * unwrap[2];
          v_tally(idx, v);
        }
      }
    }

    // variable force, wrap with clear/add
    // potential energy = evar if defined, else 0.0
    // wrap with clear/add

  } else {
    double unwrap[3];

    modify->clearstep_compute();

    if (xstyle == EQUAL) xvalue = input->variable->compute_equal(xvar);
    else if (xstyle == ATOM)
      input->variable->compute_atom(xvar, igroup, &sforce[0][0], 4, 0);
    if (ystyle == EQUAL) yvalue = input->variable->compute_equal(yvar);
    else if (ystyle == ATOM)
      input->variable->compute_atom(yvar, igroup, &sforce[0][1], 4, 0);
    if (zstyle == EQUAL) zvalue = input->variable->compute_equal(zvar);
    else if (zstyle == ATOM)
      input->variable->compute_atom(zvar, igroup, &sforce[0][2], 4, 0);
    if (estyle == ATOM)
      input->variable->compute_atom(evar, igroup, &sforce[0][3], 4, 0);

    modify->addstep_compute(update->ntimestep + 1);

#pragma omp parallel default(none) shared(region)
    for (auto iter = lmp->autopas->iterate<autopas::ownedOnly>(); iter.isValid(); ++iter) {
      auto &x = iter->getR();
      auto &f = iter->getF();
      int idx = AutoPasLMP::particle_to_index(*iter);
      if (atom->mask[idx] & groupbit) {
        if (region && !region->match(x[0], x[1], x[2])) continue;
        double unwrap[3];
        domain->unmap(x.data(), atom->image[idx], unwrap);
        if (xstyle == ATOM) xvalue = sforce[idx][0];
        if (ystyle == ATOM) yvalue = sforce[idx][1];
        if (zstyle == ATOM) zvalue = sforce[idx][2];

        if (estyle == ATOM) {
          foriginal[0] += sforce[idx][3];
        } else {
          if (xstyle) foriginal[0] -= xvalue * unwrap[0];
          if (ystyle) foriginal[0] -= yvalue * unwrap[1];
          if (zstyle) foriginal[0] -= zvalue * unwrap[2];
        }
        foriginal[1] += f[0];
        foriginal[2] += f[1];
        foriginal[3] += f[2];

        std::array<double, 3> f_add{0};

        if (xstyle) f_add[0] += xvalue;
        if (ystyle) f_add[1] += yvalue;
        if (zstyle) f_add[2] += zvalue;

        iter->addF(f_add);

        if (evflag) {
          double v[6];
          v[0] = xstyle ? xvalue * unwrap[0] : 0.0;
          v[1] = ystyle ? yvalue * unwrap[1] : 0.0;
          v[2] = zstyle ? zvalue * unwrap[2] : 0.0;
          v[3] = xstyle ? xvalue * unwrap[1] : 0.0;
          v[4] = xstyle ? xvalue * unwrap[2] : 0.0;
          v[5] = ystyle ? yvalue * unwrap[2] : 0.0;
          v_tally(idx, v);
        }
      }
    }
  }
}
