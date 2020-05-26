#include "fix_indent_autopas.h"

#include "atom.h"
#include "autopas.h"
#include "input.h"
#include "variable.h"
#include "domain.h"
#include "lattice.h"
#include "update.h"
#include "modify.h"
#include "respa.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum {
  NONE, SPHERE, CYLINDER, PLANE
};
enum {
  INSIDE, OUTSIDE
};

void FixIndentAutoPas::post_force(int vflag) {

  // indenter values, 0 = energy, 1-3 = force components
  // wrap variable evaluations with clear/add

  if (varflag) modify->clearstep_compute();

  indenter_flag = 0;
  indenter[0] = indenter[1] = indenter[2] = indenter[3] = 0.0;

  // spherical indenter

  if (istyle == SPHERE) {

    // ctr = current indenter center
    // remap into periodic box

    double ctr[3];
    if (xstr) ctr[0] = input->variable->compute_equal(xvar);
    else ctr[0] = xvalue;
    if (ystr) ctr[1] = input->variable->compute_equal(yvar);
    else ctr[1] = yvalue;
    if (zstr) ctr[2] = input->variable->compute_equal(zvar);
    else ctr[2] = zvalue;
    domain->remap(ctr);

    double radius;
    if (rstr) radius = input->variable->compute_equal(rvar);
    else radius = rvalue;


#pragma omp parallel default(none) shared(ctr, radius)
    for (auto iter = lmp->autopas->iterate<autopas::ownedOnly>(); iter.isValid(); ++iter) {
      auto &x = iter->getR();
      auto &f = iter->getF();
      int idx = AutoPasLMP::particle_to_index(*iter);
      double delx, dely, delz, r, dr, fmag, fx, fy, fz;
      if (atom->mask[idx] & groupbit) {
        delx = x[0] - ctr[0];
        dely = x[1] - ctr[1];
        delz = x[2] - ctr[2];
        domain->minimum_image(delx, dely, delz);
        r = sqrt(delx * delx + dely * dely + delz * delz);
        if (side == OUTSIDE) {
          dr = r - radius;
          fmag = k * dr * dr;
        } else {
          dr = radius - r;
          fmag = -k * dr * dr;
        }
        if (dr >= 0.0) continue;
        fx = delx * fmag / r;
        fy = dely * fmag / r;
        fz = delz * fmag / r;
        iter->addF({fx, fy, fz});
        indenter[0] -= k3 * dr * dr * dr;
        indenter[1] -= fx;
        indenter[2] -= fy;
        indenter[3] -= fz;
      }
    }

    // cylindrical indenter

  } else if (istyle == CYLINDER) {

    // ctr = current indenter axis
    // remap into periodic box
    // 3rd coord is just near box for remap(), since isn't used

    double ctr[3];
    if (cdim == 0) {
      ctr[0] = domain->boxlo[0];
      if (ystr) ctr[1] = input->variable->compute_equal(yvar);
      else ctr[1] = yvalue;
      if (zstr) ctr[2] = input->variable->compute_equal(zvar);
      else ctr[2] = zvalue;
    } else if (cdim == 1) {
      if (xstr) ctr[0] = input->variable->compute_equal(xvar);
      else ctr[0] = xvalue;
      ctr[1] = domain->boxlo[1];
      if (zstr) ctr[2] = input->variable->compute_equal(zvar);
      else ctr[2] = zvalue;
    } else {
      if (xstr) ctr[0] = input->variable->compute_equal(xvar);
      else ctr[0] = xvalue;
      if (ystr) ctr[1] = input->variable->compute_equal(yvar);
      else ctr[1] = yvalue;
      ctr[2] = domain->boxlo[2];
    }
    domain->remap(ctr);

    double radius;
    if (rstr) radius = input->variable->compute_equal(rvar);
    else radius = rvalue;


#pragma omp parallel default(none) shared(ctr, radius)
    for (auto iter = lmp->autopas->iterate<autopas::ownedOnly>(); iter.isValid(); ++iter) {
      auto &x = iter->getR();
      auto &f = iter->getF();
      int idx = AutoPasLMP::particle_to_index(*iter);
      double delx, dely, delz, r, dr, fmag, fx, fy, fz;
      if (atom->mask[idx] & groupbit) {
        if (cdim == 0) {
          delx = 0;
          dely = x[1] - ctr[1];
          delz = x[2] - ctr[2];
        } else if (cdim == 1) {
          delx = x[0] - ctr[0];
          dely = 0;
          delz = x[2] - ctr[2];
        } else {
          delx = x[0] - ctr[0];
          dely = x[1] - ctr[1];
          delz = 0;
        }
        domain->minimum_image(delx, dely, delz);
        r = sqrt(delx * delx + dely * dely + delz * delz);
        if (side == OUTSIDE) {
          dr = r - radius;
          fmag = k * dr * dr;
        } else {
          dr = radius - r;
          fmag = -k * dr * dr;
        }
        if (dr >= 0.0) continue;
        fx = delx * fmag / r;
        fy = dely * fmag / r;
        fz = delz * fmag / r;
        iter->addF({fx, fy, fz});
        indenter[0] -= k3 * dr * dr * dr;
        indenter[1] -= fx;
        indenter[2] -= fy;
        indenter[3] -= fz;
      }
    }

    // planar indenter

  } else {

    // plane = current plane position

    double plane;
    if (pstr) plane = input->variable->compute_equal(pvar);
    else plane = pvalue;

#pragma omp parallel default(none) shared(plane)
    for (auto iter = lmp->autopas->iterate<autopas::ownedOnly>(); iter.isValid(); ++iter) {
      auto &x = iter->getR();
      auto &f = iter->getF();
      int idx = AutoPasLMP::particle_to_index(*iter);
      double dr, fatom;
      if (atom->mask[idx] & groupbit) {
        dr = planeside * (plane - x[cdim]);
        if (dr >= 0.0) continue;
        fatom = -planeside * k * dr * dr;

        std::array<double, 3> f_add{0};
        f_add[cdim] += fatom;
        iter->addF(f_add);

        indenter[0] -= k3 * dr * dr * dr;
        indenter[cdim + 1] -= fatom;
      }
    }
  }

  if (varflag) modify->addstep_compute(update->ntimestep + 1);
}
