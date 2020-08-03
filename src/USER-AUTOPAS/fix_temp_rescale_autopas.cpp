#include "fix_temp_rescale_autopas.h"

#include "atom.h"
#include "autopas.h"
#include "force.h"
#include "group.h"
#include "update.h"
#include "comm.h"
#include "input.h"
#include "variable.h"
#include "modify.h"
#include "compute.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum {
  NOBIAS, BIAS
};
enum {
  CONSTANT, EQUAL
};

void FixTempRescaleAutoPas::end_of_step() {

  double t_current = temperature->compute_scalar();

  // there is nothing to do, if there are no degrees of freedom

  if (temperature->dof < 1) return;

  // protect against division by zero

  if (t_current == 0.0)
    error->all(FLERR,
               "Computed temperature for fix temp/rescale cannot be 0.0");

  double delta = update->ntimestep - update->beginstep;
  if (delta != 0.0) delta /= update->endstep - update->beginstep;

  // set current t_target
  // if variable temp, evaluate variable, wrap with clear/add

  if (tstyle == CONSTANT)
    t_target = t_start + delta * (t_stop - t_start);
  else {
    modify->clearstep_compute();
    t_target = input->variable->compute_equal(tvar);
    if (t_target < 0.0)
      error->one(FLERR,
                 "Fix temp/rescale variable returned negative temperature");
    modify->addstep_compute(update->ntimestep + nevery);
  }

  // rescale velocity of appropriate atoms if outside window
  // for BIAS:
  //   temperature is current, so do not need to re-compute
  //   OK to not test returned v = 0, since factor is multiplied by v

  if (fabs(t_current - t_target) > t_window) {
    t_target = t_current - fraction * (t_current - t_target);
    double factor = sqrt(t_target / t_current);
    double efactor = 0.5 * force->boltz * temperature->dof;

    energy += (t_current - t_target) * efactor;

#pragma omp parallel default(none) shared(factor)
    for (auto iter = lmp->autopas->iterate<autopas::ownedOnly>(); iter.isValid(); ++iter) {
      auto v = iter->getV();
      auto &f = iter->getF();
      int idx = iter->getLocalID();
      if (atom->mask[idx] & groupbit) {
        if (which == BIAS) temperature->remove_bias(idx, v.data());
        v[0] *= factor;
        v[1] *= factor;
        v[2] *= factor;
        if (which == BIAS) temperature->restore_bias(idx, v.data());
        iter->setV(v);
      }

    }
  }
}
