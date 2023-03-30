#include "compute_temp_autopas.h"

#include <mpi.h>

#include "atom.h"
#include "autopas.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "update.h"

double LAMMPS_NS::ComputeTempAutoPas::compute_scalar() {
  if (!lmp->autopas->is_initialized()) {
    return ComputeTemp::compute_scalar();
  }

  invoked_scalar = update->ntimestep;

  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;

  double t = 0.0;

#pragma omp parallel default(none) shared(mask, rmass, mass, type) reduction(+: t)
  for (auto iter = lmp->autopas->const_iterate<autopas::IteratorBehavior::owned>(); iter.isValid(); ++iter) {
    auto &v{iter->getV()};
    auto idx{iter->getLocalID()};
    if (mask[idx] & groupbit) {
      auto tmp = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);

      if (rmass) {
        tmp *= rmass[idx];
      } else {
        tmp *= mass[type[idx]];
      }

      t += tmp;
    }
  }


  MPI_Allreduce(&t, &scalar, 1, MPI_DOUBLE, MPI_SUM, world);
  if (dynamic) dof_compute();
  if (dof < 0.0 && natoms_temp > 0.0)
    error->all(FLERR, "Temperature compute degrees of freedom < 0");
  scalar *= tfactor;
  return scalar;
}

void LAMMPS_NS::ComputeTempAutoPas::compute_vector() {
  if (!lmp->autopas->is_initialized()) {
    return ComputeTemp::compute_vector();
  }

  invoked_vector = update->ntimestep;

  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;

  double massone, t[6] = {0};

#pragma omp parallel default(none) shared(mask, rmass, mass, type, massone) reduction(+:t[:6])
  for (auto iter = lmp->autopas->const_iterate<autopas::IteratorBehavior::owned>(); iter.isValid(); ++iter) {
    auto &v{iter->getV()};
    auto idx{iter->getLocalID()};

    if (mask[idx] & groupbit) {
      if (rmass) massone = rmass[idx];
      else massone = mass[type[idx]];
      t[0] += massone * v[0] * v[0];
      t[1] += massone * v[1] * v[1];
      t[2] += massone * v[2] * v[2];
      t[3] += massone * v[0] * v[1];
      t[4] += massone * v[0] * v[2];
      t[5] += massone * v[1] * v[2];
    }
  }

  MPI_Allreduce(t, vector, 6, MPI_DOUBLE, MPI_SUM, world);
  for (int i = 0; i < 6; i++) vector[i] *= force->mvv2e;
}
