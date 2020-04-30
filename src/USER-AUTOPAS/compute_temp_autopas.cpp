#include "compute_temp_autopas.h"
#include <mpi.h>
#include "atom.h"
#include "autopas.h"
#include "update.h"
#include "force.h"
#include "domain.h"
#include "group.h"
#include "error.h"

LAMMPS_NS::ComputeTempAutoPas::ComputeTempAutoPas(
    LAMMPS *lmp, int narg, char **arg) : ComputeTemp(lmp, narg, arg) {

}

double LAMMPS_NS::ComputeTempAutoPas::compute_scalar() {
  auto autopas = lmp->autopas;
  invoked_scalar = update->ntimestep;

  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;

  double t = 0.0;

#pragma omp parallel default(none) shared(autopas, mask, rmass, mass, type) reduction(+: t)
  for (auto iter = lmp->autopas->const_iterate<autopas::IteratorBehavior::ownedOnly>(); iter.isValid(); ++iter) {
    auto &v {iter->getV()};
    auto idx {AutoPasLMP::particle_to_index(*iter)};
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
  auto autopas = lmp->autopas;

  invoked_vector = update->ntimestep;

  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;

  double massone, t[6] = {0};

#pragma omp parallel default(none) shared(autopas, mask, rmass, mass, type, massone, t)
  {
    double t_private[6] = {0};
    for (auto iter = lmp->autopas->const_iterate<autopas::IteratorBehavior::ownedOnly>(); iter.isValid(); ++iter) {
      auto &v {iter->getV()};
      auto idx {AutoPasLMP::particle_to_index(*iter)};

      if (mask[idx] & groupbit) {
        if (rmass) massone = rmass[idx];
        else massone = mass[type[idx]];
        t_private[0] += massone * v[0] * v[0];
        t_private[1] += massone * v[1] * v[1];
        t_private[2] += massone * v[2] * v[2];
        t_private[3] += massone * v[0] * v[1];
        t_private[4] += massone * v[0] * v[2];
        t_private[5] += massone * v[1] * v[2];
      }
    }
#pragma omp critical
    {
      for (int i = 0; i < 6; i++) t[i] += t_private[i];
    }
  }

  MPI_Allreduce(t, vector, 6, MPI_DOUBLE, MPI_SUM, world);
  for (int i = 0; i < 6; i++) vector[i] *= force->mvv2e;
}
