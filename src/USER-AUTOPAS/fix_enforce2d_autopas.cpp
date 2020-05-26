#include "fix_enforce2d_autopas.h"

#include "atom.h"
#include "autopas.h"

using namespace LAMMPS_NS;

void FixEnforce2DAutoPas::post_force(int vflag) {

  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

#pragma omp parallel default(none) shared(mask, nlocal)
  for (auto iter = lmp->autopas->iterate<autopas::ownedOnly>(); iter.isValid(); ++iter) {
    auto v = iter->getV();
    auto f = iter->getF();
    int idx = AutoPasLMP::particle_to_index(*iter);

    if (idx < nlocal && (mask[idx] & groupbit)) {
      v[2] = 0.0;
      f[2] = 0.0;
      iter->setV(v);
      iter->setF(f);
    }
  }

  // for systems with omega/angmom/torque, zero x and y components

  if (atom->omega_flag) {
    double **omega = atom->omega;
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        omega[i][0] = 0.0;
        omega[i][1] = 0.0;
      }
  }

  if (atom->angmom_flag) {
    double **angmom = atom->angmom;
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        angmom[i][0] = 0.0;
        angmom[i][1] = 0.0;
      }
  }

  if (atom->torque_flag) {
    double **torque = atom->torque;
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        torque[i][0] = 0.0;
        torque[i][1] = 0.0;
      }
  }

  // invoke other fixes that enforce 2d
  // fix rigid variants

  for (int m = 0; m < nfixlist; m++)
    flist[m]->enforce2d();
}
