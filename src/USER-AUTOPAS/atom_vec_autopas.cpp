#include "atom_vec_autopas.h"

#include "atom.h"
#include "autopas.h"
#include "domain.h"
#include "error.h"

using namespace LAMMPS_NS;

AtomVecAutopas::AtomVecAutopas(LAMMPS_NS::LAMMPS *lmp)
    : AtomVec(lmp) {
  if (!lmp->autopas) {
    error->all(FLERR, "Cannot use AutoPas without the '-autopas on' flag");
  }
}

int
LAMMPS_NS::AtomVecAutopas::pack_comm(int n, int *list, double *buf,
                                     int pbc_flag, int *pbc) {
  int i, j, m;
  double dx, dy, dz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      auto &x = lmp->autopas->particle_by_index(j)->getR();
      buf[m++] = x[0];
      buf[m++] = x[1];
      buf[m++] = x[2];
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0] * domain->xprd;
      dy = pbc[1] * domain->yprd;
      dz = pbc[2] * domain->zprd;
    } else {
      dx = pbc[0] * domain->xprd + pbc[5] * domain->xy + pbc[4] * domain->xz;
      dy = pbc[1] * domain->yprd + pbc[3] * domain->yz;
      dz = pbc[2] * domain->zprd;
    }
    for (i = 0; i < n; i++) {
      j = list[i];
      auto &x = lmp->autopas->particle_by_index(j)->getR();
      buf[m++] = x[0] + dx;
      buf[m++] = x[1] + dy;
      buf[m++] = x[2] + dz;
    }
  }
  return m;
}

int LAMMPS_NS::AtomVecAutopas::pack_comm_vel(int n, int *list, double *buf,
                                             int pbc_flag, int *pbc) {
  int i, j, m;
  double dx, dy, dz, dvx, dvy, dvz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      auto *pj = lmp->autopas->particle_by_index(j);
      auto &x = pj->getR();
      auto &v = pj->getV();
      buf[m++] = x[0];
      buf[m++] = x[1];
      buf[m++] = x[2];
      buf[m++] = v[0];
      buf[m++] = v[1];
      buf[m++] = v[2];
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0] * domain->xprd;
      dy = pbc[1] * domain->yprd;
      dz = pbc[2] * domain->zprd;
    } else {
      dx = pbc[0] * domain->xprd + pbc[5] * domain->xy + pbc[4] * domain->xz;
      dy = pbc[1] * domain->yprd + pbc[3] * domain->yz;
      dz = pbc[2] * domain->zprd;
    }
    if (!deform_vremap) {
      for (i = 0; i < n; i++) {
        j = list[i];
        auto *pj = lmp->autopas->particle_by_index(j);
        auto &x = pj->getR();
        auto &v = pj->getV();
        buf[m++] = x[0] + dx;
        buf[m++] = x[1] + dy;
        buf[m++] = x[2] + dz;
        buf[m++] = v[0];
        buf[m++] = v[1];
        buf[m++] = v[2];
      }
    } else {
      dvx = pbc[0] * h_rate[0] + pbc[5] * h_rate[5] + pbc[4] * h_rate[4];
      dvy = pbc[1] * h_rate[1] + pbc[3] * h_rate[3];
      dvz = pbc[2] * h_rate[2];
      for (i = 0; i < n; i++) {
        j = list[i];
        auto *pj = lmp->autopas->particle_by_index(j);
        auto &x = pj->getR();
        auto &v = pj->getV();
        buf[m++] = x[0] + dx;
        buf[m++] = x[1] + dy;
        buf[m++] = x[2] + dz;
        if (mask[i] & deform_groupbit) {
          buf[m++] = v[0] + dvx;
          buf[m++] = v[1] + dvy;
          buf[m++] = v[2] + dvz;
        } else {
          buf[m++] = v[0];
          buf[m++] = v[1];
          buf[m++] = v[2];
        }
      }
    }
  }
  return m;
}

void LAMMPS_NS::AtomVecAutopas::unpack_comm(int n, int first, double *buf) {

  int m = 0;
  int last = first + n;

  for (auto iter = lmp->autopas->iterate_auto(first, last);
       iter.isValid(); ++iter) {
    auto idx{AutoPasLMP::particle_to_index(*iter)};
    if (idx >= first && idx < last) {
      AutoPasLMP::FloatVecType x_;
      x_[0] = buf[m++];
      x_[1] = buf[m++];
      x_[2] = buf[m++];
      iter->setR(x_);
    }
  }
}

void LAMMPS_NS::AtomVecAutopas::unpack_comm_vel(int n, int first, double *buf) {

  int m = 0;
  int last = first + n;

  for (auto iter = lmp->autopas->iterate_auto(first, last);
       iter.isValid(); ++iter) {
    auto idx{AutoPasLMP::particle_to_index(*iter)};
    if (idx >= first && idx < last) {
      AutoPasLMP::FloatVecType x_;
      AutoPasLMP::FloatVecType v_;
      x_[0] = buf[m++];
      x_[1] = buf[m++];
      x_[2] = buf[m++];
      v_[0] = buf[m++];
      v_[1] = buf[m++];
      v_[2] = buf[m++];
      iter->setR(x_);
      iter->setV(v_);
    }
  }
}

int LAMMPS_NS::AtomVecAutopas::pack_reverse(int n, int first, double *buf) {

  int m = 0;
  int last = first + n;

  for (auto iter = lmp->autopas->const_iterate_auto(first, last);
       iter.isValid(); ++iter) {
    auto idx{AutoPasLMP::particle_to_index(*iter)};
    if (idx >= first && idx < last) {
      auto &f_ = iter->getF();
      buf[m++] = f_[0];
      buf[m++] = f_[1];
      buf[m++] = f_[2];
    }
  }
  return m;
}

void LAMMPS_NS::AtomVecAutopas::unpack_reverse(int n, int *list, double *buf) {
  error->all(FLERR,
             "Function unpack_reverse not supported, use unpack_reverse_autopas instead");
}

void AtomVecAutopas::unpack_reverse_autopas(
    const std::vector<AutoPasLMP::ParticleType *> &list, double *buf) {
  int m = 0;
  for (auto p : list) {
    auto f_{p->getF()};
    f_[0] += buf[m++];
    f_[1] += buf[m++];
    f_[2] += buf[m++];
    p->setF(f_);
  }

}
