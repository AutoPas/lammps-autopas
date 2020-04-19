#include "domain_autopas.h"
#include "style_region.h"
#include "atom.h"
#include "atom_vec.h"
#include "molecule.h"
#include "force.h"
#include "kspace.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "fix_deform.h"
#include "lattice.h"
#include "comm.h"
#include "output.h"
#include "thermo.h"
#include "universe.h"
#include "memory.h"
#include "error.h"
#include "autopas.h"

using namespace LAMMPS_NS;

DomainAutoPas::DomainAutoPas(LAMMPS_NS::LAMMPS *lmp) : Domain(lmp) {

}

void DomainAutoPas::pbc() {
  auto &autopas = *lmp->autopas->_autopas;

  int i;
  imageint idim, otherdims;
  double *lo, *hi, *period;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  imageint *image = atom->image;

  // verify owned atoms have valid numerical coords
  // may not if computed pairwise force between 2 atoms at same location

  bool flag = true;
#pragma omp parallel default(none) shared(autopas) reduction(&& : flag)
  for (auto iter = autopas.begin(
      autopas::IteratorBehavior::ownedOnly); iter.isValid(); ++iter) {
    auto &particle = *iter;

    auto &x = particle.getR();
    flag &= std::isfinite(x[0]);
    flag &= std::isfinite(x[1]);
    flag &= std::isfinite(x[2]);
  }

  if (!flag) error->one(FLERR, "Non-numeric atom coords - simulation unstable");


  // setup for PBC checks

  if (triclinic == 0) {
    lo = boxlo;
    hi = boxhi;
    period = prd;
  } else {
    lo = boxlo_lamda;
    hi = boxhi_lamda;
    period = prd_lamda;
  }

  // apply PBC to each owned atom

#pragma omp parallel default(none) shared(autopas, hi, lo, period, mask, image) private(idim, otherdims)
  for (auto iter = autopas.begin(
      autopas::IteratorBehavior::ownedOnly); iter.isValid(); ++iter) {
    auto &particle = *iter;

    auto x = particle.getR();
    auto v = particle.getV();

    unsigned long idx = particle.getID(); // TODO Global to local index mapping

    if (xperiodic) {
      if (x[0] < lo[0]) {
        x[0] += period[0];
        if (deform_vremap && mask[idx] & deform_groupbit) v[0] += h_rate[0];
        idim = image[idx] & IMGMASK;
        otherdims = image[idx] ^ idim;
        idim--;
        idim &= IMGMASK;
        image[idx] = otherdims | idim;
      }
      if (x[0] >= hi[0]) {
        x[0] -= period[0];
        x[0] = MAX(x[0], lo[0]);
        if (deform_vremap && mask[idx] & deform_groupbit) v[0] -= h_rate[0];
        idim = image[idx] & IMGMASK;
        otherdims = image[idx] ^ idim;
        idim++;
        idim &= IMGMASK;
        image[idx] = otherdims | idim;
      }
    }

    if (yperiodic) {
      if (x[1] < lo[1]) {
        x[1] += period[1];
        if (deform_vremap && mask[idx] & deform_groupbit) {
          v[0] += h_rate[5];
          v[1] += h_rate[1];
        }
        idim = (image[idx] >> IMGBITS) & IMGMASK;
        otherdims = image[idx] ^ (idim << IMGBITS);
        idim--;
        idim &= IMGMASK;
        image[idx] = otherdims | (idim << IMGBITS);
      }
      if (x[1] >= hi[1]) {
        x[1] -= period[1];
        x[1] = MAX(x[1], lo[1]);
        if (deform_vremap && mask[idx] & deform_groupbit) {
          v[0] -= h_rate[5];
          v[1] -= h_rate[1];
        }
        idim = (image[idx] >> IMGBITS) & IMGMASK;
        otherdims = image[idx] ^ (idim << IMGBITS);
        idim++;
        idim &= IMGMASK;
        image[idx] = otherdims | (idim << IMGBITS);
      }
    }

    if (zperiodic) {
      if (x[2] < lo[2]) {
        x[2] += period[2];
        if (deform_vremap && mask[idx] & deform_groupbit) {
          v[0] += h_rate[4];
          v[1] += h_rate[3];
          v[2] += h_rate[2];
        }
        idim = image[idx] >> IMG2BITS;
        otherdims = image[idx] ^ (idim << IMG2BITS);
        idim--;
        idim &= IMGMASK;
        image[idx] = otherdims | (idim << IMG2BITS);
      }
      if (x[2] >= hi[2]) {
        x[2] -= period[2];
        x[2] = MAX(x[2], lo[2]);
        if (deform_vremap && mask[idx] & deform_groupbit) {
          v[0] -= h_rate[4];
          v[1] -= h_rate[3];
          v[2] -= h_rate[2];
        }
        idim = image[idx] >> IMG2BITS;
        otherdims = image[idx] ^ (idim << IMG2BITS);
        idim++;
        idim &= IMGMASK;
        image[idx] = otherdims | (idim << IMG2BITS);
      }
    }

    particle.setR(x);
    particle.setV(v);
  }
}

void DomainAutoPas::image_check() {
  //TODO This method is only used for molecular systems, but needs adaption
  Domain::image_check();
}

void DomainAutoPas::box_too_small_check() {
  //TODO This method is only used for molecular systems, but needs adaption
  Domain::box_too_small_check();
}

int DomainAutoPas::closest_image(int i, int j) {
  //TODO This method is only used for ntopo_* -> untested
  if (j < 0) return j;

  int *sametag = atom->sametag;
  auto &xi = lmp->autopas->particle_by_index(i)->getR();
  auto &xj = lmp->autopas->particle_by_index(j)->getR();


  int closest = j;
  double delx = xi[0] - xj[0];
  double dely = xi[1] - xj[1];
  double delz = xi[2] - xj[2];
  double rsqmin = delx * delx + dely * dely + delz * delz;
  double rsq;

  while (sametag[j] >= 0) {
    j = sametag[j];
    delx = xi[0] - xj[0];
    dely = xi[1] - xj[1];
    delz = xi[2] - xj[2];
    rsq = delx * delx + dely * dely + delz * delz;
    if (rsq < rsqmin) {
      rsqmin = rsq;
      closest = j;
    }
  }

  return closest;
}

int DomainAutoPas::closest_image(const double *const pos, int j) {
  //TODO This method is not used? -> untested

  if (j < 0) return j;

  const int *const sametag = atom->sametag;
  auto &xj = lmp->autopas->particle_by_index(j)->getR();

  int closest = j;
  double delx = pos[0] - xj[0];
  double dely = pos[1] - xj[1];
  double delz = pos[2] - xj[2];
  double rsqmin = delx * delx + dely * dely + delz * delz;
  double rsq;

  while (sametag[j] >= 0) {
    j = sametag[j];
    auto &xj_ = lmp->autopas->particle_by_index(j)->getR();
    delx = pos[0] - xj_[0];
    dely = pos[1] - xj_[1];
    delz = pos[2] - xj_[2];
    rsq = delx * delx + dely * dely + delz * delz;
    if (rsq < rsqmin) {
      rsqmin = rsq;
      closest = j;
    }
  }

  return closest;
}

void DomainAutoPas::lamda2x(int n) {

#pragma omp parallel default(none) shared(n)
  for (auto iter = lmp->autopas->_autopas->begin(
      autopas::IteratorBehavior::ownedOnly); iter.isValid(); ++iter) {

    if (iter->getID() < n) {
      auto x = iter->getR();

      x[0] = h[0] * x[0] + h[5] * x[1] + h[4] * x[2] + boxlo[0];
      x[1] = h[1] * x[1] + h[3] * x[2] + boxlo[1];
      x[2] = h[2] * x[2] + boxlo[2];

      iter->setR(x);
    }

  }
}

void DomainAutoPas::x2lamda(int n) {
  double delta[3];

#pragma omp parallel default(none) shared(n) private(delta)
  for (auto iter = lmp->autopas->_autopas->begin(
      autopas::IteratorBehavior::ownedOnly); iter.isValid(); ++iter) {

    if (iter->getID() < n) {
      auto x = iter->getR();

      delta[0] = x[0] - boxlo[0];
      delta[1] = x[1] - boxlo[1];
      delta[2] = x[2] - boxlo[2];

      x[0] = h_inv[0] * delta[0] + h_inv[5] * delta[1] + h_inv[4] * delta[2];
      x[1] = h_inv[1] * delta[1] + h_inv[3] * delta[2];
      x[2] = h_inv[2] * delta[2];

      iter->setR(x);
    }
  }
}
