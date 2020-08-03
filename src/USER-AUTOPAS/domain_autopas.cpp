#include "domain_autopas.h"

#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "error.h"
#include "fix.h"
#include "fix_deform.h"
#include "force.h"
#include "kspace.h"
#include "lattice.h"
#include "memory.h"
#include "modify.h"
#include "molecule.h"
#include "output.h"
#include "style_region.h"
#include "thermo.h"
#include "universe.h"
#include "update.h"

using namespace LAMMPS_NS;

void DomainAutoPas::pbc() {
  if (!lmp->autopas->is_initialized()) {
    return Domain::pbc();
  }

  // Update leaving particles
  auto &leavingParticles = lmp->autopas->get_leaving_particles();

  // verify owned atoms have valid numerical coords
  // may not if computed pairwise force between 2 atoms at same location

  bool flag = true;
#pragma omp parallel default(none) shared(leavingParticles) reduction(&& : flag)
  {
    for (auto iter = lmp->autopas->const_iterate<autopas::IteratorBehavior::ownedOnly>(); iter.isValid(); ++iter) {
      auto &x{iter->getR()};
      flag &= std::all_of(x.begin(), x.end(),
                          [](auto _) { return std::isfinite(_); });
    }

#pragma omp for
    for (auto iter = leavingParticles.cbegin();
         iter < leavingParticles.cend(); ++iter) {
      auto &x{iter->getR()};
      flag &= std::all_of(x.begin(), x.end(),
                          [](auto _) { return std::isfinite(_); });
    }
  }

  if (!flag) error->one(FLERR, "Non-numeric atom coords - simulation unstable");


  // setup for PBC checks

  double *lo, *hi, *period;
  if (triclinic == 0) {
    lo = boxlo;
    hi = boxhi;
    period = prd;
  } else {
    lo = boxlo_lamda;
    hi = boxhi_lamda;
    period = prd_lamda;
  }

  // apply PBC to each leaving atom
#pragma omp parallel for default(none) shared(hi, lo, period, leavingParticles)
  for (auto iter = leavingParticles.begin();
       iter < leavingParticles.end(); ++iter) {
    pbc(*iter, lo, hi, period);
  }
}

bool
DomainAutoPas::pbc(AutoPasLMP::ParticleType &particle, double *lo, double *hi,
                   double *period) {

  imageint idim, otherdims;
  int *mask = atom->mask;
  imageint *image = atom->image;

  auto x{particle.getR()};
  auto v{particle.getV()};

  auto idx{particle.getLocalID()};
  bool was_touched = false;

  if (xperiodic) {
    if (x[0] < lo[0]) {
      was_touched = true;
      x[0] += period[0];
      if (deform_vremap && mask[idx] & deform_groupbit) v[0] += h_rate[0];
      idim = image[idx] & IMGMASK;
      otherdims = image[idx] ^ idim;
      idim--;
      idim &= IMGMASK;
      image[idx] = otherdims | idim;
    }
    if (x[0] >= hi[0]) {
      was_touched = true;
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
      was_touched = true;
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
      was_touched = true;
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
      was_touched = true;
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
      was_touched = true;
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

  if (was_touched) {
    particle.setR(x);
    particle.setV(v);
  }

  return was_touched;
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
  if (!lmp->autopas->is_initialized()) {
    return Domain::closest_image(i, j);
  }

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
  if (!lmp->autopas->is_initialized()) {
    return Domain::closest_image(pos, j);
  }

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
  if (!lmp->autopas->is_initialized()) {
    return Domain::lamda2x(n);
  }

#pragma omp parallel default(none) shared(n)
  for (auto iter = lmp->autopas->iterate<autopas::IteratorBehavior::ownedOnly>(); iter.isValid(); ++iter) {

    auto idx{iter->getLocalID()};
    if (idx < n) {
      auto x{iter->getR()};

      x[0] = h[0] * x[0] + h[5] * x[1] + h[4] * x[2] + boxlo[0];
      x[1] = h[1] * x[1] + h[3] * x[2] + boxlo[1];
      x[2] = h[2] * x[2] + boxlo[2];

      iter->setR(x);
    }

  }
}

void DomainAutoPas::x2lamda(int n) {
  if (!lmp->autopas->is_initialized()) {
    return Domain::x2lamda(n);
  }

  double delta[3];

#pragma omp parallel default(none) shared(n) private(delta)
  for (auto iter = lmp->autopas->iterate<autopas::IteratorBehavior::ownedOnly>(); iter.isValid(); ++iter) {

    auto idx{iter->getLocalID()};
    if (idx < n) {
      auto x{iter->getR()};

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

void DomainAutoPas::reset_box() {
  if (!lmp->autopas->is_initialized()) {
    return Domain::reset_box();
  }

  // Set a different (larger) small value to grow domain less often
  small[0] = lmp->autopas->get_box_grow_factor() * (boxhi[0] - boxlo[0]);
  small[1] = lmp->autopas->get_box_grow_factor() * (boxhi[1] - boxlo[1]);
  small[2] = lmp->autopas->get_box_grow_factor() * (boxhi[2] - boxlo[2]);

  // Force minimum size to prevent domain shrinking
  minxlo = boxlo[0];
  minxhi = boxhi[0];

  minylo = boxlo[1];
  minyhi = boxhi[1];

  minzlo = boxlo[2];
  minzhi = boxhi[2];

  // perform shrink-wrapping
  // compute extent of atoms on this proc
  // for triclinic, this is done in lamda space

  if (nonperiodic == 2) {
    double extent[3][2], all[3][2];

    double max_x_0, max_x_1, max_x_2;
    double min_x_0, min_x_1, min_x_2;

    max_x_0 = max_x_1 = max_x_2 = std::numeric_limits<double>::min();
    min_x_0 = min_x_1 = min_x_2 = std::numeric_limits<double>::max();

#pragma omp parallel default(none) reduction(min:min_x_0,min_x_1,min_x_2) reduction(max:max_x_0,max_x_1,max_x_2)
    for (auto iter = lmp->autopas->iterate<autopas::ownedOnly>(); iter.isValid(); ++iter) {
      auto &x = iter->getR();

      min_x_0 = std::min(min_x_0, x[0]);
      max_x_0 = std::max(max_x_0, x[0]);
      min_x_1 = std::min(min_x_1, x[1]);
      max_x_1 = std::max(max_x_1, x[1]);
      min_x_2 = std::min(min_x_2, x[2]);
      max_x_2 = std::max(max_x_2, x[2]);

    }

    extent[0][0] = min_x_0;
    extent[0][1] = max_x_0;
    extent[1][0] = min_x_1;
    extent[1][1] = max_x_1;
    extent[2][0] = min_x_2;
    extent[2][1] = max_x_2;

    // compute extent across all procs
    // flip sign of MIN to do it in one Allreduce MAX

    extent[0][0] = -extent[0][0];
    extent[1][0] = -extent[1][0];
    extent[2][0] = -extent[2][0];

    MPI_Allreduce(extent, all, 6, MPI_DOUBLE, MPI_MAX, world);

    // for triclinic, convert back to box coords before changing box

    if (triclinic) lamda2x(atom->nlocal);

    // in shrink-wrapped dims, set box by atom extent
    // if minimum set, enforce min box size settings
    // for triclinic, convert lamda extent to box coords, then set box lo/hi
    // decided NOT to do the next comment - don't want to sneakily change tilt
    // for triclinic, adjust tilt factors if 2nd dim is shrink-wrapped,
    //   so that displacement in 1st dim stays the same

    if (triclinic == 0) {
      if (xperiodic == 0) {
        if (boundary[0][0] == 2) boxlo[0] = -all[0][0] - small[0];
        else if (boundary[0][0] == 3)
          boxlo[0] = std::min(-all[0][0] - small[0], minxlo);
        if (boundary[0][1] == 2) boxhi[0] = all[0][1] + small[0];
        else if (boundary[0][1] == 3)
          boxhi[0] = std::max(all[0][1] + small[0], minxhi);
        if (boxlo[0] > boxhi[0]) error->all(FLERR, "Illegal simulation box");
      }
      if (yperiodic == 0) {
        if (boundary[1][0] == 2) boxlo[1] = -all[1][0] - small[1];
        else if (boundary[1][0] == 3)
          boxlo[1] = std::min(-all[1][0] - small[1], minylo);
        if (boundary[1][1] == 2) boxhi[1] = all[1][1] + small[1];
        else if (boundary[1][1] == 3)
          boxhi[1] = std::max(all[1][1] + small[1], minyhi);
        if (boxlo[1] > boxhi[1]) error->all(FLERR, "Illegal simulation box");
      }
      if (zperiodic == 0) {
        if (boundary[2][0] == 2) boxlo[2] = -all[2][0] - small[2];
        else if (boundary[2][0] == 3)
          boxlo[2] = std::min(-all[2][0] - small[2], minzlo);
        if (boundary[2][1] == 2) boxhi[2] = all[2][1] + small[2];
        else if (boundary[2][1] == 3)
          boxhi[2] = std::max(all[2][1] + small[2], minzhi);
        if (boxlo[2] > boxhi[2]) error->all(FLERR, "Illegal simulation box");
      }

    } else {
      double lo[3], hi[3];
      if (xperiodic == 0) {
        lo[0] = -all[0][0];
        lo[1] = 0.0;
        lo[2] = 0.0;
        Domain::lamda2x(lo, lo);
        hi[0] = all[0][1];
        hi[1] = 0.0;
        hi[2] = 0.0;
        Domain::lamda2x(hi, hi);
        if (boundary[0][0] == 2) boxlo[0] = lo[0] - small[0];
        else if (boundary[0][0] == 3)
          boxlo[0] = std::min(lo[0] - small[0], minxlo);
        if (boundary[0][1] == 2) boxhi[0] = hi[0] + small[0];
        else if (boundary[0][1] == 3)
          boxhi[0] = std::max(hi[0] + small[0], minxhi);
        if (boxlo[0] > boxhi[0]) error->all(FLERR, "Illegal simulation box");
      }
      if (yperiodic == 0) {
        lo[0] = 0.0;
        lo[1] = -all[1][0];
        lo[2] = 0.0;
        Domain::lamda2x(lo, lo);
        hi[0] = 0.0;
        hi[1] = all[1][1];
        hi[2] = 0.0;
        Domain::lamda2x(hi, hi);
        if (boundary[1][0] == 2) boxlo[1] = lo[1] - small[1];
        else if (boundary[1][0] == 3)
          boxlo[1] = std::min(lo[1] - small[1], minylo);
        if (boundary[1][1] == 2) boxhi[1] = hi[1] + small[1];
        else if (boundary[1][1] == 3)
          boxhi[1] = std::max(hi[1] + small[1], minyhi);
        if (boxlo[1] > boxhi[1]) error->all(FLERR, "Illegal simulation box");
        //xy *= (boxhi[1]-boxlo[1]) / yprd;
      }
      if (zperiodic == 0) {
        lo[0] = 0.0;
        lo[1] = 0.0;
        lo[2] = -all[2][0];
        Domain::lamda2x(lo, lo);
        hi[0] = 0.0;
        hi[1] = 0.0;
        hi[2] = all[2][1];
        Domain::lamda2x(hi, hi);
        if (boundary[2][0] == 2) boxlo[2] = lo[2] - small[2];
        else if (boundary[2][0] == 3)
          boxlo[2] = std::min(lo[2] - small[2], minzlo);
        if (boundary[2][1] == 2) boxhi[2] = hi[2] + small[2];
        else if (boundary[2][1] == 3)
          boxhi[2] = std::max(hi[2] + small[2], minzhi);
        if (boxlo[2] > boxhi[2]) error->all(FLERR, "Illegal simulation box");
        //xz *= (boxhi[2]-boxlo[2]) / xprd;
        //yz *= (boxhi[2]-boxlo[2]) / yprd;
      }
    }
  }

  // if box size changed
  if (minxlo != boxlo[0] || minxhi != boxhi[0] || minylo != boxlo[1] ||
      minyhi != boxhi[1] || minzlo != boxlo[2] || minzhi != boxhi[2]) {

    // Move particles back from AutoPas to LAMMPS arrays
    lmp->autopas->move_back();

    // reset box whether shrink-wrapping or not

    set_global_box();
    set_local_box();

    // if shrink-wrapped & kspace is defined (i.e. using MSM), call setup()
    // also call init() (to test for compatibility) ?

    if (nonperiodic == 2 && force->kspace) {
      //force->kspace->init();
      force->kspace->setup();
    }

    // if shrink-wrapped & triclinic, re-convert to lamda coords for new box
    // re-invoke pbc() b/c x2lamda result can be outside [0,1] due to roundoff

    if (nonperiodic == 2 && triclinic) {
      x2lamda(atom->nlocal);
      pbc();
    }

    // Update AutoPas container and move particles back into it
    lmp->autopas->update_domain_size();
    lmp->autopas->move_into();
  }

}
