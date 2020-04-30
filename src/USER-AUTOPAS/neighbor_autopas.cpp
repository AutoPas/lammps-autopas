#include "neighbor_autopas.h"

#include <cmath>

#include "autopas.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include "modify.h"
#include "nbin.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "npair.h"
#include "nstencil.h"
#include "pair.h"
#include "update.h"


LAMMPS_NS::NeighborAutoPas::NeighborAutoPas(LAMMPS_NS::LAMMPS *lmp) : Neighbor(
    lmp) {

}

int LAMMPS_NS::NeighborAutoPas::check_distance() {

  double delx, dely, delz, rsq;
  double delta, deltasq, delta1, delta2;

  if (boxcheck) {
    if (triclinic == 0) {
      delx = bboxlo[0] - boxlo_hold[0];
      dely = bboxlo[1] - boxlo_hold[1];
      delz = bboxlo[2] - boxlo_hold[2];
      delta1 = sqrt(delx * delx + dely * dely + delz * delz);
      delx = bboxhi[0] - boxhi_hold[0];
      dely = bboxhi[1] - boxhi_hold[1];
      delz = bboxhi[2] - boxhi_hold[2];
      delta2 = sqrt(delx * delx + dely * dely + delz * delz);
      delta = 0.5 * (skin - (delta1 + delta2));
      deltasq = delta * delta;
    } else {
      domain->box_corners();
      delta1 = delta2 = 0.0;
      for (int i = 0; i < 8; i++) {
        delx = corners[i][0] - corners_hold[i][0];
        dely = corners[i][1] - corners_hold[i][1];
        delz = corners[i][2] - corners_hold[i][2];
        delta = sqrt(delx * delx + dely * dely + delz * delz);
        if (delta > delta1) delta1 = delta;
        else if (delta > delta2) delta2 = delta;
      }
      delta = 0.5 * (skin - (delta1 + delta2));
      deltasq = delta * delta;
    }
  } else deltasq = triggersq;

  int nlocal = atom->nlocal;
  if (includegroup) nlocal = atom->nfirst;

  int flag = 0;
#pragma omp parallel default(none) shared(nlocal, deltasq) private(delx, dely, delz, rsq) reduction(max: flag)

  for (auto iter = lmp->autopas->const_iterate<autopas::IteratorBehavior::ownedOnly>(); iter.isValid(); ++iter) {
    auto &x{iter->getR()};
    auto idx{AutoPasLMP::particle_to_index(*iter)};
    if (idx < nlocal) {
      delx = x[0] - xhold[idx][0];
      dely = x[1] - xhold[idx][1];
      delz = x[2] - xhold[idx][2];
      rsq = delx * delx + dely * dely + delz * delz;
      if (rsq > deltasq) flag = 1;
    }
  }

  int flagall;
  MPI_Allreduce(&flag, &flagall, 1, MPI_INT, MPI_MAX, world);
  if (flagall && ago == MAX(every, delay)) ndanger++;
  return flagall;
}

void LAMMPS_NS::NeighborAutoPas::build(int topoflag) {

  return; //TODO Do we even need the neighbor lists in lammps? Don't create them for now

  int m;

  ago = 0;
  ncalls++;
  lastcall = update->ntimestep;

  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  // check that using special bond flags will not overflow neigh lists

  if (nall > NEIGHMASK)
    error->one(FLERR, "Too many local+ghost atoms for neighbor list");

  // store current atom positions and box size if needed

  if (dist_check) {
    if (includegroup) nlocal = atom->nfirst;
    if (atom->nmax > maxhold) {
      maxhold = atom->nmax;
      memory->destroy(xhold);
      memory->create(xhold, maxhold, 3, "neigh:xhold");
    }

#pragma omp parallel default(none) shared(nlocal)
    for (auto iter = lmp->autopas->const_iterate<autopas::IteratorBehavior::ownedOnly>(); iter.isValid(); ++iter) {
      auto &x{iter->getR()};
      auto idx{AutoPasLMP::particle_to_index(*iter)};
      if (idx < nlocal) {
        xhold[idx][0] = x[0];
        xhold[idx][1] = x[1];
        xhold[idx][2] = x[2];
      }
    }
    if (boxcheck) {
      if (triclinic == 0) {
        boxlo_hold[0] = bboxlo[0];
        boxlo_hold[1] = bboxlo[1];
        boxlo_hold[2] = bboxlo[2];
        boxhi_hold[0] = bboxhi[0];
        boxhi_hold[1] = bboxhi[1];
        boxhi_hold[2] = bboxhi[2];
      } else {
        domain->box_corners();
        corners = domain->corners;
        for (int i = 0; i < 8; i++) {
          corners_hold[i][0] = corners[i][0];
          corners_hold[i][1] = corners[i][1];
          corners_hold[i][2] = corners[i][2];
        }
      }
    }
  }

  // bin atoms for all NBin instances
  // not just NBin associated with perpetual lists, also occasional lists
  // b/c cannot wait to bin occasional lists in build_one() call
  // if bin then, atoms may have moved outside of proc domain & bin extent,
  //   leading to errors or even a crash

  if (style != Neighbor::NSQ) {
    for (int i = 0; i < nbin; i++) {
      neigh_bin[i]->bin_atoms_setup(nall);
      neigh_bin[i]->bin_atoms();
    }
  }

  // build pairwise lists for all perpetual NPair/NeighList
  // grow() with nlocal/nall args so that only realloc if have to

  for (int i = 0; i < npair_perpetual; i++) {
    m = plist[i];
    if (!lists[m]->copy) lists[m]->grow(nlocal, nall);
    neigh_pair[m]->build_setup();
    neigh_pair[m]->build(lists[m]);
  }

  // build topology lists for bonds/angles/etc

  if (atom->molecular && topoflag) build_topology();
}
