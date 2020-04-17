#include "nbin_autopas.h"
#include "autopas.h"
#include "atom.h"
#include "group.h"
#include "update.h"

LAMMPS_NS::NBinAutoPas::NBinAutoPas(LAMMPS_NS::LAMMPS *lmp) : NBinStandard(
    lmp) {

}

void LAMMPS_NS::NBinAutoPas::bin_atoms() {

  int i, ibin;

  last_bin = update->ntimestep;
  for (i = 0; i < mbins; i++) binhead[i] = -1;

  // bin in reverse order so linked list will be in forward order
  // also puts ghost atoms at end of list, which is necessary

  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  if (includegroup) {
    int bitmask = group->bitmask[includegroup];
    for (i = nall - 1; i >= nlocal; i--) {
      if (mask[i] & bitmask) {
        // This are only halo
        auto x = lmp->autopas->particle_by_index(i)->getR();
        ibin = coord2bin(x.data());
        atom2bin[i] = ibin;
        bins[i] = binhead[ibin];
        binhead[ibin] = i;
      }
    }
    for (i = atom->nfirst - 1; i >= 0; i--) {
      // This are only owned
      auto x = lmp->autopas->particle_by_index(i)->getR();
      ibin = coord2bin(x.data());
      atom2bin[i] = ibin;
      bins[i] = binhead[ibin];
      binhead[ibin] = i;
    }

  } else {
    for (i = nall - 1; i >= 0; i--) {
      auto x = lmp->autopas->particle_by_index(i)->getR();
      ibin = coord2bin(x.data());
      atom2bin[i] = ibin;
      bins[i] = binhead[ibin];
      binhead[ibin] = i;
    }
  }
}
