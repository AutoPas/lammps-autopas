#include "atom_vec_atomic_autopas.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "memory.h"
#include "error.h"
#include "utils.h"

LAMMPS_NS::AtomVecAtomicAutopas::AtomVecAtomicAutopas(LAMMPS_NS::LAMMPS *lmp)
    : AtomVecAutopas(lmp) {
  molecular = 0;
  mass_type = 1;

  comm_x_only = comm_f_only = 1;
  size_forward = 3;
  size_reverse = 3;
  size_border = 6;
  size_velocity = 3;
  size_data_atom = 5;
  size_data_vel = 4;
  xcol_data = 3;
}

void LAMMPS_NS::AtomVecAtomicAutopas::grow(int n) {
  if (n == 0) grow_nmax();
  else nmax = n;
  atom->nmax = nmax;
  if (nmax < 0 || nmax > MAXSMALLINT)
    error->one(FLERR,"Per-processor system is too big");

  tag = memory->grow(atom->tag,nmax,"atom:tag");
  type = memory->grow(atom->type,nmax,"atom:type");
  mask = memory->grow(atom->mask,nmax,"atom:mask");
  image = memory->grow(atom->image,nmax,"atom:image");
  // x = memory->grow(atom->x,nmax,3,"atom:x");
  // v = memory->grow(atom->v,nmax,3,"atom:v");
  // f = memory->grow(atom->f,nmax*comm->nthreads,3,"atom:f");

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->grow_arrays(nmax);
}
void LAMMPS_NS::AtomVecAtomicAutopas::grow_reset() {
  // This is only called from special.cpp line 718
  throw "Not implemented";
}


void LAMMPS_NS::AtomVecAtomicAutopas::copy(int i, int i1, int i2) {
  throw "Not implemented";
}

int
LAMMPS_NS::AtomVecAtomicAutopas::pack_border(int i, int *pInt, double *pDouble,
                                             int i1, int *pInt1) {
  throw "Not implemented";
  return 0;
}

int LAMMPS_NS::AtomVecAtomicAutopas::pack_border_vel(int i, int *pInt,
                                                     double *pDouble, int i1,
                                                     int *pInt1) {
  throw "Not implemented";
  return 0;
}

void
LAMMPS_NS::AtomVecAtomicAutopas::unpack_border(int i, int i1, double *pDouble) {
  throw "Not implemented";

}

void LAMMPS_NS::AtomVecAtomicAutopas::unpack_border_vel(int i, int i1,
                                                        double *pDouble) {
  throw "Not implemented";

}

int LAMMPS_NS::AtomVecAtomicAutopas::pack_exchange(int i, double *pDouble) {
  throw "Not implemented";
  return 0;
}

int LAMMPS_NS::AtomVecAtomicAutopas::unpack_exchange(double *pDouble) {
  throw "Not implemented";
  return 0;
}

int LAMMPS_NS::AtomVecAtomicAutopas::size_restart() {
  throw "Not implemented";
  return 0;
}

int LAMMPS_NS::AtomVecAtomicAutopas::pack_restart(int i, double *pDouble) {
  throw "Not implemented";
  return 0;
}

int LAMMPS_NS::AtomVecAtomicAutopas::unpack_restart(double *pDouble) {
  throw "Not implemented";
  return 0;
}

void LAMMPS_NS::AtomVecAtomicAutopas::create_atom(int itype, double *coord) {

  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  tag[nlocal] = 0;
  type[nlocal] = itype;
  mask[nlocal] = 1;
  image[nlocal] = ((imageint) IMGMAX << IMG2BITS) |
                  ((imageint) IMGMAX << IMGBITS) | IMGMAX;


  AutoPasLMP::FloatVecType pos{coord[0], coord[1], coord[2]};
  AutoPasLMP::FloatVecType  vel{0};
  unsigned long moleculeId = nlocal;
  unsigned long typeId = itype;

  lmp->autopas->addParticle(AutoPasLMP::ParticleType(pos, vel, moleculeId, typeId));

  atom->nlocal++;
}

void LAMMPS_NS::AtomVecAtomicAutopas::data_atom(double *pDouble,
                                                LAMMPS_NS::imageint imageint1,
                                                char **pString) {
  throw "Not implemented";

}

void LAMMPS_NS::AtomVecAtomicAutopas::pack_data(double **pDouble) {
  throw "Not implemented";

}

void LAMMPS_NS::AtomVecAtomicAutopas::write_data(FILE *file, int i,
                                                 double **pDouble) {
  throw "Not implemented";

}

LAMMPS_NS::bigint LAMMPS_NS::AtomVecAtomicAutopas::memory_usage() {
  throw "Not implemented";
  return 0;
}
