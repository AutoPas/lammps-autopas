#include "atom_vec_atomic_autopas.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "memory.h"
#include "error.h"
#include "utils.h"

using namespace LAMMPS_NS;

AtomVecAtomicAutopas::AtomVecAtomicAutopas(LAMMPS_NS::LAMMPS *lmp)
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

void AtomVecAtomicAutopas::grow(int n) {
  if (n == 0) grow_nmax();
  else nmax = n;
  atom->nmax = nmax;
  if (nmax < 0 || nmax > MAXSMALLINT)
    error->one(FLERR, "Per-processor system is too big");

  tag = memory->grow(atom->tag, nmax, "atom:tag");
  type = memory->grow(atom->type, nmax, "atom:type");
  mask = memory->grow(atom->mask, nmax, "atom:mask");
  image = memory->grow(atom->image, nmax, "atom:image");

  if (!lmp->autopas->_autopas) { // Only used for startup
    x = memory->grow(atom->x, nmax, 3, "atom:x");
    v = memory->grow(atom->v, nmax, 3, "atom:v");
    // f = memory->grow(atom->f,nmax*comm->nthreads,3,"atom:f");
  }

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->grow_arrays(nmax);
}

void AtomVecAtomicAutopas::grow_reset() {
  // This is only called from special.cpp line 718
  throw "Not implemented";
}


void AtomVecAtomicAutopas::copy(int i, int j, int delflag) {
  tag[j] = tag[i];
  type[j] = type[i];
  mask[j] = mask[i];
  image[j] = image[i];

  // TODO Copy is only used for sorting, which we will not support, or for deleting particles -> no need to handle x and v
  /*
  auto *pj = lmp->autopas->particle_by_index(j);
  auto *pi = lmp->autopas->particle_by_index(i);

  pj->setR(pi->getR());
  pj->setV(pi->getV());
   */

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->copy_arrays(i, j, delflag);
}


int
AtomVecAtomicAutopas::pack_border(int n, int *list, double *buf,
                                             int pbc_flag, int *pbc) {
  error->all(FLERR,
             "Function pack_border not supported, use pack_border_autopas instead");
  return 0;
}

int
AtomVecAtomicAutopas::pack_border_autopas(
    const std::vector<AutoPasLMP::ParticleType *> &particles, double *buf,
    int pbc_flag, const int *pbc) {

  double dx = 0, dy = 0, dz = 0;
  if (pbc_flag != 0) {
    if (domain->triclinic == 0) {
      dx = pbc[0] * domain->xprd;
      dy = pbc[1] * domain->yprd;
      dz = pbc[2] * domain->zprd;
    } else {
      dx = pbc[0];
      dy = pbc[1];
      dz = pbc[2];
    }
  }

  int m = 0;

  for (const auto p : particles) {
    auto idx = lmp->autopas->idx(*p);
    auto &_x = p->getR();
    buf[m++] = _x[0] + dx;
    buf[m++] = _x[1] + dy;
    buf[m++] = _x[2] + dz;
    buf[m++] = ubuf(tag[idx]).d;
    buf[m++] = ubuf(type[idx]).d;
    buf[m++] = ubuf(mask[idx]).d;
  }

  if (atom->nextra_border) {
    std::vector<int> list(particles.size());
    for (const auto &p : particles) {
      list.push_back(p->getID());
    }

    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(list.size(),
                                                                list.data(),
                                                                &buf[m]);
  }

  return m;
}


int
AtomVecAtomicAutopas::pack_border_vel(int n, int *list, double *buf,
                                                 int pbc_flag, int *pbc) {
  error->all(FLERR,
             "Function pack_border_vel not supported, use pack_border_vel_autopas instead");
  return 0;
}

int
AtomVecAtomicAutopas::pack_border_vel_autopas(
    const std::vector<AutoPasLMP::ParticleType *> &particles, double *buf,
    int pbc_flag, const int *pbc) {
  double dx = 0, dy = 0, dz = 0, dvx = 0, dvy = 0, dvz = 0;

  if (pbc_flag != 0) {
    if (domain->triclinic == 0) {
      dx = pbc[0] * domain->xprd;
      dy = pbc[1] * domain->yprd;
      dz = pbc[2] * domain->zprd;
    } else {
      dx = pbc[0];
      dy = pbc[1];
      dz = pbc[2];
    }
    if (deform_vremap) {
      dvx = pbc[0] * h_rate[0] + pbc[5] * h_rate[5] + pbc[4] * h_rate[4];
      dvy = pbc[1] * h_rate[1] + pbc[3] * h_rate[3];
      dvz = pbc[2] * h_rate[2];
    }
  }
  int m = 0;

  for (const auto p : particles) {
    auto idx = lmp->autopas->idx(*p);
    auto &_x = p->getR();
    auto &_v = p->getV();
    buf[m++] = _x[0] + dx;
    buf[m++] = _x[1] + dy;
    buf[m++] = _x[2] + dz;
    buf[m++] = ubuf(tag[idx]).d;
    buf[m++] = ubuf(type[idx]).d;
    buf[m++] = ubuf(mask[idx]).d;

    if (deform_vremap && (mask[idx] & deform_groupbit)) {
      buf[m++] = _v[0] + dvx;
      buf[m++] = _v[1] + dvy;
      buf[m++] = _v[2] + dvz;
    } else {
      buf[m++] = _v[0];
      buf[m++] = _v[1];
      buf[m++] = _v[2];
    }
  }

  if (atom->nextra_border) {
    std::vector<int> list(particles.size());
    for (const auto &p : particles) {
      list.push_back(p->getID());
    }

    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(list.size(),
                                                                list.data(),
                                                                &buf[m]);
  }

  return m;
}

void
AtomVecAtomicAutopas::unpack_border(int n, int first, double *buf) {

  int i, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax) grow(0);

    AutoPasLMP::FloatVecType x;

    x[0] = buf[m++];
    x[1] = buf[m++];
    x[2] = buf[m++];
    tag[i] = (tagint) ubuf(buf[m++]).i;
    type[i] = (int) ubuf(buf[m++]).i;
    mask[i] = (int) ubuf(buf[m++]).i;

    // Always halo particles
    AutoPasLMP::ParticleType pi(x, {0, 0, 0}, static_cast<unsigned long>(i),
                                static_cast<unsigned long>(type[i]));
    lmp->autopas->add_particle</*halo*/ true>(pi);

  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
          unpack_border(n, first, &buf[m]);
}

void AtomVecAtomicAutopas::unpack_border_vel(int n, int first,
                                                        double *buf) {

  int i, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax) grow(0);

    AutoPasLMP::FloatVecType x;
    AutoPasLMP::FloatVecType v;

    x[0] = buf[m++];
    x[1] = buf[m++];
    x[2] = buf[m++];
    tag[i] = (tagint) ubuf(buf[m++]).i;
    type[i] = (int) ubuf(buf[m++]).i;
    mask[i] = (int) ubuf(buf[m++]).i;
    v[0] = buf[m++];
    v[1] = buf[m++];
    v[2] = buf[m++];

    // Always halo particles
    AutoPasLMP::ParticleType pi(x, v, static_cast<unsigned long>(i),
                                static_cast<unsigned long>(type[i]));
    lmp->autopas->add_particle</*halo*/ true>(pi);

  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
          unpack_border(n, first, &buf[m]);
}

int AtomVecAtomicAutopas::pack_exchange(int i, double *buf) {

  auto *pi = lmp->autopas->particle_by_index(i);
  auto &x = pi->getR();
  auto &v = pi->getV();

  int m = 1;
  buf[m++] = x[0];
  buf[m++] = x[1];
  buf[m++] = x[2];
  buf[m++] = v[0];
  buf[m++] = v[1];
  buf[m++] = v[2];
  buf[m++] = ubuf(tag[i]).d;
  buf[m++] = ubuf(type[i]).d;
  buf[m++] = ubuf(mask[i]).d;
  buf[m++] = ubuf(image[i]).d;

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->pack_exchange(i, &buf[m]);

  buf[0] = m;
  return m;
}

int AtomVecAtomicAutopas::unpack_exchange(double *buf) {

  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);


  AutoPasLMP::FloatVecType x;
  AutoPasLMP::FloatVecType v;

  int m = 1;
  x[0] = buf[m++];
  x[1] = buf[m++];
  x[2] = buf[m++];
  v[0] = buf[m++];
  v[1] = buf[m++];
  v[2] = buf[m++];
  tag[nlocal] = (tagint) ubuf(buf[m++]).i;
  type[nlocal] = (int) ubuf(buf[m++]).i;
  mask[nlocal] = (int) ubuf(buf[m++]).i;
  image[nlocal] = (imageint) ubuf(buf[m++]).i;

  // Always new particle from other process
  AutoPasLMP::ParticleType pi(x, v, static_cast<unsigned long>(nlocal),
                              static_cast<unsigned long>(type[nlocal]));
  lmp->autopas->add_particle(pi);

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->
          unpack_exchange(nlocal, &buf[m]);

  atom->nlocal++;
  return m;
}

int AtomVecAtomicAutopas::size_restart() {
  int i;

  int nlocal = atom->nlocal;
  int n = 11 * nlocal;

  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      for (i = 0; i < nlocal; i++)
        n += modify->fix[atom->extra_restart[iextra]]->size_restart(i);

  return n;
}

int AtomVecAtomicAutopas::pack_restart(int i, double *buf) {

  auto *pi = lmp->autopas->particle_by_index(i);
  auto &x = pi->getR();
  auto &v = pi->getV();

  int m = 1;
  buf[m++] = x[0];
  buf[m++] = x[1];
  buf[m++] = x[2];
  buf[m++] = ubuf(tag[i]).d;
  buf[m++] = ubuf(type[i]).d;
  buf[m++] = ubuf(mask[i]).d;
  buf[m++] = ubuf(image[i]).d;
  buf[m++] = v[0];
  buf[m++] = v[1];
  buf[m++] = v[2];

  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      m += modify->fix[atom->extra_restart[iextra]]->pack_restart(i, &buf[m]);

  buf[0] = m;
  return m;
}

int AtomVecAtomicAutopas::unpack_restart(double *buf) {

  // Only used from commands, thus use original arrays
  int nlocal = atom->nlocal;
  if (nlocal == nmax) {
    grow(0);
    if (atom->nextra_store)
      memory->grow(atom->extra, nmax, atom->nextra_store, "atom:extra");
  }

  int m = 1;
  x[nlocal][0] = buf[m++];
  x[nlocal][1] = buf[m++];
  x[nlocal][2] = buf[m++];
  tag[nlocal] = (tagint) ubuf(buf[m++]).i;
  type[nlocal] = (int) ubuf(buf[m++]).i;
  mask[nlocal] = (int) ubuf(buf[m++]).i;
  image[nlocal] = (imageint) ubuf(buf[m++]).i;
  v[nlocal][0] = buf[m++];
  v[nlocal][1] = buf[m++];
  v[nlocal][2] = buf[m++];

  double **extra = atom->extra;
  if (atom->nextra_store) {
    int size = static_cast<int> (buf[0]) - m;
    for (int i = 0; i < size; i++) extra[nlocal][i] = buf[m++];
  }

  atom->nlocal++;
  return m;
}

void AtomVecAtomicAutopas::create_atom(int itype, double *coord) {

  // Only used from commands, thus use original arrays
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  tag[nlocal] = 0;
  type[nlocal] = itype;
  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];
  mask[nlocal] = 1;
  image[nlocal] = ((imageint) IMGMAX << IMG2BITS) |
                  ((imageint) IMGMAX << IMGBITS) | IMGMAX;
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;

  atom->nlocal++;

  /*
  AutoPasLMP::FloatVecType pos{coord[0], coord[1], coord[2]};
  AutoPasLMP::FloatVecType  vel{0};
  unsigned long moleculeId = nlocal;
  unsigned long typeId = itype;

  lmp->autopas->addParticle(AutoPasLMP::ParticleType(pos, vel, moleculeId, typeId));*/

}

void
AtomVecAtomicAutopas::data_atom(double *coord, imageint imagetmp,
                                           char **values) {
  // Probably only used on startup / from commands -> use original arrays

  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  tag[nlocal] = utils::tnumeric(FLERR, values[0], true, lmp);
  type[nlocal] = utils::inumeric(FLERR, values[1], true, lmp);
  if (type[nlocal] <= 0 || type[nlocal] > atom->ntypes)
    error->one(FLERR, "Invalid atom type in Atoms section of data file");

  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];

  image[nlocal] = imagetmp;

  mask[nlocal] = 1;
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;

  atom->nlocal++;

}

void AtomVecAtomicAutopas::pack_data(double **buf) {

#pragma omp parallel default(none) shared(buf)
  for (auto iter = lmp->autopas->_autopas->begin(
      autopas::IteratorBehavior::ownedOnly); iter.isValid(); ++iter) {
    unsigned long idx = iter->getID();
    auto &x = iter->getR();
    buf[idx][0] = ubuf(tag[idx]).d;
    buf[idx][1] = ubuf(type[idx]).d;
    buf[idx][2] = x[0];
    buf[idx][3] = x[1];
    buf[idx][4] = x[2];
    buf[idx][5] = ubuf((image[idx] & IMGMASK) - IMGMAX).d;
    buf[idx][6] = ubuf((image[idx] >> IMGBITS & IMGMASK) - IMGMAX).d;
    buf[idx][7] = ubuf((image[idx] >> IMG2BITS) - IMGMAX).d;
  }
}

void
AtomVecAtomicAutopas::write_data(FILE *fp, int n, double **buf) {

  for (int i = 0; i < n; i++)
    fprintf(fp, TAGINT_FORMAT " %d %-1.16e %-1.16e %-1.16e %d %d %d\n",
            (tagint) ubuf(buf[i][0]).i, (int) ubuf(buf[i][1]).i,
            buf[i][2], buf[i][3], buf[i][4],
            (int) ubuf(buf[i][5]).i, (int) ubuf(buf[i][6]).i,
            (int) ubuf(buf[i][7]).i);
}

bigint AtomVecAtomicAutopas::memory_usage() {
  bigint bytes = 0;

  if (atom->memcheck("tag")) bytes += memory->usage(tag, nmax);
  if (atom->memcheck("type")) bytes += memory->usage(type, nmax);
  if (atom->memcheck("mask")) bytes += memory->usage(mask, nmax);
  if (atom->memcheck("image")) bytes += memory->usage(image, nmax);

  //TODO Memory usage from autopas
  //if (atom->memcheck("x")) bytes += memory->usage(x,nmax,3);
  //if (atom->memcheck("v")) bytes += memory->usage(v,nmax,3);
  //if (atom->memcheck("f")) bytes += memory->usage(f,nmax*comm->nthreads,3);

  return bytes;
}

int AtomVecAtomicAutopas::pack_exchange(
    const AutoPasLMP::ParticleType &p, double *buf) {
  auto idx = lmp->autopas->idx(p);
  auto &x = p.getR();
  auto &v = p.getV();

  int m = 1;
  buf[m++] = x[0];
  buf[m++] = x[1];
  buf[m++] = x[2];
  buf[m++] = v[0];
  buf[m++] = v[1];
  buf[m++] = v[2];
  buf[m++] = ubuf(tag[idx]).d;
  buf[m++] = ubuf(type[idx]).d;
  buf[m++] = ubuf(mask[idx]).d;
  buf[m++] = ubuf(image[idx]).d;

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->pack_exchange(idx, &buf[m]);

  buf[0] = m;
  return m;
}

