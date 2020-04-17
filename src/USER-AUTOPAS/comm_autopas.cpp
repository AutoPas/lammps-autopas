#include "autopas.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm_autopas.h"
#include "domain.h"

LAMMPS_NS::CommAutoPas::CommAutoPas(LAMMPS_NS::LAMMPS *lmp) : CommBrick(lmp) {

}

void LAMMPS_NS::CommAutoPas::forward_comm(int /*dummy*/) {

  int n;
  MPI_Request request;
  AtomVec *avec = atom->avec;

  double *buf;

  // exchange data with another proc
  // if other proc is self, just copy
  // if comm_x_only set, exchange or copy directly to x, don't unpack

  for (int iswap = 0; iswap < nswap; iswap++) {
    if (sendproc[iswap] != me) {
      if (comm_x_only) {
        if (size_forward_recv[iswap]) {
          auto x = lmp->autopas->particle_by_index(firstrecv[iswap])->getR();
          buf = x.data();
          MPI_Irecv(buf, size_forward_recv[iswap], MPI_DOUBLE,
                    recvproc[iswap], 0, world, &request);
        }
        n = avec->pack_comm(sendnum[iswap], sendlist[iswap],
                            buf_send, pbc_flag[iswap], pbc[iswap]);
        if (n) MPI_Send(buf_send, n, MPI_DOUBLE, sendproc[iswap], 0, world);
        if (size_forward_recv[iswap]) MPI_Wait(&request, MPI_STATUS_IGNORE);
      } else if (ghost_velocity) {
        if (size_forward_recv[iswap])
          MPI_Irecv(buf_recv, size_forward_recv[iswap], MPI_DOUBLE,
                    recvproc[iswap], 0, world, &request);
        n = avec->pack_comm_vel(sendnum[iswap], sendlist[iswap],
                                buf_send, pbc_flag[iswap], pbc[iswap]);
        if (n) MPI_Send(buf_send, n, MPI_DOUBLE, sendproc[iswap], 0, world);
        if (size_forward_recv[iswap]) MPI_Wait(&request, MPI_STATUS_IGNORE);
        avec->unpack_comm_vel(recvnum[iswap], firstrecv[iswap], buf_recv);
      } else {
        if (size_forward_recv[iswap])
          MPI_Irecv(buf_recv, size_forward_recv[iswap], MPI_DOUBLE,
                    recvproc[iswap], 0, world, &request);
        n = avec->pack_comm(sendnum[iswap], sendlist[iswap],
                            buf_send, pbc_flag[iswap], pbc[iswap]);
        if (n) MPI_Send(buf_send, n, MPI_DOUBLE, sendproc[iswap], 0, world);
        if (size_forward_recv[iswap]) MPI_Wait(&request, MPI_STATUS_IGNORE);
        avec->unpack_comm(recvnum[iswap], firstrecv[iswap], buf_recv);
      }

    } else {
      if (comm_x_only) {
        if (sendnum[iswap]) {
          auto x = lmp->autopas->particle_by_index(firstrecv[iswap])->getR();
          avec->pack_comm(sendnum[iswap], sendlist[iswap],
                          x.data(), pbc_flag[iswap], pbc[iswap]);
        }
      } else if (ghost_velocity) {
        avec->pack_comm_vel(sendnum[iswap], sendlist[iswap],
                            buf_send, pbc_flag[iswap], pbc[iswap]);
        avec->unpack_comm_vel(recvnum[iswap], firstrecv[iswap], buf_send);
      } else {
        avec->pack_comm(sendnum[iswap], sendlist[iswap],
                        buf_send, pbc_flag[iswap], pbc[iswap]);
        avec->unpack_comm(recvnum[iswap], firstrecv[iswap], buf_send);
      }
    }
  }
}

void LAMMPS_NS::CommAutoPas::reverse_comm() {
  int n;
  MPI_Request request;
  AtomVec *avec = atom->avec;
  double *buf;

  // exchange data with another proc
  // if other proc is self, just copy
  // if comm_f_only set, exchange or copy directly from f, don't pack

  for (int iswap = nswap - 1; iswap >= 0; iswap--) {
    if (sendproc[iswap] != me) {
      if (comm_f_only) {
        if (size_reverse_recv[iswap])
          MPI_Irecv(buf_recv, size_reverse_recv[iswap], MPI_DOUBLE,
                    sendproc[iswap], 0, world, &request);
        if (size_reverse_send[iswap]) {
          auto f = lmp->autopas->particle_by_index(firstrecv[iswap])->getF();
          buf = f.data();
          MPI_Send(buf, size_reverse_send[iswap], MPI_DOUBLE,
                   recvproc[iswap], 0, world);
        }
        if (size_reverse_recv[iswap]) MPI_Wait(&request, MPI_STATUS_IGNORE);
      } else {
        if (size_reverse_recv[iswap])
          MPI_Irecv(buf_recv, size_reverse_recv[iswap], MPI_DOUBLE,
                    sendproc[iswap], 0, world, &request);
        n = avec->pack_reverse(recvnum[iswap], firstrecv[iswap], buf_send);
        if (n) MPI_Send(buf_send, n, MPI_DOUBLE, recvproc[iswap], 0, world);
        if (size_reverse_recv[iswap]) MPI_Wait(&request, MPI_STATUS_IGNORE);
      }
      avec->unpack_reverse(sendnum[iswap], sendlist[iswap], buf_recv);

    } else {
      if (comm_f_only) {
        if (sendnum[iswap]) {
          auto f = lmp->autopas->particle_by_index(firstrecv[iswap])->getF();
          avec->unpack_reverse(sendnum[iswap], sendlist[iswap],
                               f.data());
        }
      } else {
        avec->pack_reverse(recvnum[iswap], firstrecv[iswap], buf_send);
        avec->unpack_reverse(sendnum[iswap], sendlist[iswap], buf_send);
      }
    }
  }
}

void LAMMPS_NS::CommAutoPas::exchange() {
  int i, m, nsend, nrecv, nrecv1, nrecv2, nlocal;
  double lo, hi, value;
  double *sublo, *subhi;
  MPI_Request request;
  AtomVec *avec = atom->avec;

  // clear global->local map for owned and ghost atoms
  // b/c atoms migrate to new procs in exchange() and
  //   new ghosts are created in borders()
  // map_set() is done at end of borders()
  // clear ghost count and any ghost bonus data internal to AtomVec

  if (map_style) atom->map_clear();
  atom->nghost = 0;
  atom->avec->clear_bonus();

  // insure send buf has extra space for a single atom
  // only need to reset if a fix can dynamically add to size of single atom

  if (maxexchange_fix_dynamic) {
    int bufextra_old = bufextra;
    init_exchange();
    if (bufextra > bufextra_old) grow_send(maxsend + bufextra, 2);
  }

  // subbox bounds for orthogonal or triclinic

  if (triclinic == 0) {
    sublo = domain->sublo;
    subhi = domain->subhi;
  } else {
    sublo = domain->sublo_lamda;
    subhi = domain->subhi_lamda;
  }

  // loop over dimensions

  int dimension = domain->dimension;

  for (int dim = 0; dim < dimension; dim++) {


    // fill buffer with atoms leaving my box, using < and >=
    // when atom is deleted, fill it in with last atom

    // TODO: Use AutoPas' leaving particle

    lo = sublo[dim];
    hi = subhi[dim];
    nlocal = atom->nlocal;
    i = nsend = 0;

    while (i < nlocal) {
      auto &x = lmp->autopas->particle_by_index(i)->getR();

      if (x[dim] < lo || x[dim] >= hi) {
        if (nsend > maxsend) grow_send(nsend, 1);
        nsend += avec->pack_exchange(i, &buf_send[nsend]);
        avec->copy(nlocal - 1, i, 1);
        nlocal--;
      } else i++;
    }
    atom->nlocal = nlocal;

    // send/recv atoms in both directions
    // send size of message first so receiver can realloc buf_recv if needed
    // if 1 proc in dimension, no send/recv
    //   set nrecv = 0 so buf_send atoms will be lost
    // if 2 procs in dimension, single send/recv
    // if more than 2 procs in dimension, send/recv to both neighbors

    if (procgrid[dim] == 1) nrecv = 0;
    else {
      MPI_Sendrecv(&nsend, 1, MPI_INT, procneigh[dim][0], 0,
                   &nrecv1, 1, MPI_INT, procneigh[dim][1], 0, world,
                   MPI_STATUS_IGNORE);
      nrecv = nrecv1;
      if (procgrid[dim] > 2) {
        MPI_Sendrecv(&nsend, 1, MPI_INT, procneigh[dim][1], 0,
                     &nrecv2, 1, MPI_INT, procneigh[dim][0], 0, world,
                     MPI_STATUS_IGNORE);
        nrecv += nrecv2;
      }
      if (nrecv > maxrecv) grow_recv(nrecv);

      MPI_Irecv(buf_recv, nrecv1, MPI_DOUBLE, procneigh[dim][1], 0,
                world, &request);
      MPI_Send(buf_send, nsend, MPI_DOUBLE, procneigh[dim][0], 0, world);
      MPI_Wait(&request, MPI_STATUS_IGNORE);

      if (procgrid[dim] > 2) {
        MPI_Irecv(&buf_recv[nrecv1], nrecv2, MPI_DOUBLE, procneigh[dim][0], 0,
                  world, &request);
        MPI_Send(buf_send, nsend, MPI_DOUBLE, procneigh[dim][1], 0, world);
        MPI_Wait(&request, MPI_STATUS_IGNORE);
      }
    }

    // check incoming atoms to see if they are in my box
    // if so, add to my list
    // box check is only for this dimension,
    //   atom may be passed to another proc in later dims

    m = 0;
    while (m < nrecv) {
      value = buf_recv[m + dim + 1];
      if (value >= lo && value < hi) m += avec->unpack_exchange(&buf_recv[m]);
      else m += static_cast<int> (buf_recv[m]);
    }
  }

  if (atom->firstgroupname) atom->first_reorder();
}

void LAMMPS_NS::CommAutoPas::borders() {

  int i, n, itype, iswap, dim, ineed, twoneed;
  int nsend, nrecv, sendflag, nfirst, nlast, ngroup;
  double lo, hi;
  int *type;
  double *buf, *mlo, *mhi;
  MPI_Request request;
  AtomVec *avec = atom->avec;

  // do swaps over all 3 dimensions

  iswap = 0;
  smax = rmax = 0;

  for (dim = 0; dim < 3; dim++) {
    nlast = 0;
    twoneed = 2 * maxneed[dim];
    for (ineed = 0; ineed < twoneed; ineed++) {

      // find atoms within slab boundaries lo/hi using <= and >=
      // check atoms between nfirst and nlast
      //   for first swaps in a dim, check owned and ghost
      //   for later swaps in a dim, only check newly arrived ghosts
      // store sent atom indices in sendlist for use in future timesteps

      if (mode == Comm::SINGLE) {
        lo = slablo[iswap];
        hi = slabhi[iswap];
      } else {
        type = atom->type;
        mlo = multilo[iswap];
        mhi = multihi[iswap];
      }
      if (ineed % 2 == 0) {
        nfirst = nlast;
        nlast = atom->nlocal + atom->nghost;
      }

      nsend = 0;

      // sendflag = 0 if I do not send on this swap
      // sendneed test indicates receiver no longer requires data
      // e.g. due to non-PBC or non-uniform sub-domains

      if (ineed / 2 >= sendneed[dim][ineed % 2]) sendflag = 0;
      else sendflag = 1;

      // find send atoms according to SINGLE vs MULTI
      // all atoms eligible versus only atoms in bordergroup
      // can only limit loop to bordergroup for first sends (ineed < 2)
      // on these sends, break loop in two: owned (in group) and ghost

      if (sendflag) {
        if (!bordergroup || ineed >= 2) {
          if (mode == Comm::SINGLE) {
            for (i = nfirst; i < nlast; i++) {
              auto &x = lmp->autopas->particle_by_index(i)->getR();
              if (x[dim] >= lo && x[dim] <= hi) {
                if (nsend == maxsendlist[iswap]) grow_list(iswap, nsend);
                sendlist[iswap][nsend++] = i;
              }
            }
          } else {
            for (i = nfirst; i < nlast; i++) {
              itype = type[i];
              auto &x = lmp->autopas->particle_by_index(i)->getR();
              if (x[dim] >= mlo[itype] && x[dim] <= mhi[itype]) {
                if (nsend == maxsendlist[iswap]) grow_list(iswap, nsend);
                sendlist[iswap][nsend++] = i;
              }
            }
          }

        } else {
          if (mode == Comm::SINGLE) {
            ngroup = atom->nfirst;
            for (i = 0; i < ngroup; i++) {
              auto &x = lmp->autopas->particle_by_index(i)->getR();
              if (x[dim] >= lo && x[dim] <= hi) {
                if (nsend == maxsendlist[iswap]) grow_list(iswap, nsend);
                sendlist[iswap][nsend++] = i;
              }
            }
            for (i = atom->nlocal; i < nlast; i++) {
              auto &x = lmp->autopas->particle_by_index(i)->getR();
              if (x[dim] >= lo && x[dim] <= hi) {
                if (nsend == maxsendlist[iswap]) grow_list(iswap, nsend);
                sendlist[iswap][nsend++] = i;
              }
            }
          } else {
            ngroup = atom->nfirst;
            for (i = 0; i < ngroup; i++) {
              itype = type[i];
              auto &x = lmp->autopas->particle_by_index(i)->getR();
              if (x[dim] >= mlo[itype] && x[dim] <= mhi[itype]) {
                if (nsend == maxsendlist[iswap]) grow_list(iswap, nsend);
                sendlist[iswap][nsend++] = i;
              }
            }
            for (i = atom->nlocal; i < nlast; i++) {
              itype = type[i];
              auto &x = lmp->autopas->particle_by_index(i)->getR();
              if (x[dim] >= mlo[itype] && x[dim] <= mhi[itype]) {
                if (nsend == maxsendlist[iswap]) grow_list(iswap, nsend);
                sendlist[iswap][nsend++] = i;
              }
            }
          }
        }
      }

      // pack up list of border atoms

      if (nsend * size_border > maxsend) grow_send(nsend * size_border, 0);
      if (ghost_velocity)
        n = avec->pack_border_vel(nsend, sendlist[iswap], buf_send,
                                  pbc_flag[iswap], pbc[iswap]);
      else
        n = avec->pack_border(nsend, sendlist[iswap], buf_send,
                              pbc_flag[iswap], pbc[iswap]);

      // swap atoms with other proc
      // no MPI calls except SendRecv if nsend/nrecv = 0
      // put incoming ghosts at end of my atom arrays
      // if swapping with self, simply copy, no messages

      if (sendproc[iswap] != me) {
        MPI_Sendrecv(&nsend, 1, MPI_INT, sendproc[iswap], 0,
                     &nrecv, 1, MPI_INT, recvproc[iswap], 0, world,
                     MPI_STATUS_IGNORE);
        if (nrecv * size_border > maxrecv) grow_recv(nrecv * size_border);
        if (nrecv)
          MPI_Irecv(buf_recv, nrecv * size_border, MPI_DOUBLE,
                    recvproc[iswap], 0, world, &request);
        if (n) MPI_Send(buf_send, n, MPI_DOUBLE, sendproc[iswap], 0, world);
        if (nrecv) MPI_Wait(&request, MPI_STATUS_IGNORE);
        buf = buf_recv;
      } else {
        nrecv = nsend;
        buf = buf_send;
      }

      // unpack buffer

      if (ghost_velocity)
        avec->unpack_border_vel(nrecv, atom->nlocal + atom->nghost, buf);
      else
        avec->unpack_border(nrecv, atom->nlocal + atom->nghost, buf);

      // set all pointers & counters

      smax = MAX(smax, nsend);
      rmax = MAX(rmax, nrecv);
      sendnum[iswap] = nsend;
      recvnum[iswap] = nrecv;
      size_forward_recv[iswap] = nrecv * size_forward;
      size_reverse_send[iswap] = nrecv * size_reverse;
      size_reverse_recv[iswap] = nsend * size_reverse;
      firstrecv[iswap] = atom->nlocal + atom->nghost;
      atom->nghost += nrecv;
      iswap++;
    }
  }

  // insure send/recv buffers are long enough for all forward & reverse comm

  int max = MAX(maxforward * smax, maxreverse * rmax);
  if (max > maxsend) grow_send(max, 0);
  max = MAX(maxforward * rmax, maxreverse * smax);
  if (max > maxrecv) grow_recv(max);

  // reset global->local map

  if (map_style) atom->map_set();
}