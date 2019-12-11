/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include <cmath>
#include <domain.h>
#include "pair_lj_cut_autopas.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"

#include "suffix.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairLJCutAutoPas::PairLJCutAutoPas(LAMMPS *lmp) :
    PairLJCut(lmp) {
  suffix_flag |= Suffix::AUTOPAS;
  respa_enable = 0;
  cut_respa = NULL;
  init_autopas();
}

/* ---------------------------------------------------------------------- */

void PairLJCutAutoPas::compute(int eflag, int vflag) {

  printf("Simulating computation in AutoPas\n");

  ev_init(eflag, vflag);

  const int nall = atom->nlocal + atom->nghost;
  const int inum = list->inum;


  // Copy to AutoPas
  for(int i=0; i<atom->nlocal; ++i){

    floatVecType pos{atom->x[i][0], atom->x[i][1], atom->x[i][2]};
    floatVecType vel{atom->v[i][0], atom->v[i][1], atom->v[i][2]};
    unsigned long moleculeId = i;
    unsigned long typeId = atom->type[i];

    auto particle = ParticleType(pos, vel, moleculeId, typeId);
    _autopas.addParticle(particle);

    // TODO_AUTOPAS Needs rvalue ref support in AutoPas
    // _autopas.addParticle(ParticleType(pos, vel, moleculeId, typeId))
  }

  // Copy from AutoPas
#ifdef AUTOPAS_OPENMP
#pragma omp parallel default(none)
#endif
  for(auto iter=_autopas.begin(autopas::IteratorBehavior::ownedOnly); iter.isValid(); ++iter){
    auto pos = iter->getR();
    auto vel = iter->getV();
    auto moleculeId = iter->getId();
  }

  for(int i=0; i<atom->nlocal; ++i){
    _autopas.i
  }
  _autopas.deleteAllParticles()

  /*
#if defined(_OPENMP)
#pragma omp parallel default(none) shared(eflag,vflag)
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, thr);

    if (evflag) {
      if (eflag) {
        if (force->newton_pair) eval<1,1,1>(ifrom, ito, thr);
        else eval<1,1,0>(ifrom, ito, thr);
      } else {
        if (force->newton_pair) eval<1,0,1>(ifrom, ito, thr);
        else eval<1,0,0>(ifrom, ito, thr);
      }
    } else {
      if (force->newton_pair) eval<0,0,1>(ifrom, ito, thr);
      else eval<0,0,0>(ifrom, ito, thr);
    }

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region
  */
}

/*
template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
void PairLJCutAutoPas::eval(int iifrom, int iito, ThrData * const thr)
{
  int i,j,ii,jj,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r2inv,r3inv,r6inv,forcelj,factor_lj;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;

  const dbl3_t * _noalias const x = (dbl3_t *) atom->x[0];
  dbl3_t * _noalias const f = (dbl3_t *) thr->get_f()[0];
  const int * _noalias const type = atom->type;
  const int nlocal = atom->nlocal;
  const double * _noalias const special_lj = force->special_lj;
  double fxtmp,fytmp,fztmp;

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = iifrom; ii < iito; ++ii) {

    i = ilist[ii];
    xtmp = x[i].x;
    ytmp = x[i].y;
    ztmp = x[i].z;
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    fxtmp=fytmp=fztmp=0.0;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j].x;
      dely = ytmp - x[j].y;
      delz = ztmp - x[j].z;
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r2inv = 1.0/rsq;
        r6inv = r2inv*r2inv*r2inv;
        r3inv = sqrt(r6inv);

        forcelj = r6inv * (lj1[itype][jtype]*r3inv - lj2[itype][jtype]);
        fpair = factor_lj*forcelj*r2inv;

        fxtmp += delx*fpair;
        fytmp += dely*fpair;
        fztmp += delz*fpair;
        if (NEWTON_PAIR || j < nlocal) {
          f[j].x -= delx*fpair;
          f[j].y -= dely*fpair;
          f[j].z -= delz*fpair;
        }

        if (EFLAG) {
          evdwl = r6inv*(lj3[itype][jtype]*r3inv-lj4[itype][jtype])
            - offset[itype][jtype];
          evdwl *= factor_lj;
        }

        if (EVFLAG) ev_tally_thr(this,i,j,nlocal,NEWTON_PAIR,
                                 evdwl,0.0,fpair,delx,dely,delz,thr);
      }
    }
    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;
  }
}*/

/* ---------------------------------------------------------------------- */

double PairLJCutAutoPas::memory_usage() {
  double bytes = 0; // memory_usage_thr(); //TODO_AUTOPAS Get memory usage from AutoPas
  bytes += PairLJCut::memory_usage();

  return bytes;
}

void PairLJCutAutoPas::init_autopas() {

  // Initialize particle properties
  for (int i = 1; i <= atom->ntypes; ++i) {
    _particlePropertiesLibrary->addType(
        i, epsilon[i][i], sigma[i][i], atom->mass[i], true
    );
  }

  // _autopas.setAllowedCellSizeFactors(*cellSizeFactors);
  //_autopas.setAllowedContainers(containerChoice);
  //_autopas.setAllowedDataLayouts(dataLayoutOptions);
  // _autopas.setAllowedNewton3Options(newton3Options);
  //_autopas.setAllowedTraversals(traversalOptions);

  floatVecType boxMax{}, boxMin{};
  std::copy(std::begin(domain->boxhi), std::end(domain->boxhi), boxMax.begin());
  std::copy(std::begin(domain->boxlo), std::end(domain->boxlo), boxMin.begin());
  _autopas.setBoxMax(boxMax);
  _autopas.setBoxMin(boxMin);

  _autopas.setCutoff(cut_global); // TODO_AUTOPAS Test: cut_global (PairLJCut) or cutforce (Pair)
  //_autopas.setNumSamples(tuningSamples);
  //_autopas.setSelectorStrategy(selectorStrategy);
  //_autopas.setTuningInterval(tuningInterval);
  //_autopas.setTuningStrategyOption(tuningStrategy);
  //_autopas.setVerletClusterSize(_config->verletClusterSize);
  //_autopas.setVerletRebuildFrequency(_config->verletRebuildFrequency);
  //_autopas.setVerletRebuildFrequency(verletRebuildFrequency);
  //_autopas.setVerletSkin(verletSkinRadius);
  //autopas::Logger::get()->set_level(logLevel);
  _autopas.init();

}
