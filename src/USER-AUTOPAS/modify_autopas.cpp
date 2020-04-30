#include "modify_autopas.h"

#include "atom_autopas.h"
#include "autopas.h"
#include "compute.h"
#include "fix.h"
#include "update.h"

using namespace LAMMPS_NS;

ModifyAutoPas::ModifyAutoPas(LAMMPS *lmp) : Modify(lmp) {
  // atomKK = (AtomKokkos *) atom;
}

/* ----------------------------------------------------------------------
   setup for run, calls setup() of all fixes and computes
   called from Verlet, RESPA, Min
------------------------------------------------------------------------- */

void ModifyAutoPas::setup(int vflag) {
  throw "Not implemented";
}

/* ----------------------------------------------------------------------
   setup pre_exchange call, only for fixes that define pre_exchange
   called from Verlet, RESPA, Min, and WriteRestart with whichflag = 0
------------------------------------------------------------------------- */

void ModifyAutoPas::setup_pre_exchange() {
  throw "Not implemented";
}

/* ----------------------------------------------------------------------
   setup pre_neighbor call, only for fixes that define pre_neighbor
   called from Verlet, RESPA
------------------------------------------------------------------------- */

void ModifyAutoPas::setup_pre_neighbor() {
  throw "Not implemented";
}

/* ----------------------------------------------------------------------
   setup pre_force call, only for fixes that define pre_force
   called from Verlet, RESPA, Min
------------------------------------------------------------------------- */

void ModifyAutoPas::setup_pre_force(int vflag) {
  throw "Not implemented";
}

/* ----------------------------------------------------------------------
   setup pre_reverse call, only for fixes that define pre_reverse
   called from Verlet, RESPA, Min
------------------------------------------------------------------------- */

void ModifyAutoPas::setup_pre_reverse(int eflag, int vflag) {
  throw "Not implemented";
}

/* ----------------------------------------------------------------------
   1st half of integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyAutoPas::initial_integrate(int vflag) {
  throw "Not implemented";
}

/* ----------------------------------------------------------------------
   post_integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyAutoPas::post_integrate() {
  throw "Not implemented";
}

/* ----------------------------------------------------------------------
   pre_exchange call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyAutoPas::pre_exchange() {
  throw "Not implemented";
}

/* ----------------------------------------------------------------------
   pre_neighbor call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyAutoPas::pre_neighbor() {
  throw "Not implemented";
}

/* ----------------------------------------------------------------------
   pre_force call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyAutoPas::pre_force(int vflag) {
  throw "Not implemented";
}

/* ----------------------------------------------------------------------
   pre_reverse call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyAutoPas::pre_reverse(int eflag, int vflag) {
  throw "Not implemented";
}

/* ----------------------------------------------------------------------
   post_force call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyAutoPas::post_force(int vflag) {
  throw "Not implemented";
}

/* ----------------------------------------------------------------------
   2nd half of integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyAutoPas::final_integrate() {
  throw "Not implemented";
}

/* ----------------------------------------------------------------------
   end-of-timestep call, only for relevant fixes
   only call fix->end_of_step() on timesteps that are multiples of nevery
------------------------------------------------------------------------- */

void ModifyAutoPas::end_of_step() {
  throw "Not implemented";
}

/* ----------------------------------------------------------------------
   thermo energy call, only for relevant fixes
   called by Thermo class
   compute_scalar() is fix call to return energy
------------------------------------------------------------------------- */

double ModifyAutoPas::thermo_energy() {
  throw "Not implemented";
  return 0;
}

/* ----------------------------------------------------------------------
   post_run call
------------------------------------------------------------------------- */

void ModifyAutoPas::post_run() {
  throw "Not implemented";
}

/* ----------------------------------------------------------------------
   setup rRESPA pre_force call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyAutoPas::setup_pre_force_respa(int vflag, int ilevel) {
  throw "Not implemented";
}

/* ----------------------------------------------------------------------
   1st half of rRESPA integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyAutoPas::initial_integrate_respa(int vflag, int ilevel, int iloop) {
  throw "Not implemented";
}

/* ----------------------------------------------------------------------
   rRESPA post_integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyAutoPas::post_integrate_respa(int ilevel, int iloop) {
  throw "Not implemented";
}

/* ----------------------------------------------------------------------
   rRESPA pre_force call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyAutoPas::pre_force_respa(int vflag, int ilevel, int iloop) {
  throw "Not implemented";
}

/* ----------------------------------------------------------------------
   rRESPA post_force call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyAutoPas::post_force_respa(int vflag, int ilevel, int iloop) {
  throw "Not implemented";
}

/* ----------------------------------------------------------------------
   2nd half of rRESPA integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyAutoPas::final_integrate_respa(int ilevel, int iloop) {
  throw "Not implemented";
}

/* ----------------------------------------------------------------------
   minimizer pre-exchange call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyAutoPas::min_pre_exchange() {
  throw "Not implemented";
}

/* ----------------------------------------------------------------------
   minimizer pre-neighbor call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyAutoPas::min_pre_neighbor() {
  throw "Not implemented";
}

/* ----------------------------------------------------------------------
   minimizer pre-force call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyAutoPas::min_pre_force(int vflag) {
  throw "Not implemented";
}

/* ----------------------------------------------------------------------
   minimizer pre-reverse call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyAutoPas::min_pre_reverse(int eflag, int vflag) {
  throw "Not implemented";
}

/* ----------------------------------------------------------------------
   minimizer force adjustment call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyAutoPas::min_post_force(int vflag) {
  throw "Not implemented";
}

/* ----------------------------------------------------------------------
   minimizer energy/force evaluation, only for relevant fixes
   return energy and forces on extra degrees of freedom
------------------------------------------------------------------------- */

double ModifyAutoPas::min_energy(double *fextra) {
  throw "Not implemented";
  return 0;
}

/* ----------------------------------------------------------------------
   store current state of extra dof, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyAutoPas::min_store() {
  throw "Not implemented";
}

/* ----------------------------------------------------------------------
   mange state of extra dof on a stack, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyAutoPas::min_clearstore() {
  throw "Not implemented";
}

void ModifyAutoPas::min_pushstore() {
  throw "Not implemented";
}

void ModifyAutoPas::min_popstore() {
  throw "Not implemented";
}

/* ----------------------------------------------------------------------
   displace extra dof along vector hextra, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyAutoPas::min_step(double alpha, double *hextra) {
  throw "Not implemented";
}

/* ----------------------------------------------------------------------
   compute max allowed step size along vector hextra, only for relevant fixes
------------------------------------------------------------------------- */

double ModifyAutoPas::max_alpha(double *hextra) {
  throw "Not implemented";
  return 0;
}

/* ----------------------------------------------------------------------
   extract extra dof for minimization, only for relevant fixes
------------------------------------------------------------------------- */

int ModifyAutoPas::min_dof() {
  throw "Not implemented";
  return 0;
}

/* ----------------------------------------------------------------------
   reset reference state of fix, only for relevant fixes
------------------------------------------------------------------------- */

int ModifyAutoPas::min_reset_ref() {
  throw "Not implemented";
  return 0;
}


