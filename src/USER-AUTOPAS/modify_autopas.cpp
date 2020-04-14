#include "modify_autopas.h"
#include "atom_autopas.h"
#include "update.h"
#include "fix.h"
#include "compute.h"
#include "autopas.h"

using namespace LAMMPS_NS;

ModifyAutoPas::ModifyAutoPas(LAMMPS *lmp) : Modify(lmp)
{
  // atomKK = (AtomKokkos *) atom;
}
/* ----------------------------------------------------------------------
   setup for run, calls setup() of all fixes and computes
   called from Verlet, RESPA, Min
------------------------------------------------------------------------- */

void ModifyAutoPas::setup(int vflag)
{
}

/* ----------------------------------------------------------------------
   setup pre_exchange call, only for fixes that define pre_exchange
   called from Verlet, RESPA, Min, and WriteRestart with whichflag = 0
------------------------------------------------------------------------- */

void ModifyAutoPas::setup_pre_exchange()
{
}

/* ----------------------------------------------------------------------
   setup pre_neighbor call, only for fixes that define pre_neighbor
   called from Verlet, RESPA
------------------------------------------------------------------------- */

void ModifyAutoPas::setup_pre_neighbor()
{
}

/* ----------------------------------------------------------------------
   setup pre_force call, only for fixes that define pre_force
   called from Verlet, RESPA, Min
------------------------------------------------------------------------- */

void ModifyAutoPas::setup_pre_force(int vflag)
{
}

/* ----------------------------------------------------------------------
   setup pre_reverse call, only for fixes that define pre_reverse
   called from Verlet, RESPA, Min
------------------------------------------------------------------------- */

void ModifyAutoPas::setup_pre_reverse(int eflag, int vflag)
{
}

/* ----------------------------------------------------------------------
   1st half of integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyAutoPas::initial_integrate(int vflag)
{
}

/* ----------------------------------------------------------------------
   post_integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyAutoPas::post_integrate()
{
}

/* ----------------------------------------------------------------------
   pre_exchange call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyAutoPas::pre_exchange()
{
}

/* ----------------------------------------------------------------------
   pre_neighbor call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyAutoPas::pre_neighbor()
{
}

/* ----------------------------------------------------------------------
   pre_force call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyAutoPas::pre_force(int vflag)
{
}

/* ----------------------------------------------------------------------
   pre_reverse call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyAutoPas::pre_reverse(int eflag, int vflag)
{
}

/* ----------------------------------------------------------------------
   post_force call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyAutoPas::post_force(int vflag)
{
}

/* ----------------------------------------------------------------------
   2nd half of integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyAutoPas::final_integrate()
{
}

/* ----------------------------------------------------------------------
   end-of-timestep call, only for relevant fixes
   only call fix->end_of_step() on timesteps that are multiples of nevery
------------------------------------------------------------------------- */

void ModifyAutoPas::end_of_step()
{
}

/* ----------------------------------------------------------------------
   thermo energy call, only for relevant fixes
   called by Thermo class
   compute_scalar() is fix call to return energy
------------------------------------------------------------------------- */

double ModifyAutoPas::thermo_energy()
{
  return 0;
}

/* ----------------------------------------------------------------------
   post_run call
------------------------------------------------------------------------- */

void ModifyAutoPas::post_run()
{
}

/* ----------------------------------------------------------------------
   setup rRESPA pre_force call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyAutoPas::setup_pre_force_respa(int vflag, int ilevel)
{
}

/* ----------------------------------------------------------------------
   1st half of rRESPA integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyAutoPas::initial_integrate_respa(int vflag, int ilevel, int iloop)
{
}

/* ----------------------------------------------------------------------
   rRESPA post_integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyAutoPas::post_integrate_respa(int ilevel, int iloop)
{
}

/* ----------------------------------------------------------------------
   rRESPA pre_force call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyAutoPas::pre_force_respa(int vflag, int ilevel, int iloop)
{
}

/* ----------------------------------------------------------------------
   rRESPA post_force call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyAutoPas::post_force_respa(int vflag, int ilevel, int iloop)
{
}

/* ----------------------------------------------------------------------
   2nd half of rRESPA integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyAutoPas::final_integrate_respa(int ilevel, int iloop)
{
}

/* ----------------------------------------------------------------------
   minimizer pre-exchange call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyAutoPas::min_pre_exchange()
{
}

/* ----------------------------------------------------------------------
   minimizer pre-neighbor call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyAutoPas::min_pre_neighbor()
{
}

/* ----------------------------------------------------------------------
   minimizer pre-force call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyAutoPas::min_pre_force(int vflag)
{
}

/* ----------------------------------------------------------------------
   minimizer pre-reverse call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyAutoPas::min_pre_reverse(int eflag, int vflag)
{
}

/* ----------------------------------------------------------------------
   minimizer force adjustment call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyAutoPas::min_post_force(int vflag)
{
}

/* ----------------------------------------------------------------------
   minimizer energy/force evaluation, only for relevant fixes
   return energy and forces on extra degrees of freedom
------------------------------------------------------------------------- */

double ModifyAutoPas::min_energy(double *fextra)
{
  return 0;
}

/* ----------------------------------------------------------------------
   store current state of extra dof, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyAutoPas::min_store()
{
}

/* ----------------------------------------------------------------------
   mange state of extra dof on a stack, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyAutoPas::min_clearstore()
{
}

void ModifyAutoPas::min_pushstore()
{
}

void ModifyAutoPas::min_popstore()
{
}

/* ----------------------------------------------------------------------
   displace extra dof along vector hextra, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyAutoPas::min_step(double alpha, double *hextra)
{
  
}

/* ----------------------------------------------------------------------
   compute max allowed step size along vector hextra, only for relevant fixes
------------------------------------------------------------------------- */

double ModifyAutoPas::max_alpha(double *hextra)
{
  return 0;
}

/* ----------------------------------------------------------------------
   extract extra dof for minimization, only for relevant fixes
------------------------------------------------------------------------- */

int ModifyAutoPas::min_dof()
{
  return 0;
}

/* ----------------------------------------------------------------------
   reset reference state of fix, only for relevant fixes
------------------------------------------------------------------------- */

int ModifyAutoPas::min_reset_ref()
{
  return 0;
}


