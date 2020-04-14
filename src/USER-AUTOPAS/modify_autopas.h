#pragma once

#include "modify.h"

namespace LAMMPS_NS {

class ModifyAutoPas : public Modify {
 public:
  explicit ModifyAutoPas(class LAMMPS *);
  ~ModifyAutoPas() override = default;

  void setup(int) override;
  void setup_pre_exchange() override;
  void setup_pre_neighbor() override;
  void setup_pre_force(int) override;
  void setup_pre_reverse(int, int) override;
  void initial_integrate(int) override;
  void post_integrate() override;
  void pre_decide();
  void pre_exchange() override;
  void pre_neighbor() override;
  void pre_force(int) override;
  void pre_reverse(int,int) override;
  void post_force(int) override;
  void final_integrate() override;
  void end_of_step() override;
  double thermo_energy() override;
  void post_run() override;

  void setup_pre_force_respa(int, int) override;
  void initial_integrate_respa(int, int, int) override;
  void post_integrate_respa(int, int) override;
  void pre_force_respa(int, int, int) override;
  void post_force_respa(int, int, int) override;
  void final_integrate_respa(int, int) override;

  void min_pre_exchange() override;
  void min_pre_neighbor() override;
  void min_pre_force(int) override;
  void min_pre_reverse(int,int) override;
  void min_post_force(int) override;

  double min_energy(double *) override;
  void min_store() override;
  void min_step(double, double *) override;
  void min_clearstore() override;
  void min_pushstore() override;
  void min_popstore() override;
  double max_alpha(double *) override;
  int min_dof() override;
  int min_reset_ref() override;

 protected:

};

}
