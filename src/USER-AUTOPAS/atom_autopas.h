#include "atom.h"
// #include "autopas_type.h"  // TODO_AP: Same types as lammps

#pragma once

namespace LAMMPS_NS {

class AtomAutoPas : public Atom {
 public:
  explicit AtomAutoPas(class LAMMPS *);
  ~AtomAutoPas() override;

  //void allocate_type_arrays() override;
  //void sync(const ExecutionSpace space, unsigned int mask);
  //void modified(const ExecutionSpace space, unsigned int mask);
  //void sync_overlapping_device(const ExecutionSpace space, unsigned int mask);
  void sort() override;
  //virtual void grow(unsigned int mask);
  //int add_custom(const char *, int) override;
  //void remove_custom(int, int) override;
  //virtual void deallocate_topology();
  //void sync_modify(ExecutionSpace, unsigned int, unsigned int) override;
};

}