#include "atom.h"
// #include "autopas_type.h"  // TODO_AP: Same types as lammps

#pragma once

namespace LAMMPS_NS {

class AtomAutoPas : public Atom {
 public:
  AtomAutoPas(class LAMMPS *);
  ~AtomAutoPas();

  virtual void allocate_type_arrays();
  void sync(const ExecutionSpace space, unsigned int mask);
  void modified(const ExecutionSpace space, unsigned int mask);
  void sync_overlapping_device(const ExecutionSpace space, unsigned int mask);
  virtual void sort();
  virtual void grow(unsigned int mask);
  int add_custom(const char *, int);
  void remove_custom(int, int);
  virtual void deallocate_topology();
  void sync_modify(ExecutionSpace, unsigned int, unsigned int);
};

}