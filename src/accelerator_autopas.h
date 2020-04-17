#ifndef LMP_ACCELERATOR_AUTOPAS_H
#define LMP_ACCELERATOR_AUTOPAS_H

// true interface to AUTOPAS
// used when AUTOPAS is installed

#ifdef LMP_AUTOPAS

#include "autopas.h"             // IWYU pragma: export
#include "atom_autopas.h"        // IWYU pragma: export
#include "comm_autopas.h"         // IWYU pragma: export
// #include "comm_tiled_autopas.h"   // IWYU pragma: export
#include "domain_autopas.h"       // IWYU pragma: export
#include "neighbor_autopas.h"     // IWYU pragma: export
#include "memory_autopas.h"       // IWYU pragma: export
#include "modify_autopas.h"       // IWYU pragma: export

#else

// dummy interface to AUTOPAS
// needed for compiling when AUTOPAS is not installed

#include "atom.h"
#include "comm_brick.h"
#include "comm_tiled.h"
#include "domain.h"
#include "neighbor.h"
#include "memory.h"
#include "modify.h"

namespace LAMMPS_NS {

class AutoPasLMP {
 public:
  int autopas_exists;

  AutoPasLMP(class LAMMPS *, int, char **) {autopas_exists = 0;}
  ~AutoPasLMP() {}
  void accelerator(int, char **) {}
  int neigh_list_autopas(int) {return 0;}
  int neigh_count(int) {return 0;}
};

class AtomAutoPas : public Atom {
 public:
  tagint **k_special;
  AtomAutoPas(class LAMMPS *lmp) : Atom(lmp) {}
  ~AtomAutoPas() {}
  void sync(const ExecutionSpace /*space*/, unsigned int /*mask*/) {}
  void modified(const ExecutionSpace /*space*/, unsigned int /*mask*/) {}
};

class CommAutoPas : public CommBrick {
 public:
  CommAutoPas(class LAMMPS *lmp) : CommBrick(lmp) {}
  ~CommAutoPas() {}
};

class CommTiledAutoPas : public CommTiled {
 public:
  CommTiledAutoPas(class LAMMPS *lmp) : CommTiled(lmp) {}
  CommTiledAutoPas(class LAMMPS *lmp, Comm *oldcomm) : CommTiled(lmp,oldcomm) {}
  ~CommTiledAutoPas() {}
};

class DomainAutoPas : public Domain {
 public:
  DomainAutoPas(class LAMMPS *lmp) : Domain(lmp) {}
  ~DomainAutoPas() {}
};

class NeighborAutoPas : public Neighbor {
 public:
  NeighborAutoPas(class LAMMPS *lmp) : Neighbor(lmp) {}
  ~NeighborAutoPas() {}
};

class MemoryAutoPas : public Memory {
 public:
  MemoryAutoPas(class LAMMPS *lmp) : Memory(lmp) {}
  ~MemoryAutoPas() {}
  void grow_autopas(tagint **, tagint **, int, int, const char*) {}
};

class ModifyAutoPas : public Modify {
 public:
  ModifyAutoPas(class LAMMPS *lmp) : Modify(lmp) {}
  ~ModifyAutoPas() {}
};

}

#endif
#endif
