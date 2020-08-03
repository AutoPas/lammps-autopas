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
#include "output_autopas.h"       // IWYU pragma: export


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
};

class AtomAutoPas : public Atom {
 public:
  AtomAutoPas(class LAMMPS *lmp) : Atom(lmp) {}
};

class CommAutoPas : public CommBrick {
 public:
  CommAutoPas(class LAMMPS *lmp) : CommBrick(lmp) {}
};

class CommTiledAutoPas : public CommTiled {
 public:
  CommTiledAutoPas(class LAMMPS *lmp) : CommTiled(lmp) {}
  CommTiledAutoPas(class LAMMPS *lmp, Comm *oldcomm) : CommTiled(lmp,oldcomm) {}
};

class DomainAutoPas : public Domain {
 public:
  DomainAutoPas(class LAMMPS *lmp) : Domain(lmp) {}
};

class NeighborAutoPas : public Neighbor {
 public:
  NeighborAutoPas(class LAMMPS *lmp) : Neighbor(lmp) {}
};

class MemoryAutoPas : public Memory {
 public:
  MemoryAutoPas(class LAMMPS *lmp) : Memory(lmp) {}
};

class ModifyAutoPas : public Modify {
 public:
  ModifyAutoPas(class LAMMPS *lmp) : Modify(lmp) {}
};

class OutputAutoPas : public Output {
 public:
  OutputAutoPas(class LAMMPS *lmp) : Output(lmp) {}
};

}

#endif
#endif
