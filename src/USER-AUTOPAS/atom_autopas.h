#include "atom.h"

#ifndef LMP_ATOM_AUTOPAS_H
#define LMP_ATOM_AUTOPAS_H

namespace LAMMPS_NS {

class AtomAutoPas : public Atom {
public:
  using Atom::Atom;

  void sort() override;
};

}
#endif
