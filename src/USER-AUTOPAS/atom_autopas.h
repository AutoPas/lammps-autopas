#include "atom.h"

#ifndef LMP_ATOM_AUTOPAS_H
#define LMP_ATOM_AUTOPAS_H

namespace LAMMPS_NS {

class AtomAutoPas : public Atom {
public:
  explicit AtomAutoPas(class LAMMPS *);

  void sort() override;
};

}
#endif
