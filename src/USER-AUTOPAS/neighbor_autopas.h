#ifndef LMP_NEIGHBOR_AUTOPAS_H
#define LMP_NEIGHBOR_AUTOPAS_H

#include "neighbor.h"

namespace LAMMPS_NS {
class NeighborAutoPas : public Neighbor {
public:
  explicit NeighborAutoPas(class LAMMPS *);

  int check_distance() override;

  void build(int i) override;
};
}

#endif
