#ifndef LMP_NEIGHBOR_AUTOPAS_H
#define LMP_NEIGHBOR_AUTOPAS_H

#include "neighbor.h"

namespace LAMMPS_NS {
class NeighborAutoPas : public Neighbor {
public:
  using Neighbor::Neighbor;

  int check_distance() override;

  int decide() override;

  void build(int i) override;
};
}

#endif
