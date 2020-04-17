#pragma once

#include "neighbor.h"

namespace LAMMPS_NS {
class NeighborAutoPas : public Neighbor {
public:
  explicit NeighborAutoPas(class LAMMPS *);

  ~NeighborAutoPas() override = default;

public:
  int check_distance() override;

  void build(int i) override;
};
}