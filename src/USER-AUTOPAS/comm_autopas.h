#pragma once

#include "comm_brick.h"

namespace LAMMPS_NS {

class CommAutoPas : public CommBrick {
public:
  explicit CommAutoPas(class LAMMPS *);
  ~CommAutoPas() override = default;

  void forward_comm(int) override;

  void reverse_comm() override;

  void exchange() override;

  void borders() override;
};
}