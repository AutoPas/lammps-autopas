#pragma once

#include "output.h"

namespace LAMMPS_NS {

class OutputAutoPas : public Output {

public:
  explicit OutputAutoPas(class LAMMPS *);
  ~OutputAutoPas() override = default;

  void setup(int memflag) override;

  void write(bigint) override;

  void write_dump(bigint) override;

  void write_restart(bigint) override;

};

}