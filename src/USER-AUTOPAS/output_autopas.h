#ifndef LMP_OUTPUT_AUTOPAS_H
#define LMP_OUTPUT_AUTOPAS_H

#include "output.h"

namespace LAMMPS_NS {

class OutputAutoPas : public Output {
public:
  explicit OutputAutoPas(class LAMMPS *);

  void setup(int memflag) override;

  void write(bigint) override;

  void write_dump(bigint) override;

  void write_restart(bigint) override;

};

}

#endif
