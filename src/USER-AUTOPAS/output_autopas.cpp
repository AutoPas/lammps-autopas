#include "output_autopas.h"

#include "autopas.h"

using namespace LAMMPS_NS;

OutputAutoPas::OutputAutoPas(LAMMPS *lmp) : Output(lmp) {

}

void OutputAutoPas::setup(int memflag) {
  lmp->autopas->copy_back();
  Output::setup(memflag);
}

void OutputAutoPas::write(bigint ntimestep) {
  lmp->autopas->copy_back();
  Output::write(ntimestep);
}

void OutputAutoPas::write_dump(bigint ntimestep) {
  lmp->autopas->copy_back();
  Output::write_dump(ntimestep);
}

void OutputAutoPas::write_restart(bigint ntimestep) {
  lmp->autopas->copy_back();
  Output::write_restart(ntimestep);
}
