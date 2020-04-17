#include "comm_autopas.h"

LAMMPS_NS::CommAutoPas::CommAutoPas(LAMMPS_NS::LAMMPS *lmp) : CommBrick(lmp){

}

void LAMMPS_NS::CommAutoPas::forward_comm(int dummy) {
  CommBrick::forward_comm(dummy);
}

void LAMMPS_NS::CommAutoPas::reverse_comm() {
  CommBrick::reverse_comm();
}

void LAMMPS_NS::CommAutoPas::exchange() {
  CommBrick::exchange();
}

void LAMMPS_NS::CommAutoPas::borders() {
  CommBrick::borders();
}
