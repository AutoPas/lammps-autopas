#pragma once

#include <autopas/molecularDynamics/MoleculeLJ.h>

namespace LAMMPS_NS {

template<typename floatType = double>
class MoleculeLJLammps : public autopas::MoleculeLJ<floatType> {
public:

  /**
   * Constructor of lennard jones molecule with initialization of typeID and localID.
   * @param pos Position of the molecule.
   * @param v Velocity of the molecule.
   * @param moleculeId Global Id of the molecule.
   * @param localId Local Id of the molecule.
   * @param typeId TypeId of the molecule.
   */
  MoleculeLJLammps(std::array<floatType, 3> pos, std::array<floatType, 3> v,
                   unsigned long moleculeId, int localId,
                   unsigned long typeId = 0)
      : autopas::MoleculeLJ<floatType>(pos, v, moleculeId, typeId),
        _localId(localId) {}

  [[nodiscard]] int getLocalID() {
    return _localId;
  }

  void setLocalID(int localId) {
    _localId = localId;
  }

private:
  int _localId = -1;
};

}
