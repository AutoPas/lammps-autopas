#include "domain_autopas.h"

using namespace LAMMPS_NS;

DomainAutoPas::DomainAutoPas(LAMMPS_NS::LAMMPS *lmp) : Domain(lmp) {

}

void DomainAutoPas::pbc() {
  Domain::pbc();
}

void DomainAutoPas::image_check() {
  Domain::image_check();
}

void DomainAutoPas::box_too_small_check() {
  Domain::box_too_small_check();
}

int DomainAutoPas::closest_image(int i, int i1) {
  return Domain::closest_image(i, i1);
}

int DomainAutoPas::closest_image(const double *const pDouble, int i) {
  return Domain::closest_image(pDouble, i);
}

void DomainAutoPas::lamda2x(int i) {
  Domain::lamda2x(i);
}

void DomainAutoPas::x2lamda(int i) {
  Domain::x2lamda(i);
}
