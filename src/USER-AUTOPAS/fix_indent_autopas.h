#ifdef FIX_CLASS

FixStyle(indent/autopas,FixIndentAutoPas)

#else

#ifndef LMP_FIX_INDENT_AUTOPAS_H
#define LMP_FIX_INDENT_AUTOPAS_H

#include "fix_indent.h"

namespace LAMMPS_NS {

class FixIndentAutoPas : public FixIndent {
public:
  using FixIndent::FixIndent;

  void post_force(int i) override;

};
};

#endif
#endif
