/**
 * PYTHON Wrapper for Splines
 * 
 * License MIT - See LICENSE file
 * 2019 Matteo Ragni, Claudio Kerov Ghiglianovich, Enrico Bertolazzi
 */

#ifndef PY_PCHIP_SPLINES_HH
#define PY_PCHIP_SPLINES_HH

#include <pybind11/pybind11.h>

namespace pySpline {
  using pybind11::module;

  void python_register_pchip_splines_class(module & m);  
}

#endif /* PY_PCHIP_SPLINES_HH */