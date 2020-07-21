/**
 * PYTHON Wrapper for Splines
 * 
 * License MIT - See LICENSE file
 * 2019 Matteo Ragni, Claudio Kerov Ghiglianovich, Enrico Bertolazzi
 */

#ifndef PY_AKIMA_SPLINES_HH
#define PY_AKIMA_SPLINES_HH

#include <pybind11/pybind11.h>

namespace pySpline {
  using pybind11::module;

  void python_register_akima_splines_class(module & m);  
}

#endif /* PY_AKIMA_SPLINES_HH */