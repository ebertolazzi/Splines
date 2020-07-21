/**
 * PYTHON Wrapper for Splines
 * 
 * License MIT - See LICENSE file
 * 2019 Matteo Ragni, Claudio Kerov Ghiglianovich, Enrico Bertolazzi
 */

#ifndef PY_BESSEL_SPLINES_HH
#define PY_BESSEL_SPLINES_HH

#include <pybind11/pybind11.h>

namespace pySpline {
  using pybind11::module;

  void python_register_bessel_splines_class(module & m);  
}

#endif /* PY_BESSEL_SPLINES_HH */