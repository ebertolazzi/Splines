/**
 * PYTHON Wrapper for Splines
 * 
 * License MIT - See LICENSE file
 * 2019 Matteo Ragni, Claudio Kerov Ghiglianovich, Enrico Bertolazzi
 */

#include <string>

#include <Splines.hh>
#include "py-BesselSplines.hh"

namespace pySpline {
  namespace py = pybind11;

  using Splines::Spline;
  using Splines::BesselSpline;

  void python_register_bessel_splines_class(py::module & m) {
    py::class_<BesselSpline, Spline>(m, "BesselSpline")
      .def(py::init<std::string const &>(), py::arg("name") = "BesselSpline");
  }
  
}