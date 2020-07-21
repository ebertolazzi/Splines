/**
 * PYTHON Wrapper for Splines
 * 
 * License MIT - See LICENSE file
 * 2019 Matteo Ragni, Claudio Kerov Ghiglianovich, Enrico Bertolazzi
 */

#include <string>

#include <Splines.hh>
#include "py-PchipSplines.hh"

namespace pySpline {
  namespace py = pybind11;

  using Splines::Spline;
  using Splines::PchipSpline;

  void python_register_pchip_splines_class(py::module & m) {
    py::class_<PchipSpline, Spline>(m, "PchipSpline")
      .def(py::init<std::string const &>(), py::arg("name") = "PchipSpline");
  }
  
}