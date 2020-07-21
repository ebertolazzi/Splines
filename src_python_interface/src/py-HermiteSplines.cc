/**
 * PYTHON Wrapper for Splines
 * 
 * License MIT - See LICENSE file
 * 2019 Matteo Ragni, Claudio Kerov Ghiglianovich, Enrico Bertolazzi
 */

#include <string>

#include <Splines.hh>
#include "py-HermiteSplines.hh"

namespace pySpline {
  namespace py = pybind11;

  using Splines::Spline;
  using Splines::HermiteSpline;

  void python_register_hermite_splines_class(py::module & m) {
    py::class_<HermiteSpline, Spline>(m, "HermiteSpline")
      .def(py::init<std::string const &>(), py::arg("name") = "HermiteSpline");
  }
  
}