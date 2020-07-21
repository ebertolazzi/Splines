/**
 * PYTHON Wrapper for Splines
 * 
 * License MIT - See LICENSE file
 * 2019 Matteo Ragni, Claudio Kerov Ghiglianovich, Enrico Bertolazzi
 */

#include <string>

#include <Splines.hh>
#include "py-ConstantSplines.hh"

namespace pySpline {
  namespace py = pybind11;

  using Splines::Spline;
  using Splines::ConstantSpline;

  void python_register_constant_spline_class(py::module & m) {
    py::class_<ConstantSpline, Spline>(m, "ConstantSpline")
      .def(py::init<std::string const &>(), py::arg("name") = "ConstantSpline");
  }
  
}