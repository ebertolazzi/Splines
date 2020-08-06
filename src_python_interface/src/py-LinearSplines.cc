/**
 * PYTHON Wrapper for Splines
 * 
 * License MIT - See LICENSE file
 * 2019 Matteo Ragni, Claudio Kerov Ghiglianovich, Enrico Bertolazzi
 */

#include <string>

#include <Splines.hh>
#include "py-LinearSplines.hh"

namespace pySpline {
  namespace py = pybind11;

  using Splines::Spline;
  using Splines::LinearSpline;

  void python_register_linear_splines_class(py::module & m) {
    py::class_<LinearSpline, Spline>(m, "LinearSpline")
      .def(py::init<std::string const &>(), py::arg("name") = "LinearSpline");
  }
  
}