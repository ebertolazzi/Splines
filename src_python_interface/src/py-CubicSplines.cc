/**
 * PYTHON Wrapper for Splines
 * 
 * License MIT - See LICENSE file
 * 2019 Matteo Ragni, Claudio Kerov Ghiglianovich, Enrico Bertolazzi
 */

#include <string>

#include <Splines.hh>
#include "py-CubicSplines.hh"

namespace pySpline {
  namespace py = pybind11;

  using Splines::real_type;
  using Splines::integer;

  using Splines::Spline;
  using Splines::CubicSplineBase;
  using Splines::CubicSpline;

  void python_register_cubic_splines_base_class(py::module & m) {
    py::class_<CubicSplineBase, Spline>(m, "CubicSplineBase")
      // .def(py::init<std::string const &>(), py::arg("name") = "CubicSplineBase")
      .def("copySpline", &CubicSplineBase::copySpline)
      .def("ypNode", &CubicSplineBase::ypNode);
  }
  
  void python_register_cubic_splines_class(module & m) {
    py::class_<CubicSpline, CubicSplineBase>(m, "CubicSpline")
      .def(py::init<std::string const &>(), py::arg("name") = "CubicSpline")
      .def("setInitialBC", &CubicSpline::setInitialBC)
      .def("setFinalBC", &CubicSpline::setFinalBC);
  }
  
}