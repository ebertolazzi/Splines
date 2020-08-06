/**
 * PYTHON Wrapper for Splines
 * 
 * License MIT - See LICENSE file
 * 2019 Matteo Ragni, Claudio Kerov Ghiglianovich, Enrico Bertolazzi
 */

#include <string>

#include <Splines.hh>
#include "py-AkimaSplines.hh"

namespace pySpline {
  namespace py = pybind11;

  using Splines::Spline;
  using Splines::AkimaSpline;

  void python_register_akima_splines_class(py::module & m) {
    py::class_<AkimaSpline, Spline>(m, "AkimaSpline")
      .def(py::init<std::string const &>(), py::arg("name") = "AkimaSpline");
  }
  
}