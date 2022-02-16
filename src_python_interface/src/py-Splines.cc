/**
 * PYTHON Wrapper for Splines
 *
 * License MIT - See LICENSE file
 * 2019 Matteo Ragni, Claudio Kerov Ghiglianovich, Enrico Bertolazzi
 */

#include "py-Splines.hh"
#include "py-CubicSplines.hh"
#include "py-ConstantSplines.hh"
#include "py-AkimaSplines.hh"
#include "py-BesselSplines.hh"
#include "py-HermiteSplines.hh"
#include "py-LinearSplines.hh"
#include "py-PchipSplines.hh"
#include "py-QuinticSplines.hh"

#include <pybind11/stl.h>
#include <tuple>

namespace pySpline {
  using Splines::real_type;
  using Splines::integer;

  using Splines::Spline;
  using Splines::SplineType1D;
  using Splines::SplineType2D;
  using Splines::CUBIC_SPLINE_TYPE_BC;
  using Splines::QUINTIC_SPLINE_TYPE;

  namespace py = pybind11;

  void python_register_splines_class(py::module & m) {
    py::enum_<SplineType1D>(m, "SplineType1D")
      .value("CONSTANT_TYPE", SplineType1D::CONSTANT_TYPE)
      .value("LINEAR_TYPE", SplineType1D::LINEAR_TYPE)
      .value("CUBIC_TYPE", SplineType1D::CUBIC_TYPE)
      .value("AKIMA_TYPE", SplineType1D::AKIMA_TYPE)
      .value("BESSEL_TYPE", SplineType1D::BESSEL_TYPE)
      .value("PCHIP_TYPE", SplineType1D::PCHIP_TYPE)
      .value("QUINTIC_TYPE", SplineType1D::QUINTIC_TYPE)
      .value("HERMITE_TYPE", SplineType1D::HERMITE_TYPE)
      .value("SPLINE_SET_TYPE", SplineType1D::SPLINE_SET_TYPE)
      .value("SPLINE_VEC_TYPE", SplineType1D::SPLINE_VEC_TYPE)
      .export_values();

    py::enum_<SplineType2D>(m, "SplineType2D")
      .value("BILINEAR_TYPE", SplineType2D::BILINEAR_TYPE)
      .value("BICUBIC_TYPE", SplineType2D::BICUBIC_TYPE)
      .value("BIQUINTIC_TYPE", SplineType2D::BIQUINTIC_TYPE)
      .value("AKIMA2D_TYPE", SplineType2D::AKIMA2D_TYPE)
      .export_values();

    py::enum_<CUBIC_SPLINE_TYPE_BC>(m, "CUBIC_SPLINE_TYPE_BC")
      .value("EXTRAPOLATE_BC", CUBIC_SPLINE_TYPE_BC::EXTRAPOLATE_BC)
      .value("NATURAL_BC", CUBIC_SPLINE_TYPE_BC::NATURAL_BC)
      .value("PARABOLIC_RUNOUT_BC", CUBIC_SPLINE_TYPE_BC::PARABOLIC_RUNOUT_BC)
      .value("NOT_A_KNOT", CUBIC_SPLINE_TYPE_BC::NOT_A_KNOT)
      .export_values();

    py::enum_<QUINTIC_SPLINE_TYPE>(m, "QUINTIC_SPLINE_TYPE")
      .value("CUBIC_QUINTIC", QUINTIC_SPLINE_TYPE::CUBIC_QUINTIC)
      .value("PCHIP_QUINTIC", QUINTIC_SPLINE_TYPE::PCHIP_QUINTIC)
      .value("AKIMA_QUINTIC", QUINTIC_SPLINE_TYPE::AKIMA_QUINTIC)
      .value("BESSEL_QUINTIC", QUINTIC_SPLINE_TYPE::BESSEL_QUINTIC)
      .export_values();

    py::class_<Spline, PythonicSpline>(m, "Spline")
      .def(py::init<std::string const &>(), py::arg("name") = "Spline")
      .def("name", &Spline::name)
      .def("is_closed", &Spline::is_closed)
      .def("make_closed", &Spline::make_closed)
      .def("make_opened", &Spline::make_opened)
      .def("is_bounded", &Spline::is_bounded)
      .def("make_unbounded", &Spline::make_unbounded)
      .def("make_bounded", &Spline::make_bounded)
      .def("num_points", &Spline::num_points)
      .def("xNode", &Spline::xNode)
      .def("yNode", &Spline::yNode)
      .def("xBegin", &Spline::xBegin)
      .def("yBegin", &Spline::yBegin)
      .def("xEnd", &Spline::xEnd)
      .def("reserve", &Spline::reserve)
      .def("pushBack", &Spline::pushBack)
      .def("dropBack", &Spline::dropBack)
      .def("build", py::overload_cast<>(&Spline::build))
      .def("build", py::overload_cast< std::vector<real_type> const &, std::vector<real_type> const & >(&Spline::build))
      .def("clear", &Spline::clear)
      .def("xMin", &Spline::xMin)
      .def("xMax", &Spline::xMax)
      .def("yMin", &Spline::yMin)
      .def("yMax", &Spline::yMax)
      .def("type", &Spline::type)
      .def("setOrigin", &Spline::setOrigin)
      .def("setRange", &Spline::setRange)
      .def("__call__", &Spline::operator())
      .def("D", &Spline::D)
      .def("DD", &Spline::DD)
      .def("DDD", &Spline::DDD)
      .def("DDDD", &Spline::DDDD)
      .def("DDDDD", &Spline::DDDDD)
      .def("eval", &Spline::eval)
      .def("eval_D", &Spline::eval_D)
      .def("eval_DD", &Spline::eval_DD)
      .def("eval_DDD", &Spline::eval_DDD)
      .def("eval_DDDD", &Spline::eval_DDDD)
      .def("eval_DDDDD", &Spline::eval_DDDDD)
      .def("coeffs", [](Spline & self, bool transpose)
        -> std::tuple< std::vector<real_type>, std::vector<real_type> > {
        integer n_nodes = self.num_points() - 1;
        integer n_order = self.order();
        std::vector<real_type> cfs(n_order * n_nodes);
        std::vector<real_type> nodes(n_nodes);
        self.coeffs(cfs.data(), nodes.data(), transpose);
        return std::make_tuple(cfs, nodes);
      }, py::arg("transpose") = false)
      .def("order", &Spline::order)
      .def("typre_name", &Spline::type_name)
      .def("__str__", [](const Spline & self)
        -> std::string {
        std::ostringstream str;
        self.info(str);
        return str.str();
      })
      .def("writeToString", [](const Spline & self)
        -> std::string {
        std::ostringstream str;
        self.writeToStream(str);
        return str.str();
      })
      .def("dump", [](const Spline & self, integer nintervals, std::string header)
        -> std::string {
        std::ostringstream str;
        self.dump(str, nintervals, header.c_str());
        return str.str();
      }, py::arg("nintervals"), py::arg("header") = "x\ty");
  }

  void python_register_hermite_functions(module & m) {
    m.def("Hermite3", [](real_type const x, real_type const H)
    -> std::vector<real_type> {
      std::vector<real_type> result(4);
      Splines::Hermite3(x, H, result.data());
      return result;
    })
    .def("Hermite3_D", [](real_type const x, real_type const H)
    -> std::vector<real_type> {
      std::vector<real_type> result(4);
      Splines::Hermite3_D(x, H, result.data());
      return result;
    })
    .def("Hermite3_DD", [](real_type const x, real_type const H)
    -> std::vector<real_type> {
      std::vector<real_type> result(4);
      Splines::Hermite3_DD(x, H, result.data());
      return result;
    })
    .def("Hermite3_DDD", [](real_type const x, real_type const H)
    -> std::vector<real_type> {
      std::vector<real_type> result(4);
      Splines::Hermite3_DDD(x, H, result.data());
      return result;
    })
    .def("Hermite5", [](real_type const x, real_type const H)
    -> std::vector<real_type> {
      std::vector<real_type> result(6);
      Splines::Hermite5(x, H, result.data());
      return result;
    })
    .def("Hermite5_D", [](real_type const x, real_type const H)
    -> std::vector<real_type> {
      std::vector<real_type> result(6);
      Splines::Hermite5_D(x, H, result.data());
      return result;
    })
    .def("Hermite5_DD", [](real_type const x, real_type const H)
    -> std::vector<real_type> {
      std::vector<real_type> result(6);
      Splines::Hermite5_DD(x, H, result.data());
      return result;
    })
    .def("Hermite5_DDD", [](real_type const x, real_type const H)
    -> std::vector<real_type> {
      std::vector<real_type> result(6);
      Splines::Hermite5_DDD(x, H, result.data());
      return result;
    })
    .def("Hermite5_DDDD", [](real_type const x, real_type const H)
    -> std::vector<real_type> {
      std::vector<real_type> result(6);
      Splines::Hermite5_DDDD(x, H, result.data());
      return result;
    })
    .def("Hermite5_DDDDD", [](real_type const x, real_type const H)
    -> std::vector<real_type> {
      std::vector<real_type> result(6);
      Splines::Hermite5_DDDDD(x, H, result.data());
      return result;
    });
  }

}

PYBIND11_MODULE(Splines, m) {
  pySpline::python_register_splines_class(m);
  pySpline::python_register_cubic_splines_base_class(m);
  pySpline::python_register_cubic_splines_class(m);
  pySpline::python_register_constant_splines_class(m);
  pySpline::python_register_akima_splines_class(m);
  pySpline::python_register_bessel_splines_class(m);
  pySpline::python_register_hermite_splines_class(m);
  pySpline::python_register_linear_splines_class(m);
  pySpline::python_register_pchip_splines_class(m);
  pySpline::python_register_quintic_splines_base_class(m);
  pySpline::python_register_quintic_splines_class(m);
  pySpline::python_register_hermite_functions(m);
}
