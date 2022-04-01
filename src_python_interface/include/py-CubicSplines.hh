/**
 * PYTHON Wrapper for Splines
 *
 * License MIT - See LICENSE file
 * 2019 Matteo Ragni, Claudio Kerov Ghiglianovich, Enrico Bertolazzi
 */

#ifndef PY_CUBIC_SPLINES_HH
#define PY_CUBIC_SPLINES_HH

#include <Splines.hh>
#include <pybind11/pybind11.h>

namespace pySpline {
  using Splines::real_type;
  using Splines::integer;
  using Splines::ostream_type;

  using Splines::CubicSplineBase;

  using pybind11::module;

  class PythonicCubicSplineBase : public CubicSplineBase {
    public:

    PythonicCubicSplineBase(std::string const & name = "CubicSplineBase") : CubicSplineBase(name) {}

    real_type operator()(real_type x) const override {
      PYBIND11_OVERLOAD_PURE(real_type, CubicSplineBase, operator(), x);
    }

    real_type D(real_type x) const override {
      PYBIND11_OVERLOAD_PURE(real_type, CubicSplineBase, D, x);
    }

    real_type DD(real_type x) const override {
      PYBIND11_OVERLOAD_PURE(real_type, CubicSplineBase, DD, x);
    }

    real_type DDD(real_type x) const override {
      PYBIND11_OVERLOAD_PURE(real_type, CubicSplineBase, DDD, x);
    }

    integer coeffs(real_type cfs[], real_type nodes[], bool transpose = false) const override {
      PYBIND11_OVERLOAD_PURE(integer, CubicSplineBase, coeffs, cfs, nodes, transpose);
    }

    integer order() const override {
      PYBIND11_OVERLOAD_PURE(integer, CubicSplineBase, order);
    }

    void write_to_stream(ostream_type & s) const override {
      PYBIND11_OVERLOAD_PURE(void, CubicSplineBase, write_to_stream, s);
    }

    unsigned type() const override {
      PYBIND11_OVERLOAD_PURE(unsigned, CubicSplineBase, type);
    }

    void build(real_type const x[], integer incx, real_type const y[], integer incy, integer n) override {
      PYBIND11_OVERLOAD_PURE(void, CubicSplineBase, build, x, incx, y, incy, n);
    }

    void reserve(integer npts) override {
      PYBIND11_OVERLOAD_PURE(void, CubicSplineBase, reserve, npts);
    }

    void build() override {
      PYBIND11_OVERLOAD_PURE(void, CubicSplineBase, build);
    }
  };

  void python_register_cubic_splines_base_class(module & m);
  void python_register_cubic_splines_class(module & m);

}

#endif /* PY_CUBIC_SPLINES_HH */
