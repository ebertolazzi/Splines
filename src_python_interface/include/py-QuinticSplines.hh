/**
 * PYTHON Wrapper for Splines
 *
 * License MIT - See LICENSE file
 * 2019 Matteo Ragni, Claudio Kerov Ghiglianovich, Enrico Bertolazzi
 */

#ifndef PY_QUINTIC_SPLINES_HH
#define PY_QUINTIC_SPLINES_HH

#include <Splines.hh>
#include <pybind11/pybind11.h>

namespace pySpline {
  using Splines::real_type;
  using Splines::integer;
  using Splines::ostream_type;

  using Splines::QuinticSplineBase;

  using pybind11::module;

  class PythonicQuinticSplineBase : public QuinticSplineBase {
    public:

    PythonicQuinticSplineBase string_view name = "QuinticSplineBase" ) : QuinticSplineBase(name) {}

    real_type operator()(real_type x) const override {
      PYBIND11_OVERLOAD_PURE(real_type, QuinticSplineBase, operator(), x);
    }

    real_type D(real_type x) const override {
      PYBIND11_OVERLOAD_PURE(real_type, QuinticSplineBase, D, x);
    }

    real_type DD(real_type x) const override {
      PYBIND11_OVERLOAD_PURE(real_type, QuinticSplineBase, DD, x);
    }

    real_type DDD(real_type x) const override {
      PYBIND11_OVERLOAD_PURE(real_type, QuinticSplineBase, DDD, x);
    }

    real_type DDDD(real_type x) const override {
      PYBIND11_OVERLOAD_PURE(real_type, QuinticSplineBase, DDDD, x);
    }

    real_type DDDDD(real_type x) const override {
      PYBIND11_OVERLOAD_PURE(real_type, QuinticSplineBase, DDDDD, x);
    }

    integer coeffs(real_type cfs[], real_type nodes[], bool transpose = false) const override {
      PYBIND11_OVERLOAD_PURE(integer, QuinticSplineBase, coeffs, cfs, nodes, transpose);
    }

    integer order() const override {
      PYBIND11_OVERLOAD_PURE(integer, QuinticSplineBase, order);
    }

    void write_to_stream(ostream_type & s) const override {
      PYBIND11_OVERLOAD_PURE(void, QuinticSplineBase, write_to_stream, s);
    }

    unsigned type() const override {
      PYBIND11_OVERLOAD_PURE(unsigned, QuinticSplineBase, type);
    }

    void build(real_type const x[], integer incx, real_type const y[], integer incy, integer n) override {
      PYBIND11_OVERLOAD_PURE(void, QuinticSplineBase, build, x, incx, y, incy, n);
    }

    void reserve(integer npts) override {
      PYBIND11_OVERLOAD_PURE(void, QuinticSplineBase, reserve, npts);
    }

    void build() override {
      PYBIND11_OVERLOAD_PURE(void, QuinticSplineBase, build);
    }
  };

  void python_register_quintic_splines_base_class(module & m);
  void python_register_quintic_splines_class(module & m);

}

#endif /* PY_CUBIC_SPLINES_HH */
