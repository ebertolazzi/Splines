/**
 * PYTHON Wrapper for Splines
 * 
 * License MIT - See LICENSE file
 * 2019 Matteo Ragni, Claudio Kerov Ghiglianovich, Enrico Bertolazzi
 */

#ifndef PY_SPLINES_HH
#define PY_SPLINES_HH

#include <string>

#include <Splines.hh>
#include <pybind11/pybind11.h>

namespace pySpline {
  using Splines::real_type;
  using Splines::integer;
  using Splines::ostream_type;

  using Splines::Spline;
  using Splines::real_type;
  using Splines::integer;

  using pybind11::module;
  
  using GenericContainerNamespace::GenericContainer;

  class PythonicSpline : public Spline {
    public:
    
    PythonicSpline(std::string const & name = "Spline") : Spline(name) {}
    
    void reserve(integer npts) override {
      PYBIND11_OVERLOAD_PURE(void, Spline, reserve, npts);
    }

    void build() override {
      PYBIND11_OVERLOAD_PURE(void, Spline, build);
    }    

    void setup(GenericContainer const & gc) override {
      PYBIND11_OVERLOAD_PURE(void, Spline, setup, gc);
    }

    void build(real_type const x[], integer incx, real_type const y[], integer incy, integer n) override {
      PYBIND11_OVERLOAD_PURE(void, Spline, build, x, incx, y, incy, n);
    }

    void clear() override {
      PYBIND11_OVERLOAD_PURE(void, Spline, clear);
    }

    unsigned type() const override {
      PYBIND11_OVERLOAD_PURE(unsigned, Spline, type);
    }

    real_type operator()(real_type x) const override {
      PYBIND11_OVERLOAD_PURE(real_type, Spline, operator(), x);
    }

    real_type D(real_type x) const override {
      PYBIND11_OVERLOAD_PURE(real_type, Spline, D, x);
    }

    real_type DD(real_type x) const override {
      PYBIND11_OVERLOAD_PURE(real_type, Spline, DD, x);
    }

    real_type DDD(real_type x) const override {
      PYBIND11_OVERLOAD_PURE(real_type, Spline, DDD, x);
    }

    real_type DDDD(real_type x) const override {
      PYBIND11_OVERLOAD_PURE(real_type, Spline, DDDD, x);
    }

    real_type DDDDD(real_type x) const override {
      PYBIND11_OVERLOAD_PURE(real_type, Spline, DDDDD, x);
    }

    integer coeffs(real_type cfs[], real_type nodes[], bool transpose = false) const override {
      PYBIND11_OVERLOAD_PURE(integer, Spline, coeffs, cfs, nodes, transpose);
    }

    integer order() const override {
      PYBIND11_OVERLOAD_PURE(integer, Spline, order);
    }

    void writeToStream(ostream_type & s) const override {
      PYBIND11_OVERLOAD_PURE(void, Spline, writeToStream, s);
    }
  };

  void python_register_splines_class(module & m);
  void python_register_hermite_functions(module & m);
}

#endif /* PY_SPLINES_HH */
