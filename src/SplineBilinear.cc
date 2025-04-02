/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2016                                                      |
 |                                                                          |
 |         , __                 , __                                        |
 |        /|/  \               /|/  \                                       |
 |         | __/ _   ,_         | __/ _   ,_                                |
 |         |   \|/  /  |  |   | |   \|/  /  |  |   |                        |
 |         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                       |
 |                           /|                   /|                        |
 |                           \|                   \|                        |
 |                                                                          |
 |      Enrico Bertolazzi                                                   |
 |      Dipartimento di Ingegneria Industriale                              |
 |      Università degli Studi di Trento                                    |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wpoison-system-directories"
#pragma clang diagnostic ignored "-Wundefined-func-template"
#endif

#include "Splines.hh"
#include "Utils_fmt.hh"

#ifndef DOXYGEN_SHOULD_SKIP_THIS
using namespace std; // load standard namespace
#endif

namespace Splines {

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BilinearSpline::eval( real_type x, real_type y ) const {
    integer   const i   { this->search_x( x ) };
    integer   const j   { this->search_y( y ) };
    real_type const DX  { m_X[i+1] - m_X[i] };
    real_type const DY  { m_Y[j+1] - m_Y[j] };
    real_type const u   { (x-m_X[i])/DX };
    real_type const v   { (y-m_Y[j])/DY };
    real_type const u1  { 1-u };
    real_type const v1  { 1-v };
    real_type const Z00 { m_Z[this->ipos_C(i,j)] };
    real_type const Z01 { m_Z[this->ipos_C(i,j+1)] };
    real_type const Z10 { m_Z[this->ipos_C(i+1,j)] };
    real_type const Z11 { m_Z[this->ipos_C(i+1,j+1)] };
    return u1 * ( Z00 * v1 + Z01 * v ) +
           u  * ( Z10 * v1 + Z11 * v );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BilinearSpline::Dx( real_type x, real_type y ) const {
    integer   const i   { this->search_x( x ) };
    integer   const j   { this->search_y( y ) };
    real_type const DX  { m_X[i+1] - m_X[i] };
    real_type const DY  { m_Y[j+1] - m_Y[j] };
    real_type const v   { (y-m_Y[j])/DY };
    real_type const Z00 { m_Z[ipos_C(i,j)] };
    real_type const Z01 { m_Z[ipos_C(i,j+1)] };
    real_type const Z10 { m_Z[ipos_C(i+1,j)] };
    real_type const Z11 { m_Z[ipos_C(i+1,j+1)] };
    return ( (Z10-Z00) * (1-v) + (Z11-Z01) * v ) / DX;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BilinearSpline::Dy( real_type x, real_type y ) const {
    integer   const i   { this->search_x( x ) };
    integer   const j   { this->search_y( y ) };
    real_type const DX  { m_X[i+1] - m_X[i] };
    real_type const DY  { m_Y[j+1] - m_Y[j] };
    real_type const u   { (x-m_X[i])/DX };
    real_type const Z00 { m_Z[ipos_C(i,j)] };
    real_type const Z01 { m_Z[ipos_C(i,j+1)] };
    real_type const Z10 { m_Z[ipos_C(i+1,j)] };
    real_type const Z11 { m_Z[ipos_C(i+1,j+1)] };
    return ( ( Z01-Z00 ) * (1-u) + ( Z11-Z10 ) * u ) / DY;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BilinearSpline::D( real_type x, real_type y, real_type d[3] ) const {
    integer   const i   { this->search_x( x ) };
    integer   const j   { this->search_y( y ) };
    real_type const DX  { m_X[i+1] - m_X[i] };
    real_type const DY  { m_Y[j+1] - m_Y[j] };
    real_type const u   { (x-m_X[i])/DX };
    real_type const v   { (y-m_Y[j])/DY };
    real_type const u1  { 1-u };
    real_type const v1  { 1-v };
    real_type const Z00 { m_Z[this->ipos_C(i,j)] };
    real_type const Z01 { m_Z[this->ipos_C(i,j+1)] };
    real_type const Z10 { m_Z[this->ipos_C(i+1,j)] };
    real_type const Z11 { m_Z[this->ipos_C(i+1,j+1)] };
    d[0] = u1 * ( Z00 * v1 + Z01 * v ) + u * ( Z10 * v1 + Z11 * v );
    d[1] = v1 * (Z10-Z00) + v * (Z11-Z01); d[1] /= DX;
    d[2] = u1 * (Z01-Z00) + u * (Z11-Z10); d[2] /= DY;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BilinearSpline::DD( real_type const x, real_type const y, real_type dd[6] ) const {
    this->D( x, y, dd );
    dd[3] = dd[4] = dd[5] = 0; // second derivative are 0
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BilinearSpline::write_to_stream( ostream_type & s ) const {
    fmt::print( s, "Nx={} Ny={}\n", m_nx, m_ny );
    for ( integer i{1}; i < m_nx; ++i ) {
      for ( integer j{1}; j < m_ny; ++j ) {
        integer const i00 { this->ipos_C(i-1,j-1) };
        integer const i10 { this->ipos_C(i,j-1) };
        integer const i01 { this->ipos_C(i-1,j) };
        integer const i11 { this->ipos_C(i,j) };
        fmt::print( s,
          "patch ({},{})\n"
          "  DX  = {:<12.4}  DY  = {:<12.4}\n"
          "  Z00 = {:<12.4}  Z01 = {:<12.4}  Z10 = {:<12.4}  Z11 = {:<12.4}\n",
          i, j,
          m_X[i]-m_X[i-1],
          m_Y[j]-m_Y[j-1],
          m_Z[i00], m_Z[i01], m_Z[i10], m_Z[i11]
        );
      }
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  char const *
  BilinearSpline::type_name() const
  { return "bilinear"; }

}
