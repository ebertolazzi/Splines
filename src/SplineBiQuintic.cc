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
using namespace std; // load standard namspace
#endif

namespace Splines {
#ifdef AUTODIFF_SUPPORT
  //!
  //! \name Autodiff
  //!
  autodiff::dual1st
  BiQuinticSplineBase::eval( autodiff::dual1st const & x, autodiff::dual1st const & y ) const {
    using autodiff::dual1st;
    using autodiff::detail::val;

    real_type dd[3];
    D( val(x), val(y), dd );

    dual1st res{ dd[0] };
    res.grad = dd[1] * x.grad + dd[2] * y.grad;

    return res;
  }

  autodiff::dual2nd
  BiQuinticSplineBase::eval( autodiff::dual2nd const & x, autodiff::dual2nd const & y ) const {
    using autodiff::dual2nd;
    using autodiff::derivative;

    real_type dd[6], dx{ val(x.grad) }, dy{ val(y.grad) }, ddx{ x.grad.grad }, ddy{ y.grad.grad };
    DD( val(x), val(y), dd );

    dual2nd res{ dd[0] };
    res.grad = dd[1] * dx + dd[2] * dy;
    res.grad.grad = dd[3]*dx*dx + 2*dx*dy*dd[4]+ dy*dy*dd[5] + ddx*dd[1] + ddy*dd[2];
    return res;
  }
  #endif
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiQuinticSpline::make_spline() {

    integer const dim{ m_nx*m_ny };
    m_mem_biquintic.reallocate( 8*dim );
    m_DX    = m_mem_biquintic( dim );
    m_DY    = m_mem_biquintic( dim );
    m_DXY   = m_mem_biquintic( dim );
    m_DXX   = m_mem_biquintic( dim );
    m_DYY   = m_mem_biquintic( dim );
    m_DXYY  = m_mem_biquintic( dim );
    m_DXXY  = m_mem_biquintic( dim );
    m_DXXYY = m_mem_biquintic( dim );

    make_derivative_x( m_Z, m_DX );
    make_derivative_y( m_Z, m_DY );
    make_derivative_xy( m_DX, m_DY, m_DXY );

    make_derivative_x( m_DX, m_DXX );
    make_derivative_y( m_DY, m_DYY );

    make_derivative_y( m_DXX, m_DXXY );
    make_derivative_x( m_DYY, m_DXYY );

    make_derivative_xy( m_DXXY, m_DXYY, m_DXXYY );

    m_search_x.must_reset();
    m_search_y.must_reset();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiQuinticSpline::write_to_stream( ostream_type & s ) const {
    fmt::print( s, "Nx = {} Ny = {}\n", m_nx, m_ny );
    for ( integer i{1}; i < m_nx; ++i ) {
      real_type dx{ m_X[i]-m_X[i-1] };
      for ( integer j{1}; j < m_ny; ++j ) {
        integer   const i00 { ipos_C(i-1,j-1) };
        integer   const i10 { ipos_C(i,j-1) };
        integer   const i01 { ipos_C(i-1,j) };
        integer   const i11 { ipos_C(i,j) };
        real_type const dy  { m_Y[j]-m_Y[j-1] };
        fmt::print( s,
          "patch ({},{})\n"
          "  DX    = {:<12.4}  DY    = {:<12.4}\n"
          "  Z00   = {:<12.4}  Z10   = {:<12.4}\n"
          "  Z01   = {:<12.4}  Z11   = {:<12.4}\n"
          "  Dx00  = {:<12.4}  Dx10  = {:<12.4}\n"
          "  Dx01  = {:<12.4}  Dx11  = {:<12.4}\n"
          "  Dy00  = {:<12.4}  Dy10  = {:<12.4}\n"
          "  Dy01  = {:<12.4}  Dy11  = {:<12.4}\n"
          "  Dxy00 = {:<12.4}  Dxy10 = {:<12.4}\n"
          "  Dxy01 = {:<12.4}  Dxy11 = {:<12.4}\n",
          i, j, dx, dy,
          m_Z[i00],   m_Z[i10],   m_Z[i01],   m_Z[i11],
          m_DX[i00],  m_DX[i10],  m_DX[i01],  m_DX[i11],
          m_DY[i00],  m_DY[i10],  m_DY[i01],  m_DY[i11],
          m_DXY[i00], m_DXY[i10], m_DXY[i01], m_DXY[i11]
        );
      }
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  char const *
  BiQuinticSpline::type_name() const
  { return "BiQuintic"; }

}
