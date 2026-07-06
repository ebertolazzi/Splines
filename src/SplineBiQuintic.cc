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

#include "Splines/Splines.hh"

namespace Splines
{

  /*\
   |   ____  _  ___        _       _   _      ____        _ _
   |  | __ )(_)/ _ \ _   _(_)_ __ | |_(_) ___/ ___| _ __ | (_)_ __   ___
   |  |  _ \| | | | | | | | | '_ \| __| |/ __\___ \| '_ \| | | '_ \ / _ \
   |  | |_) | | |_| | |_| | | | | | |_| | (__ ___) | |_) | | | | | |  __/
   |  |____/|_|\__\_\\__,_|_|_| |_|\__|_|\___|____/| .__/|_|_|_| |_|\___|
   |                                               |_|
  \*/
  void BiQuinticSpline::make_spline()
  {
    mDX.resize( m_nx, m_ny );
    mDY.resize( m_nx, m_ny );
    mDXY.resize( m_nx, m_ny );
    mDXX.resize( m_nx, m_ny );
    mDYY.resize( m_nx, m_ny );
    mDXYY.resize( m_nx, m_ny );
    mDXXY.resize( m_nx, m_ny );
    mDXXYY.resize( m_nx, m_ny );

    CubicSpline   cs;
    AkimaSpline   ak;
    VanLeerSpline vl;
    PchipSpline   pc;

    CubicSplineBase * spline{ &pc };
    switch ( m_sub_type )
    {
      case Spline_sub_type::CUBIC: spline = &cs; break;
      case Spline_sub_type::AKIMA: spline = &ak; break;
      case Spline_sub_type::VANLEER: spline = &vl; break;
      case Spline_sub_type::PCHIP: spline = &pc; break;
      default: UTILS_ERROR( "Unknown Spline_sub_type value\n" );
    }

    make_derivative_x( spline, mZ, mDX );
    make_derivative_y( spline, mZ, mDY );
    make_derivative_xy( spline, mDX, mDY, mDXY );

    make_derivative_x( spline, mDX, mDXX );
    make_derivative_y( spline, mDY, mDYY );

    make_derivative_y( spline, mDXX, mDXXY );
    make_derivative_x( spline, mDYY, mDXYY );

    make_derivative_xy( spline, mDXXY, mDXYY, mDXXYY );

    m_search_x.must_reset();
    m_search_y.must_reset();
  }

  void BiQuinticSpline::write_to_stream( ostream_type & s ) const
  {
    fmt::print( s, "Nx = {} Ny = {}\n", m_nx, m_ny );
    for ( integer i = 1; i < m_nx; ++i )
    {
      integer   i1 = i - 1;
      real_type dx = mX.coeff( i ) - mX.coeff( i1 );
      for ( integer j = 1; j < m_ny; ++j )
      {
        integer         j1 = j - 1;
        real_type const dy = mY.coeff( j ) - mY.coeff( j1 );
        fmt::print(
          s,
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
          i,
          j,
          dx,
          dy,
          mZ.coeff( i1, j1 ),
          mZ.coeff( i, j1 ),
          mZ.coeff( i1, j ),
          mZ.coeff( i, j ),
          mDX.coeff( i1, j1 ),
          mDX.coeff( i, j1 ),
          mDX.coeff( i1, j ),
          mDX.coeff( i, j ),
          mDY.coeff( i1, j1 ),
          mDY.coeff( i, j1 ),
          mDY.coeff( i1, j ),
          mDY.coeff( i, j ),
          mDXY.coeff( i1, j1 ),
          mDXY.coeff( i, j1 ),
          mDXY.coeff( i1, j ),
          mDXY.coeff( i, j ) );
      }
    }
  }

}  // namespace Splines

//
// EOF: SplineBiQuintic.cc
//
