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

#define AUTODIFF_SUPPORT
#include "Splines/Splines.hh"

namespace Splines
{

  /*\
   |   ____  _  ____      _     _      ____        _ _
   |  | __ )(_)/ ___|   _| |__ (_) ___/ ___| _ __ | (_)_ __   ___
   |  |  _ \| | |  | | | | '_ \| |/ __\___ \| '_ \| | | '_ \ / _ \
   |  | |_) | | |__| |_| | |_) | | (__ ___) | |_) | | | | | |  __/
   |  |____/|_|\____\__,_|_.__/|_|\___|____/| .__/|_|_|_| |_|\___|
   |                                        |_|
  \*/


  void BiCubicSpline::make_spline()
  {
    mDX.resize( m_nx, m_ny );
    mDY.resize( m_nx, m_ny );
    mDXY.resize( m_nx, m_ny );

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

    m_search_x.must_reset();
    m_search_y.must_reset();
  }


  void BiCubicSpline::write_to_stream( ostream_type & s ) const
  {
    // Stampa intestazione
    fmt::print( s, "Nx = {} Ny = {}\n", m_nx, m_ny );

    for ( integer i = 1; i < m_nx; ++i )
    {
      // Accesso diretto al vettore mappato
      real_type const dx = mX.coeff( i ) - mX.coeff( i - 1 );

      for ( integer j = 1; j < m_ny; ++j )
      {
        real_type const dy = mY.coeff( j ) - mY.coeff( j - 1 );

        // Indici base per la patch corrente
        integer const r0 = i - 1;
        integer const r1 = i;
        integer const c0 = j - 1;
        integer const c1 = j;

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
          // Z values
          mZ( r0, c0 ),
          mZ( r1, c0 ),
          mZ( r0, c1 ),
          mZ( r1, c1 ),
          // DX values
          mDX( r0, c0 ),
          mDX( r1, c0 ),
          mDX( r0, c1 ),
          mDX( r1, c1 ),
          // DY values
          mDY( r0, c0 ),
          mDY( r1, c0 ),
          mDY( r0, c1 ),
          mDY( r1, c1 ),
          // DXY values
          mDXY( r0, c0 ),
          mDXY( r1, c0 ),
          mDXY( r0, c1 ),
          mDXY( r1, c1 ) );
      }
    }
  }

}  // namespace Splines

//
// EOF: SplineBiCubic.cc
//
