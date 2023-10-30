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
 |      Universita` degli Studi di Trento                                   |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#include "Splines.hh"

#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
#pragma clang diagnostic ignored "-Wpoison-system-directories"
#pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
#endif

namespace Splines {

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Spline2D::new_spline( SplineType2D tp ) {
    if ( m_spline_2D == nullptr ) {
      delete m_spline_2D;
      m_spline_2D = nullptr;
    }
    switch ( tp ) {
    case SplineType2D::BILINEAR:  m_spline_2D = new BilinearSpline(m_name);  break;
    case SplineType2D::BICUBIC:   m_spline_2D = new BiCubicSpline(m_name);   break;
    case SplineType2D::BIQUINTIC: m_spline_2D = new BiQuinticSpline(m_name); break;
    case SplineType2D::AKIMA2D:   m_spline_2D = new Akima2Dspline(m_name);   break;
//    default:
//      UTILS_ERROR( "new_spline, type `{}` unknown\n", tp );
    }
  }

}
