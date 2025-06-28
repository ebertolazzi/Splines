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

namespace Splines {

  #ifndef DOXYGEN_SHOULD_SKIP_THIS

  void
  Spline2D::new_spline( SplineType2D const tp ) {
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

  #endif

  void
  Spline2D::setup( GenericContainer const & gc ) {
    string const   where{ fmt::format("Spline2D[{}]::setup( gc ):", m_name ) };
    string const & type{ gc.get_map_string("spline_type",where) };
    new_spline( string_to_splineType2D( type ) );
    m_spline_2D->setup( gc );
  }

}
