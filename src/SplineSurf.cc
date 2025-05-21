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

  void
  SplineSurf::make_derivative_x( real_type const z[], real_type dx[] ) {
    PchipSpline pchip_work;
    for ( integer j{0}; j < m_ny; ++j ) {
      pchip_work.build( m_X, 1, z + ipos_C(0,j), m_ny, m_nx );
      for ( integer i{0}; i < m_nx; ++i ) dx[ipos_C(i,j)] = pchip_work.yp_node(i);
    }
  }

  void
  SplineSurf::make_derivative_y( real_type const z[], real_type dy[] ) {
    PchipSpline pchip_work;
    for ( integer i{0}; i < m_nx; ++i ) {
      pchip_work.build( m_Y, 1, z + ipos_C(i,0), 1, m_ny );
      for ( integer j{0}; j < m_ny; ++j ) dy[ipos_C(i,j)] = pchip_work.yp_node(j);
    }
  }

  void
  SplineSurf::make_derivative_xy( real_type const dx[], real_type const dy[], real_type dxy[] ) {
    PchipSpline pchip_work;

    auto minmod = [] ( real_type a, real_type b ) -> real_type {
      if ( a*b <= 0 ) return 0;
      if ( a > 0    ) return std::min(a,b);
      return std::max(a,b);
    };

    for ( integer j{0}; j < m_ny; ++j ) {
      pchip_work.build( m_X, 1, dy + ipos_C(0,j), m_ny, m_nx );
      for ( integer i{0}; i < m_nx; ++i ) dxy[ipos_C(i,j)] = pchip_work.yp_node(i);
    }

    for ( integer i{0}; i < m_nx; ++i ) {
      pchip_work.build( m_Y, 1, dx + ipos_C(i,0), 1, m_ny );
      for ( integer j{0}; j < m_ny; ++j ) {
        integer const ij{ipos_C(i,j)};
        dxy[ij] = minmod( dxy[ij], pchip_work.yp_node(j) );
      }
    }
  }

}
