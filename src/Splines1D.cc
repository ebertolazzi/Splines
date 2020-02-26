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

/*!
 |
 | \date     October 28, 2015
 | \version  5.2
 | \note     first release Jan 12, 1998
 |
 | \author   Enrico Bertolazzi
 |
 | \par      Affiliation:
 |           Department of Industrial Engineering<br>
 |           University of Trento<br>
 |           Via Sommarive 9, I -- 38123 Trento, Italy <br>
 |           enrico.bertolazzi@unitn.it
 |
\*/

#include "Splines.hh"

//! Various kind of splines
namespace Splines {

  static
  Spline *
  new_Spline1D( string const & _name, SplineType tp ) {
    switch ( tp ) {
    case CONSTANT_TYPE: return new ConstantSpline(_name);
    case LINEAR_TYPE:   return new LinearSpline(_name);
    case CUBIC_TYPE:    return new CubicSpline(_name);
    case AKIMA_TYPE:    return new AkimaSpline(_name);
    case BESSEL_TYPE:   return new BesselSpline(_name);
    case PCHIP_TYPE:    return new PchipSpline(_name);
    case QUINTIC_TYPE:  return new QuinticSpline(_name);
    default:            break;
    }
    return nullptr;
  }

  void
  Spline1D::build(
    SplineType tp,
    real_type const x[], integer incx,
    real_type const y[], integer incy,
    integer n
  ) {
    pSpline = new_Spline1D(_name,tp);
    SPLINE_ASSERT( pSpline != nullptr, "Spline1D::build, failed" );
    pSpline->build( x, incx, y, incy, n );
  }

  /*
  //    ____  ____   ____                               _
  //   / ___|/ ___| / ___| _   _ _ __  _ __   ___  _ __| |_
  //  | |  _| |     \___ \| | | | '_ \| '_ \ / _ \| '__| __|
  //  | |_| | |___   ___) | |_| | |_) | |_) | (_) | |  | |_
  //   \____|\____| |____/ \__,_| .__/| .__/ \___/|_|   \__|
  //                            |_|   |_|
  */

  #ifndef SPLINES_DO_NOT_USE_GENERIC_CONTAINER

  using GenericContainerNamespace::GC_VEC_REAL;
  using GenericContainerNamespace::vec_real_type;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Spline1D::setup( GenericContainer const & gc ) {
    /*
    // gc["xdata"]
    // gc["ydata"]
    //
    */
    std::string spl_type = gc("spline_type").get_string(
      "Spline1D::setup, spline_type expected to be a string\n"
    );
    SplineType tp;
    if ( spl_type == "constant" ) {
      tp = CONSTANT_TYPE;
    } else if ( spl_type == "linear" ) {
      tp = LINEAR_TYPE;
    } else if ( spl_type == "cubic" ) {
      tp = CUBIC_TYPE;
    } else if ( spl_type == "akima" ) {
      tp = AKIMA_TYPE;
    } else if ( spl_type == "bessel" ) {
      tp = BESSEL_TYPE;
    } else if ( spl_type == "pchip" ) {
      tp = PCHIP_TYPE;
    } else if ( spl_type == "quintic" ) {
      tp = QUINTIC_TYPE;
    } else {
      SPLINE_DO_ERROR(
       "Spline1D::setup[" << _name << "] unknown type " << spl_type <<
       ", not in [constant,linear,cubic,akima,bessel,pchip,quintic]"
      );
    }
    this->pSpline = new_Spline1D( _name, tp );
    this->pSpline->build( gc );
  }
  #endif

}
