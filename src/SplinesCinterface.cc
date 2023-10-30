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

//
// file: SplineCinterface.cc
//

#ifndef DOXYGEN_SHOULD_SKIP_THIS

//!
//! \file SplinesCinterface.cc
//! This file contains the sources for the C interface to `Splines`
//!

#include "Splines.hh"
#include "SplinesCinterface.h"

#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
#pragma clang diagnostic ignored "-Wpoison-system-directories"
#endif

using namespace SplinesLoad;

#include <vector>
#include <map>
#include <string>
#include <deque>

#include <string.h>

#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#endif

using namespace std; // load standard namspace

extern "C" {

  typedef std::map<std::string,Spline*> MAP_SPLINE;

  static std::map<std::string,Spline*> spline_stored;
  static Spline * head = nullptr;

  int
  SPLINE_new(
    char const * id,
    char const * type
  ) {
    fmt::print( "SPLINE_new, id = {} type = {}\n", id, type );
    MAP_SPLINE::iterator it = spline_stored.find(id);
    if ( it != spline_stored.end() ) delete it->second;
    int ok = 0;
    switch ( Splines::string_to_splineType1D(type) ) {
    case SplineType1D::AKIMA:    head = new AkimaSpline;    break;
    case SplineType1D::BESSEL:   head = new BesselSpline;   break;
    case SplineType1D::PCHIP:    head = new PchipSpline;    break;
    case SplineType1D::CUBIC:    head = new CubicSpline;    break;
    case SplineType1D::LINEAR:   head = new LinearSpline;   break;
    case SplineType1D::CONSTANT: head = new ConstantSpline; break;
    case SplineType1D::QUINTIC:  head = new QuinticSpline;  break;
    default:
      head = nullptr;
      ok   = -1;
      break;
    }
    spline_stored[id] = head;
    return ok;
  }

  int
  SPLINE_select( char const * id ) {
    MAP_SPLINE::iterator it = spline_stored.find(id);
    if ( it != spline_stored.end() ) {
      head = it->second;
    } else {
      return -1; // spline non trovata
    }
    return 0;
  }

  int
  SPLINE_delete( char const * id ) {
    MAP_SPLINE::iterator it = spline_stored.find(id);
    if ( it != spline_stored.end() ) {
      delete it->second;
      spline_stored.erase(it);
      head = nullptr;
    } else {
      return -1; // spline non trovata
    }
    return 0;
  }

  int
  SPLINE_print() {
    if ( head != nullptr ) {
      head->write_to_stream( cout );
      return 0;
    } else {
      cout << "No Spline!\n";
      return -1;
    }
  }

  char const *
  SPLINE_get_type_name() {
    if ( head == nullptr ) return "NOTYPE - head = nullptr";
    return head->type_name();
  }

  void *
  SPLINE_mem_ptr( char const * id ) {
    // check if exists ?
    return static_cast<void*>(&spline_stored[id]);
  }

  int
  SPLINE_init() {
    if ( head != nullptr ) {
      head->clear();
      return 0;
    } else {
      return -1;
    }
  }

  int
  SPLINE_push( double x, double y ) {
    if ( head != nullptr ) {
      head->push_back(x,y);
      return 0;
    } else {
      return -1;
    }
  }

  int
  SPLINE_build() {
    if ( head != nullptr ) {
      head->build();
      return 0;
    } else {
      return -1;
    }
  }

  int
  SPLINE_build2(
    double const * x,
    double const * y,
    int            n
  ) {
    if ( head != nullptr ) {
      head->build( x, y, n );
      return 0;
    } else {
      return -1;
    }
  }

  double
  SPLINE_eval( double x ) {
    if ( head != nullptr ) {
      return head->operator()(x);
    } else {
      return 0;
    }
  }

  double
  SPLINE_eval_D( double x ) {
    if ( head != nullptr ) {
      return head->D(x);
    } else {
      return 0;
    }
  }

  double
  SPLINE_eval_DD( double x ) {
    if ( head != nullptr ) {
      return head->DD(x);
    } else {
      return 0;
    }
  }

  double
  SPLINE_eval_DDD( double x ) {
    if ( head != nullptr ) {
      return head->DDD(x);
    } else {
      return 0;
    }
  }

  double
  SPLINE_eval_DDDD( double x ) {
    if ( head != nullptr ) {
      return head->DDDD(x);
    } else {
      return 0;
    }
  }

  double
  SPLINE_eval_DDDDD( double x ) {
    if ( head != nullptr ) {
      return head->DDDDD(x);
    } else {
      return 0;
    }
  }

}

#endif

//
// eof: SplineCinterface.cc
//
