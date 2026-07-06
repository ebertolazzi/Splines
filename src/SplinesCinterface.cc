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

//
// file: SplineCinterface.cc
//

//!
//! \file SplinesCinterface.cc
//! This file contains the sources for the C interface to `Splines`
//!
#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wpoison-system-directories"
#pragma clang diagnostic ignored "-Wundefined-func-template"
#pragma clang diagnostic ignored "-Wexit-time-destructors"
#endif

#include "Splines/Splines.hh"
#include "Splines/SplinesCinterface.h"
#include "Utils_fmt.hh"

using namespace SplinesLoad;

#include <map>
#include <mutex>
#include <string>

#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#endif

namespace
{
  // This C interface keeps a single process-wide registry (`spline_stored`)
  // plus a "currently selected" `head`, all mutated by the entry points
  // below. Every entry point funnels through one of the c_api_call_* wrappers,
  // so taking one lock there makes each C-API call atomic -- no data races on
  // the registry/head and no use-after-free between a delete on one thread and
  // an eval on another. The API remains logically stateful (two threads
  // sharing `head` is a logic hazard, not a memory-safety one); heavy
  // concurrent evaluation should use the C++ API directly. A function-local
  // static sidesteps any static-init-order concerns.
  std::mutex & c_api_mutex() noexcept
  {
    static std::mutex m;
    return m;
  }

  template <typename F> int
  c_api_call_int( F const & fn ) noexcept
  {
    try
    {
      std::lock_guard<std::mutex> const lock( c_api_mutex() );
      return fn();
    }
    catch ( ... )
    {
      return -1;
    }
  }

  template <typename F> double
  c_api_call_real( F const & fn ) noexcept
  {
    try
    {
      std::lock_guard<std::mutex> const lock( c_api_mutex() );
      return fn();
    }
    catch ( ... )
    {
      return 0;
    }
  }

  template <typename F> char const *
  c_api_call_cstr( F const & fn ) noexcept
  {
    try
    {
      std::lock_guard<std::mutex> const lock( c_api_mutex() );
      return fn();
    }
    catch ( ... )
    {
      return "ERROR";
    }
  }

  template <typename F> void *
  c_api_call_ptr( F const & fn ) noexcept
  {
    try
    {
      std::lock_guard<std::mutex> const lock( c_api_mutex() );
      return fn();
    }
    catch ( ... )
    {
      return nullptr;
    }
  }

  Spline *
  make_spline( SplineType1D const type )
  {
    switch ( type )
    {
      case SplineType1D::AKIMA: return new AkimaSpline;
      case SplineType1D::VANLEER: return new VanLeerSpline;
      case SplineType1D::PCHIP: return new PchipSpline;
      case SplineType1D::CUBIC: return new CubicSpline;
      case SplineType1D::LINEAR: return new LinearSpline;
      case SplineType1D::CONSTANT: return new ConstantSpline;
      case SplineType1D::QUINTIC_CUBIC: return new QuinticSpline( Splines::Spline_sub_type::CUBIC );
      case SplineType1D::QUINTIC_AKIMA: return new QuinticSpline( Splines::Spline_sub_type::AKIMA );
      case SplineType1D::QUINTIC_VANLEER: return new QuinticSpline( Splines::Spline_sub_type::VANLEER );
      case SplineType1D::QUINTIC_PCHIP: return new QuinticSpline( Splines::Spline_sub_type::PCHIP );
      default: return nullptr;
    }
  }
}  // namespace

extern "C"
{
  typedef std::map<std::string, Spline *> MAP_SPLINE;

  static std::map<std::string, Spline *> spline_stored;
  static Spline *                        head = nullptr;

  int SPLINE_new( char const id[], char const type[] )
  {
    return c_api_call_int(
      [&]() -> int
      {
        if ( id == nullptr || type == nullptr ) return -1;

        if ( auto const it = spline_stored.find( id ); it != spline_stored.end() )
        {
          if ( head == it->second ) head = nullptr;
          delete it->second;
          spline_stored.erase( it );
        }

        Spline * spline = make_spline( Splines::string_to_splineType1D( type ) );
        if ( spline == nullptr ) return -1;

        head             = spline;
        spline_stored[id] = spline;
        return 0;
      } );
  }

  int SPLINE_select( char const id[] )
  {
    return c_api_call_int(
      [&]() -> int
      {
        if ( id == nullptr ) return -1;
        if ( auto const it = spline_stored.find( id ); it != spline_stored.end() )
        {
          head = it->second;
          return 0;
        }
        return -1;
      } );
  }

  int SPLINE_delete( char const id[] )
  {
    return c_api_call_int(
      [&]() -> int
      {
        if ( id == nullptr ) return -1;
        if ( auto const it = spline_stored.find( id ); it != spline_stored.end() )
        {
          if ( head == it->second ) head = nullptr;
          delete it->second;
          spline_stored.erase( it );
          return 0;
        }
        return -1;
      } );
  }

  int SPLINE_print()
  {
    return c_api_call_int(
      [&]() -> int
      {
        if ( head != nullptr )
        {
          head->write_to_stream( std::cout );
          return 0;
        }
        std::cout << "No Spline!\n";
        return -1;
      } );
  }

  char const * SPLINE_get_type_name()
  {
    return c_api_call_cstr(
      [&]() -> char const *
      {
        if ( head == nullptr ) return "NOTYPE - head = nullptr";
        return head->type_name();
      } );
  }

  void * SPLINE_mem_ptr( char const id[] )
  {
    return c_api_call_ptr(
      [&]() -> void *
      {
        if ( id == nullptr ) return nullptr;
        if ( auto const it = spline_stored.find( id ); it != spline_stored.end() )
          return static_cast<void *>( it->second );
        return nullptr;
      } );
  }

  int SPLINE_init()
  {
    return c_api_call_int(
      [&]() -> int
      {
        if ( head != nullptr )
        {
          head->clear();
          return 0;
        }
        return -1;
      } );
  }

  int SPLINE_push( double const x, double const y )
  {
    return c_api_call_int(
      [&]() -> int
      {
        if ( head != nullptr )
        {
          head->push_back( x, y );
          return 0;
        }
        return -1;
      } );
  }

  int SPLINE_build()
  {
    return c_api_call_int(
      [&]() -> int
      {
        if ( head != nullptr )
        {
          head->build();
          return 0;
        }
        return -1;
      } );
  }

  int SPLINE_build2( double const x[], double const y[], int const n )
  {
    return c_api_call_int(
      [&]() -> int
      {
        if ( head != nullptr )
        {
          head->build( x, y, n );
          return 0;
        }
        return -1;
      } );
  }

  double SPLINE_eval( double const x )
  {
    return c_api_call_real(
      [&]() -> double
      {
        if ( head != nullptr ) return head->eval( x );
        return 0;
      } );
  }

  double SPLINE_eval_D( double const x )
  {
    return c_api_call_real(
      [&]() -> double
      {
        if ( head != nullptr ) return head->D( x );
        return 0;
      } );
  }

  double SPLINE_eval_DD( double const x )
  {
    return c_api_call_real(
      [&]() -> double
      {
        if ( head != nullptr ) return head->DD( x );
        return 0;
      } );
  }

  double SPLINE_eval_DDD( double const x )
  {
    return c_api_call_real(
      [&]() -> double
      {
        if ( head != nullptr ) return head->DDD( x );
        return 0;
      } );
  }

  double SPLINE_eval_DDDD( double const x )
  {
    return c_api_call_real(
      [&]() -> double
      {
        if ( head != nullptr ) return head->DDDD( x );
        return 0;
      } );
  }

  double SPLINE_eval_DDDDD( double const x )
  {
    return c_api_call_real(
      [&]() -> double
      {
        if ( head != nullptr ) return head->DDDDD( x );
        return 0;
      } );
  }
}

//
// eof: SplineCinterface.cc
//
