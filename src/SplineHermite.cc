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

  #ifndef DOXYGEN_SHOULD_SKIP_THIS

  /*\
   |   _   _                     _ _
   |  | | | | ___ _ __ _ __ ___ (_) |_ ___
   |  | |_| |/ _ \ '__| '_ ` _ \| | __/ _ \
   |  |  _  |  __/ |  | | | | | | | ||  __/
   |  |_| |_|\___|_|  |_| |_| |_|_|\__\___|
  \*/

  void
  Hermite3( real_type const x, real_type const H, real_type base[4] ) {
    real_type const X{x/H};
    base[1] = X*X*(3-2*X);
    base[0] = 1-base[1];
    base[2] = x*(X*(X-2)+1);
    base[3] = x*X*(X-1);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Hermite3_D( real_type const x, real_type const H, real_type base_D[4] ) {
    real_type const X{x/H};
    base_D[0] = 6.0*X*(X-1.0)/H;
    base_D[1] = -base_D[0];
    base_D[2] = ((3*X-4)*X+1);
    base_D[3] = X*(3*X-2);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Hermite3_DD( real_type const x, real_type const H, real_type base_DD[4] ) {
    real_type const X{x/H};
    base_DD[0] = (12*X-6)/(H*H);
    base_DD[1] = -base_DD[0];
    base_DD[2] = (6*X-4)/H;
    base_DD[3] = (6*X-2)/H;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Hermite3_DDD( real_type, real_type const H, real_type base_DDD[4] ) {
    base_DDD[0] = 12/(H*H*H);
    base_DDD[1] = -base_DDD[0];
    base_DDD[2] = 6/(H*H);
    base_DDD[3] = base_DDD[2];
  }

  // --------------------------------------------------------------------------

  void
  Hermite5( real_type const x, real_type const H, real_type base[6] ) {
    real_type const t1  = H*H;
    real_type const t4  = x*x;
    real_type const t7  = H-x;
    real_type const t8  = t7*t7;
    real_type const t9  = t8*t7;
    real_type const t11 = t1*t1;
    real_type const t2  = 1/t11;
    real_type const t3  = 1/H;
    real_type const t13 = t3*t2;
    real_type const t14 = t4*x;
    real_type const t17 = t4*t4;
    base[0] = t13*t9*(3.0*x*H+t1+6.0*t4);
    base[1] = t13*(-15.0*H*t17+6.0*t17*x+10.0*t1*t14);
    base[2] = t2*t9*x*(H+3*x);
    base[3] = t2*(3*x-4*H)*t7*t14;
    real_type const t36 = t3/t1/2;
    base[4] = t36*t9*t4;
    base[5] = t36*t8*t14;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Hermite5_D( real_type const x, real_type const H, real_type base_D[6] ) {
    real_type const t1 = H-x;
    real_type const t2 = t1*t1;
    real_type const t3 = x*x;
    real_type const t5 = H*H;
    real_type const t6 = t5*t5;
    real_type const t4 = 1/t6;
    real_type const t7 = 1/H;
    real_type const t10 = 30.0*t3*t2*t7*t4;
    real_type const t11 = 5.0*x;
    real_type const t23 = t3*t3;
    real_type const t30 = t7/t5/2;
    base_D[0] = -t10;
    base_D[1] = t10;
    base_D[2] = t4*(H-3.0*x)*(H+t11)*t2;
    base_D[3] = t4*(t3*(28*x*H-12*t5)-15*t23);
    base_D[4] = t30*(2*H-t11)*t2*x;
    base_D[5] = t30*(3*H-t11)*t3*t1;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Hermite5_DD( real_type const x, real_type const H, real_type base_DD[6] ) {
    real_type const t1 = H-x;
    real_type const t2 = x*t1;
    real_type const t5 = H*H;
    real_type const t6 = t5*t5;
    real_type const t3 = 1/t6;
    real_type const t4 = 1/H;
    real_type const t11 = 60*(H-2*x)*t2*t4*t3;
    real_type const t26 = x*x;
    real_type const t31 = t4/t5;
    base_DD[0] = -t11;
    base_DD[1] = t11;
    base_DD[2] = 12*t3*t1*(5*x-3*H)*x;
    base_DD[3] = 12*t3*t2*(5*x-2*H);
    base_DD[4] = t31*(10*t26+t5-8*H*x)*t1;
    base_DD[5] = t31*(t26*(10*x-12*H)+3*x*t5);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Hermite5_DDD( real_type x, real_type H, real_type base_DDD[6] ) {
    real_type const t1  = H*H;
    real_type const t3  = H*x;
    real_type const t5  = x*x;
    real_type const t7  = 360*t3-60*t1-360*t5;
    real_type const t8  = t1*t1;
    real_type const t9  = 1/t8;
    real_type const t11 = 1/H;
    real_type const t10 = t11*t9;
    real_type const t14 = 180.0*t5;
    real_type const t22 = 30.0*t5;
    real_type const t25 = t11/t1;
    base_DDD[0] = t7*t10;
    base_DDD[1] = -base_DDD[0];
    base_DDD[2] = t9*(192*t3-36*t1-t14);
    base_DDD[3] = t9*(168*t3-24*t1-t14);
    base_DDD[4] = t25*(36*t3-9*t1-t22);
    base_DDD[5] = t25*(3*t1-24*t3+t22);

  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Hermite5_DDDD( real_type const x, real_type const H, real_type base_DDDD[6] ) {
    real_type const t3 = 360.0*H-720.0*x;
    real_type const t4 = H*H;
    real_type const t5 = t4*t4;
    real_type const t6 = 1/t5;
    real_type const t8 = 1/H;
    real_type const t7 = t8*t6;
    real_type const t10 = 360*x;
    real_type const t16 = 60*x;
    real_type const t19 = t8/t4;
    base_DDDD[0] = t7*t3;
    base_DDDD[1] = -base_DDDD[0];
    base_DDDD[2] = t6*(192*H-t10);
    base_DDDD[3] = t6*(168*H-t10);
    base_DDDD[4] = t19*(36*H-t16);
    base_DDDD[5] = t19*(t16-24*H);

  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Hermite5_DDDDD( real_type, real_type const H, real_type base_DDDDD[6] ) {
    real_type const t1  = H*H;
    real_type const t2  = t1*t1;
    real_type const t3  = 1/t2;
    real_type const t4  = 1/H;
    real_type const t5  = 720.0*t4*t3;
    real_type const t10 = 60.0*t4/t1;
    base_DDDDD[0] = -t5;
    base_DDDDD[1] = t5;
    base_DDDDD[2] = -360.0*t3;
    base_DDDDD[3] = base_DDDDD[2];
    base_DDDDD[4] = -t10;
    base_DDDDD[5] = t10;
  }

  /*
  //   ____  _ _ _
  //  | __ )(_) (_)_ __   ___  __ _ _ __
  //  |  _ \| | | | '_ \ / _ \/ _` | '__|
  //  | |_) | | | | | | |  __/ (_| | |
  //  |____/|_|_|_|_| |_|\___|\__,_|_|
  */

  real_type
  bilinear3(
    real_type const p[4],
    real_type const M[4][4],
    real_type const q[4]
  ) {
    return p[0] * ( M[0][0]*q[0] + M[0][1]*q[1] + M[0][2]*q[2] + M[0][3]*q[3] ) +
           p[1] * ( M[1][0]*q[0] + M[1][1]*q[1] + M[1][2]*q[2] + M[1][3]*q[3] ) +
           p[2] * ( M[2][0]*q[0] + M[2][1]*q[1] + M[2][2]*q[2] + M[2][3]*q[3] ) +
           p[3] * ( M[3][0]*q[0] + M[3][1]*q[1] + M[3][2]*q[2] + M[3][3]*q[3] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  bilinear5(
    real_type const p[6],
    real_type const M[6][6],
    real_type const q[6]
  ) {
    return p[0] * ( M[0][0]*q[0] + M[0][1]*q[1] + M[0][2]*q[2] + M[0][3]*q[3] + M[0][4]*q[4] + M[0][5]*q[5]) +
           p[1] * ( M[1][0]*q[0] + M[1][1]*q[1] + M[1][2]*q[2] + M[1][3]*q[3] + M[1][4]*q[4] + M[1][5]*q[5]) +
           p[2] * ( M[2][0]*q[0] + M[2][1]*q[1] + M[2][2]*q[2] + M[2][3]*q[3] + M[2][4]*q[4] + M[2][5]*q[5]) +
           p[3] * ( M[3][0]*q[0] + M[3][1]*q[1] + M[3][2]*q[2] + M[3][3]*q[3] + M[3][4]*q[4] + M[3][5]*q[5]) +
           p[4] * ( M[4][0]*q[0] + M[4][1]*q[1] + M[4][2]*q[2] + M[4][3]*q[3] + M[4][4]*q[4] + M[4][5]*q[5]) +
           p[5] * ( M[5][0]*q[0] + M[5][1]*q[1] + M[5][2]*q[2] + M[5][3]*q[3] + M[5][4]*q[4] + M[5][5]*q[5]);
  }

  #endif

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  HermiteSpline::build(
    real_type const [], integer,
    real_type const [], integer,
    integer
  ) {
    UTILS_ERROR(
      "HermiteSpline[{}]::build(x,incx,y,incy,n) cannot be used\n", m_name
    );
  }

  using GC_namespace::GC_type;
  using GC_namespace::vec_real_type;

  //!
  //!
  //! Setup a spline using a `GenericContainer`
  //!
  //! - gc("xdata")  vector with the `x` coordinate of the data
  //! - gc("ydata")  vector with the `y` coordinate of the data
  //! - gc("ypdata") vector with the `y` derivative of the data
  //!
  void
  HermiteSpline::setup( GenericContainer const & gc ) {
    /*
    // gc["xdata"]
    // gc["ydata"]
    //
    */
    string const where{ fmt::format("HermiteSpline[{}]::setup( gc ):", m_name ) };
    GenericContainer const & gc_x  { gc("xdata",where) };
    GenericContainer const & gc_y  { gc("ydata",where) };
    GenericContainer const & gc_yp { gc("ypdata",where) };

    vec_real_type x, y, yp;
    {
      string const ff{ fmt::format( "{}, field `xdata'", where ) };
      gc_x.copyto_vec_real( x, ff );
    }
    {
      string const ff{ fmt::format( "{}, field `ydata'", where ) };
      gc_y.copyto_vec_real( y, ff );
    }
    {
      string const ff{ fmt::format( "{}, field `ypdata'", where ) };
      gc_yp.copyto_vec_real( yp, ff );
    }
    this->build( x, y, yp );
  }

}
