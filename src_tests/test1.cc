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

#include <fstream>

#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat"
#endif

using namespace SplinesLoad;
using namespace std;
using Splines::real_type;
using Splines::integer;

// Test problem for Akima interpolation
// Ref. : Hiroshi Akima, Journal of the ACM, Vol. 17, No. 4, October 1970, pages 589-602.

static real_type xx0[] = {  0,  1,  2,  3,  4,  5,    6,  7,  8,  9, 10 };
static real_type yy0[] = { 10, 10, 10, 10, 10, 10, 10.5, 15, 50, 60, 85 };

static real_type xx1[] = {  0,  1,  3,  4,  6,  7,    9, 10, 12, 13, 15 };
static real_type yy1[] = { 10, 10, 10, 10, 10, 10, 10.5, 15, 50, 60, 85 };

static real_type xx2[] = {  0,  2,  3,  5,  6,  8,    9, 11, 12, 14, 15 };
static real_type yy2[] = { 10, 10, 10, 10, 10, 10, 10.5, 15, 50, 60, 85 };

// RPN 14
static real_type xx3[] = { 7.99, 8.09,       8.19,       8.7,      9.2,      10,       12,       15,       20       };
static real_type yy3[] = { 0,    2.76429e-5, 4.37498e-2, 0.169183, 0.469428, 0.943740, 0.998636, 0.999919, 0.999994 };

// Titanium
static real_type xx4[] = { 595,   635,   695,   795,   855,   875,   895,   915,   935,   985,   1035,  1075  };
static real_type yy4[] = { 0.644, 0.652, 0.644, 0.694, 0.907, 1.336, 2.169, 1.598, 0.916, 0.607, 0.603, 0.608 };

// toolpath
static real_type xx5[] = { 0.11,   0.12,   0.15,   0.16   };
static real_type yy5[] = { 0.0003, 0.0003, 0.0004, 0.0004 };

static integer nn[] = { 11, 11, 11, 9, 12, 4 };

int
main() {
  fmt::print("\n\nTEST N.1\n\n");

  LinearSpline   li;
  ConstantSpline co;
  AkimaSpline    ak;
  CubicSpline    cs;
  BesselSpline   be;
  PchipSpline    pc;
  QuinticSpline  qs;
  ofstream       file_li;
  ofstream       file_co;
  ofstream       file_ak;
  ofstream       file_cs;
  ofstream       file_be;
  ofstream       file_pc;
  ofstream       file_qs;

  for ( integer k{0}; k < 6; ++ k ) {
    fmt::print( "\n\nk = {}\n", k );
    real_type * xx = nullptr, * yy = nullptr;
    switch ( k ) {
      case 0: xx = xx0; yy = yy0; break;
      case 1: xx = xx1; yy = yy1; break;
      case 2: xx = xx2; yy = yy2; break;
      case 3: xx = xx3; yy = yy3; break;
      case 4: xx = xx4; yy = yy4; break;
      case 5: xx = xx5; yy = yy5; break;
    }
    string fname;
    fname = fmt::format( "out/Linear{}.txt",   k); file_li.open(fname.data());
    fname = fmt::format( "out/Constant{}.txt", k); file_co.open(fname.data());
    fname = fmt::format( "out/Akima{}.txt",    k); file_ak.open(fname.data());
    fname = fmt::format( "out/Cubic{}.txt",    k); file_cs.open(fname.data());
    fname = fmt::format( "out/Bessel{}.txt",   k); file_be.open(fname.data());
    fname = fmt::format( "out/Pchip{}.txt",    k); file_pc.open(fname.data());
    fname = fmt::format( "out/Quintic{}.txt",  k); file_qs.open(fname.data());
    real_type xmin = xx[0];
    real_type xmax = xx[nn[k]-1];

    #define SAVE(NAME,S)                                                    \
    fmt::print( "\n\n\n{}\n\n\n", NAME );                                   \
    fmt::print( #S": n[k] = {}\n", nn[k] );                                 \
    S.clear();                                                              \
    S.reserve(nn[k]);                                                       \
    for ( integer i{0}; i < integer(nn[k]); ++i )                           \
      S.push_back(xx[i],yy[i]);                                             \
    S.build(); /*( xx, yy, n[k] );*/                                        \
    {                                                                       \
      integer   i_min_pos;                                                  \
      real_type x_min_pos;                                                  \
      real_type y_min;                                                      \
      integer   i_max_pos;                                                  \
      real_type x_max_pos;                                                  \
      real_type y_max;                                                      \
      S.y_min_max(                                                          \
        i_min_pos, x_min_pos, y_min, i_max_pos, x_max_pos, y_max            \
      );                                                                    \
      fmt::print( "i_min_pos = {}\n", i_min_pos );                          \
      fmt::print( "x_min_pos = {}\n", x_min_pos );                          \
      fmt::print( "y_min     = {}\n", y_min );                              \
      fmt::print( "i_max_pos = {}\n", i_max_pos );                          \
      fmt::print( "x_max_pos = {}\n", x_max_pos );                          \
      fmt::print( "y_max     = {}\n", y_max );                              \
    }                                                                       \
    fmt::print( #S": xMin    = {}\n", S.x_min() );                          \
    fmt::print( #S": xMax    = {}\n", S.x_max() );                          \
    fmt::print( #S": xx[0]   = {}\n", xx[0] );                              \
    fmt::print( #S": xx[end] = {}\n", xx[nn[k]-1] );                        \
    file_##S << "x\ty\tDy\tDDy\n";                                          \
    for ( real_type x = xmin; x <= xmax; x += (xmax-xmin)/1000 )            \
      fmt::print( file_##S, "{}\t{}\t{}\t{}\n", x, S(x), S.D(x), S.DD(x) ); \
    file_##S.close()

    SAVE("LinearSpline",   li);
    SAVE("ConstantSpline", co);
    SAVE("AkimaSpline",    ak);
    SAVE("CubicSpline",    cs);
    SAVE("BesselSpline",   be);
    SAVE("PchipSpline",    pc);
    SAVE("QuinticSpline",  qs);
  }
  cout << "\nALL DONE!\n\n";
}
