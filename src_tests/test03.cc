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
static real_type xx5[] = { 0.11, 0.12, 0.15, 0.16 };
static real_type yy5[] = { 0.0003, 0.0003, 0.0004, 0.0004 };

// monotone
static real_type xx6[] = { 0, 1, 2, 3, 4 };
static real_type yy6[] = { 0, 1, 1.1, 2.0, 2.1 };

static integer nn[] = { 11, 11, 11, 9, 12, 4, 5 };

int
main() {

  cout << "\n\nTEST N.3\n\n";

  SplineSet ss;
  ofstream  file, file_D;

  for ( integer k{0}; k < 7; ++ k ) {
    real_type * xx{ nullptr };
    real_type * yy{ nullptr };
    switch ( k ) {
      case 0: xx = xx0; yy = yy0; break;
      case 1: xx = xx1; yy = yy1; break;
      case 2: xx = xx2; yy = yy2; break;
      case 3: xx = xx3; yy = yy3; break;
      case 4: xx = xx4; yy = yy4; break;
      case 5: xx = xx5; yy = yy5; break;
      case 6: xx = xx6; yy = yy6; break;
    }
    string fname;
    fname = fmt::format( "out/SplineSet{}.txt", k );
    file.open(fname.data());

    fname = fmt::format( "out/SplineSet{}_D.txt", k );
    file_D.open(fname.data());

    real_type xmin { xx[0]       };
    real_type xmax { xx[nn[k]-1] };
    integer   nspl { 7           };
    integer   npts { nn[k]       };

    char const *headers[] = {
      "SPLINE_CONSTANT",
      "SPLINE_LINEAR",
      "SPLINE_CUBIC",
      "SPLINE_AKIMA",
      "SPLINE_BESSEL",
      "SPLINE_PCHIP",
      "SPLINE_QUINTIC"
    };

    constexpr SplineType1D stype[] = {
       SplineType1D::CONSTANT,
       SplineType1D::LINEAR,
       SplineType1D::CUBIC,
       SplineType1D::AKIMA,
       SplineType1D::BESSEL,
       SplineType1D::PCHIP,
       SplineType1D::QUINTIC
    };

    //bool const rp_policy[] = { true, true, true, true, true, true, true };
    real_type const *Y[]{ yy, yy, yy, yy, yy, yy, yy, yy };

    ss.build( nspl, npts, headers, stype, xx, Y );
    ss.info(cout);

    file   << "x";
    file_D << "x";
    for ( integer i{0}; i < nspl; ++i ) {
      file   << '\t' << ss.header(i);
      file_D << '\t' << ss.header(i);
    }
    file   << '\n';
    file_D << '\n';
    for ( real_type x = xmin; x <= xmax; x += (xmax-xmin)/1000 ) {
      file   << x;
      file_D << x;
      for ( integer i{0}; i < nspl; ++i ) {
        file   << '\t' << ss(x,i);
        file_D << '\t' << ss.D(x,i);
      }
      file   << '\n';
      file_D << '\n';
    }
    file.close();
    file_D.close();
  }

  cout << "\nALL DONE!\n\n";
}
