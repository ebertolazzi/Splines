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
#include <fstream>

#if defined(__GCC__) || defined(__GNUC__)
#pragma GCC diagnostic ignored "-Wc++98-compat"
#endif
#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat"
#endif

using namespace SplinesLoad;
using namespace std;
using Splines::real_type;
using Splines::integer;

// monotone
static real_type xx[] = { 0, 0.9, 2.1, 3, 4.5 };
static real_type yy[] = { 0, 1, 1.1, 2.0, 2.1 };

static integer n = 5;

int
main() {

  cout << "\n\nTEST N.5\n\n";

  SplineSet ss;
  ofstream  file, file_D, fileR, fileR_D;

  file.open("out/SplineSet.txt");
  file_D.open("out/SplineSet_D.txt");
  fileR.open("out/SplineSetR.txt");
  fileR_D.open("out/SplineSetR_D.txt");
  
  real_type xmin = xx[0];
  real_type xmax = xx[n-1];

  integer  nspl = 7;
  integer  npts = n;
  real_type val[8], val_D[8];

  char const *headers[] = {
    "SPLINE_CONSTANT",
    "SPLINE_LINEAR",
    "SPLINE_CUBIC",
    "SPLINE_AKIMA",
    "SPLINE_BESSEL",
    "SPLINE_PCHIP",
    "SPLINE_QUINTIC"
  };
    
  SplineType const stype[] = {
    Splines::CONSTANT_TYPE,
    Splines::LINEAR_TYPE,
    Splines::CUBIC_TYPE,
    Splines::AKIMA_TYPE,
    Splines::BESSEL_TYPE,
    Splines::PCHIP_TYPE,
    Splines::QUINTIC_TYPE
  };

  std::vector<real_type> YpZero;
  YpZero.resize( size_t(npts) );
  std::fill(YpZero.begin(), YpZero.end(), 0 );

  real_type const *Y[]  = { yy, yy, yy, yy, yy, yy, yy, yy };

  ss.build( nspl, npts, headers, stype, xx, Y );
  ss.info(cout);

  file   << "x";
  file_D << "x";
  for ( integer i = 0; i < nspl; ++i ) {
    file   << '\t' << ss.header(i);
    file_D << '\t' << ss.header(i);
  }
  file   << '\n';
  file_D << '\n';
  for ( real_type x = xmin-(xmax-xmin)*0.01; x <= xmax+(xmax-xmin)*0.01; x += (xmax-xmin)/1000 ) {
    file   << x;
    file_D << x;
    ss.eval( x, val );
    ss.eval_D( x, val_D );
    for ( integer i = 0; i < nspl; ++i ) {
      file   << '\t' << val[i];
      file_D << '\t' << val_D[i];
    }
    file   << '\n';
    file_D << '\n';
  }
  file.close();
  file_D.close();


  xmin = yy[0];
  xmax = yy[n-1];

  fileR   << "x";
  fileR_D << "x";
  for ( integer i = 0; i < nspl; ++i ) {
    fileR   << '\t' << ss.header(i);
    fileR_D << '\t' << ss.header(i);
  }
  fileR   << '\n';
  fileR_D << '\n';

  for ( real_type x = xmin; x <= xmax; x += (xmax-xmin)/1000 ) {
    fileR   << x;
    fileR_D << x;
    ss.eval2( 5, x, val );
    ss.eval2_D( 5, x, val_D );
    for ( integer i = 0; i < nspl; ++i ) {
      fileR   << '\t' << val[i];
      fileR_D << '\t' << val_D[i];
    }
    fileR   << '\n';
    fileR_D << '\n';
  }
  fileR.close();
  fileR_D.close();

  cout << "\nALL DONE!\n\n";
}
