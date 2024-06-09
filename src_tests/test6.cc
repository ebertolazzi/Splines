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
#include "Utils_fmt.hh"

#include <GenericContainer/GenericContainer.hh>
#include <fstream>

#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat"
#endif

using namespace SplinesLoad;
using namespace std;
using Splines::real_type;
using Splines::integer;

// monotone
static real_type xx[] = { 0, 0.9, 2.1, 3, 4.5 };
static real_type yy[] = { 0, 1, 1.99, 2.0, 2.1 };

static integer npt = 5;

int
main() {

  cout << "\n\nTEST N.6\n\n";

  SplineSet ss;
  ofstream  file, file_D, fileR, fileR_D;

  file.open("out/SplineSet.txt");
  file_D.open("out/SplineSet_D.txt");
  fileR.open("out/SplineSetR.txt");
  fileR_D.open("out/SplineSetR_D.txt");

  real_type xmin = xx[0];
  real_type xmax = xx[npt-1];

  integer  nspl = 7;
  integer  npts = npt;
  real_type val[7], val_D[7];

  char const *headers[] = {
    "constant",
    "linear",
    "cubic",
    "akima",
    "bessel",
    "pchip",
    "quintic"
  };

  real_type const *Y[] = { yy, yy, yy, yy, yy, yy, yy, yy };

  GC::GenericContainer gc;

  GC::vec_string_type & t = gc["spline_type"].set_vec_string();
  GC::vec_string_type & h = gc["headers"].set_vec_string();
  t.resize( size_t(nspl) );
  h.resize( size_t(nspl) );
  std::copy_n( headers, nspl, h.begin() );
  std::copy_n( headers, nspl, t.begin() );

  GC::vector_type & data = gc["ydata"].set_vector();
  data.resize( size_t(nspl) );
  for ( integer i = 0; i < nspl; ++i ) {
    GC::GenericContainer & di = data[size_t(i)];
    GC::vec_real_type    & v  = di.set_vec_real();
    if ( i == 0 ) {
      // spline constante ha 1 punto in meno
      v.resize( size_t(npts-1) );
      std::copy_n( Y[i], npts-1, v.begin() );
    } else {
      v.resize( size_t(npts) );
      std::copy_n( Y[i], npts, v.begin() );
    }
  }

  GC::vec_real_type & xdata = gc["xdata"].set_vec_real();
  xdata.resize( size_t(npts) );
  std::copy_n( xx, npts, xdata.begin() );

  gc.print(cout);
  ss.build( gc ); // nspl, npts, headers, stype, xx, Y, Yp );
  ss.info(cout);

  file   << "x";
  file_D << "x";
  for ( integer i = 0; i < nspl; ++i ) {
    file   << '\t' << ss.header(i);
    file_D << '\t' << ss.header(i);
  }
  file   << '\n';
  file_D << '\n';
  for ( real_type x = xmin; x <= xmax; x += (xmax-xmin)/1000 ) {
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
  xmax = yy[npt-1];

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
