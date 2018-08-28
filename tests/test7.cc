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

#include "GenericContainer.hh"
#include "Splines.hh"
#include <fstream>

#ifdef __GCC__
#pragma GCC diagnostic ignored "-Wc++98-compat"
#pragma GCC diagnostic ignored "-Wvla-extension"
#endif
#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wvla-extension"
#endif

using namespace SplinesLoad;
using namespace std;
using Splines::real_type;
using Splines::integer;

// monotone
//static real_type xx[] = { 595, 635, 695, 795, 855, 875, 895, 915, 935, 985, 1035, 1075 };
//static real_type yy[] = { 0.644, 0.652, 0.644, 0.694, 0.907, 1.336, 2.169, 1.598, 0.916, 0.607, 0.603, 0.608 };
//static real_type xx[] = { 0, 1, 2, 3, 5, 6, 7, 9, 10, 11, 12 };
//static real_type yy[] = { 0, 0, 0, 0, 1, 2, 3, 4, 4,   4,  4 };
//static real_type xx[] = { -10, -9, -6, -1, 2, 3, 5, 6, 7, 9, 10, 11, 12 };
//static real_type yy[] = { 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 4, 4, 4 };
//static real_type yy[] = {  595, 635, 695, 795, 855, 875, 895, 915, 935, 985, 1035, 1075 };
static real_type xx[] = { 0, 1, 2, 3, 4 };
static real_type yy[] = { 0, 1, 2, 3, 4 };
//static real_type yy[] = { 1, 1, 1, 1 };
//static real_type yy[] = { 0, 1, 4, 9, 16, 64, 100 };
static integer  n = sizeof(xx)/sizeof(xx[0]);

#define SAVE(S) \
cout << #S": n = " << n << '\n'; \
S.clear(); \
S.reserve(n); \
for ( integer i = 0; i < n; ++i ) S.pushBack(xx[i],yy[i]); \
S.build(); \
cout << #S": xMin    = " << S.xMin() << '\n'; \
cout << #S": xMax    = " << S.xMax() << '\n'; \
cout << #S": xx[0]   = " << xx[0]    << '\n'; \
cout << #S": xx[end] = " << xx[n-1]  << '\n'; \
file_##S << "x\ty\tDy\tDDy\n"; \
for ( real_type x = xmin-(xmax-xmin)*0.01; x <= xmax+(xmax-xmin)*0.01; x += (xmax-xmin)/1000 ) \
  file_##S << x << '\t' << S(x) << '\t' << S.D(x) << '\t' << S.DD(x) << '\n'; \
file_##S.close()

int
main() {

  cout << "\n\nTEST N.7\n\n";

  real_type xmin = xx[0];
  real_type xmax = xx[n-1];

  PchipSpline pc;
  //CubicSpline pc;
  BSpline<3>  bs;
  ofstream file_pc("out/PchipSpline.txt");
  ofstream file_bs("out/BSpline.txt");
  ofstream file_bs_base("out/BSpline_base.txt");
  ofstream file_bs_base_D("out/BSpline_base_D.txt");
  ofstream file_bs_base_DD("out/BSpline_base_DD.txt");
  ofstream file_bs_base_DDD("out/BSpline_base_DDD.txt");
  ofstream file_bs2("out/DATA.txt");

  SAVE(pc);
  SAVE(bs);

  file_bs_base     << "x";
  file_bs_base_D   << "x";
  file_bs_base_DD  << "x";
  file_bs_base_DDD << "x";
  for ( integer i = 0; i < n; ++i ) {
    file_bs_base     << "\tb" << i;
    file_bs_base_D   << "\tb" << i;
    file_bs_base_DD  << "\tb" << i;
    file_bs_base_DDD << "\tb" << i;
  }
  file_bs_base     << '\n';
  file_bs_base_D   << '\n';
  file_bs_base_DD  << '\n';
  file_bs_base_DDD << '\n';

  for ( real_type x = xmin-(xmax-xmin)*0.5;
        x <= xmax+(xmax-xmin)*0.5;
        x += (xmax-xmin)/1000 ) {
    file_bs_base     << x;
    file_bs_base_D   << x;
    file_bs_base_DD  << x;
    file_bs_base_DDD << x;
    vector<real_type> vals, vals_D, vals_DD, vals_DDD;
    vals.resize( size_t(n) );
    vals_D.resize( size_t(n) );
    vals_DD.resize( size_t(n) );
    vals_DDD.resize( size_t(n) );
    bs.bases     ( x, &vals.front() );
    bs.bases_D   ( x, &vals_D.front() );
    bs.bases_DD  ( x, &vals_DD.front() );
    bs.bases_DDD ( x, &vals_DDD.front() );
    for ( size_t i = 0; i < size_t(n); ++i ) {
      file_bs_base     << '\t' << vals[i];
      file_bs_base_D   << '\t' << vals_D[i];
      file_bs_base_DD  << '\t' << vals_DD[i];
      file_bs_base_DDD << '\t' << vals_DDD[i];
    }
    file_bs_base     << '\n';
    file_bs_base_D   << '\n';
    file_bs_base_DD  << '\n';
    file_bs_base_DDD << '\n';
  }
  file_bs_base.close();
  file_bs_base_D.close();
  file_bs_base_DD.close();
  file_bs_base_DDD.close();

  file_bs2 << "x\ty\tz\n";
  for ( integer i = 0; i < n; ++i )
    file_bs2 << xx[i] << '\t' << yy[i] << "\t0\n";
  file_bs2.close();

  cout << "\nALL DONE!\n\n";
}
