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

#include <fstream>

#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat"
#endif

using namespace SplinesLoad;
using namespace std;
using Splines::real_type;
using Splines::integer;

static real_type x[] = {  0,  1,  2,  3,  4,  5 };
static real_type y[] = {  0,  1,  2,  3,  4,  5 };
static real_type z[] = {  10, 10, 10,   10,  10, 11,
                          10, 10, 10.5, 11,  9,  9,
                          11, 12, 13,   8,   7,  9,
                          11, 12, 13,   8,   7,  9,
                          11, 12, 13,   8,   7,  9,
                          11, 12, 13,   8,   7,  9 };

int
main() {
  cout << "\n\nTEST N.2\n\n";

  BiCubicSpline   bc;
  BiQuinticSpline bq;
  BilinearSpline  bl;
  Akima2Dspline   ak;

  real_type X[6], Y[6], Z[6*6];

  std::copy_n( x, 6,   X );
  std::copy_n( y, 6,   Y );
  std::copy_n( z, 6*6, Z );

  //bc.build( x, 1, y, 1, z, 6, 4, 6 );
  //bl.build( x, 1, y, 1, z, 6, 4, 6 );
  //ak.build( x, 1, y, 1, z, 6, 4, 6 );
  bc.build( X, 1, Y, 1, Z, 6, 6, 6 );
  bq.build( X, 1, Y, 1, Z, 6, 6, 6 );
  bl.build( X, 1, Y, 1, Z, 6, 6, 6 );
  ak.build( X, 1, Y, 1, Z, 6, 6, 6 );

  //bl.write_to_stream( cout );

  {
    ofstream file_bl("out/bilinear.txt");
    ofstream file_bc("out/bicubic.txt");
    ofstream file_bq("out/biquintic.txt");
    ofstream file_ak("out/akima2d.txt");

    for ( int i = 0; i <= 100; ++i ) {
      real_type xxx = bc.x_min() + (bc.x_max()-bc.x_min())*i/100.0;
      for ( int j = 0; j <= 100; ++j ) {
        real_type yyy = bc.y_min() + (bc.y_max()-bc.y_min())*j/100.0;
        file_bc << bc(xxx,yyy) << '\t';
        file_bq << bq(xxx,yyy) << '\t';
        file_bl << bl(xxx,yyy) << '\t';
        file_ak << ak(xxx,yyy) << '\t';
      }
      file_bc << '\n';
      file_bq << '\n';
      file_bl << '\n';
      file_ak << '\n';
    }

    file_bc.close();
    file_bq.close();
    file_bl.close();
    file_ak.close();
  }

  {
    ofstream file_bl("out/bilinear_Dx.txt");
    ofstream file_bc("out/bicubic_Dx.txt");
    ofstream file_bq("out/biquintic_Dx.txt");
    ofstream file_ak("out/akima2d_Dx.txt");

    for ( int i = 0; i <= 100; ++i ) {
      real_type xxx = bc.x_min() + (bc.x_max()-bc.x_min())*i/100.0;
      for ( int j = 0; j <= 100; ++j ) {
        real_type yyy = bc.y_min() + (bc.y_max()-bc.y_min())*j/100.0;
        file_bc << bc.Dx(xxx,yyy) << '\t';
        file_bq << bq.Dx(xxx,yyy) << '\t';
        file_bl << bl.Dx(xxx,yyy) << '\t';
        file_ak << ak.Dx(xxx,yyy) << '\t';
      }
      file_bc << '\n';
      file_bq << '\n';
      file_bl << '\n';
      file_ak << '\n';
    }

    file_bc.close();
    file_bq.close();
    file_bl.close();
    file_ak.close();
  }

  {
    ofstream file_bl("out/bilinear_Dxx.txt");
    ofstream file_bc("out/bicubic_Dxx.txt");
    ofstream file_bq("out/biquintic_Dxx.txt");
    ofstream file_ak("out/akima2d_Dxx.txt");

    for ( int i = 0; i <= 100; ++i ) {
      real_type xxx = bc.x_min() + (bc.x_max()-bc.x_min())*i/100.0;
      for ( int j = 0; j <= 100; ++j ) {
        real_type yyy = bc.y_min() + (bc.y_max()-bc.y_min())*j/100.0;
        file_bc << bc.Dxx(xxx,yyy) << '\t';
        file_bq << bq.Dxx(xxx,yyy) << '\t';
        file_bl << bl.Dxx(xxx,yyy) << '\t';
        file_ak << ak.Dxx(xxx,yyy) << '\t';
      }
      file_bc << '\n';
      file_bq << '\n';
      file_bl << '\n';
      file_ak << '\n';
    }

    file_bc.close();
    file_bq.close();
    file_bl.close();
    file_ak.close();
  }

  {
    ofstream file_bl("out/bilinear_Dy.txt");
    ofstream file_bc("out/bicubic_Dy.txt");
    ofstream file_bq("out/biquintic_Dy.txt");
    ofstream file_ak("out/akima2d_Dy.txt");

    for ( int i = 0; i <= 100; ++i ) {
      real_type xxx = bc.x_min() + (bc.x_max()-bc.x_min())*i/100.0;
      for ( int j = 0; j <= 100; ++j ) {
        real_type yyy = bc.y_min() + (bc.y_max()-bc.y_min())*j/100.0;
        file_bc << bc.Dy(xxx,yyy) << '\t';
        file_bq << bq.Dy(xxx,yyy) << '\t';
        file_bl << bl.Dy(xxx,yyy) << '\t';
        file_ak << ak.Dy(xxx,yyy) << '\t';
      }
      file_bc << '\n';
      file_bq << '\n';
      file_bl << '\n';
      file_ak << '\n';
    }

    file_bc.close();
    file_bq.close();
    file_bl.close();
    file_ak.close();
  }

  {
    ofstream file_bl("out/bilinear_Dyy.txt");
    ofstream file_bc("out/bicubic_Dyy.txt");
    ofstream file_bq("out/biquintic_Dyy.txt");
    ofstream file_ak("out/akima2d_Dyy.txt");

    for ( int i = 0; i <= 100; ++i ) {
      real_type xxx = bc.x_min() + (bc.x_max()-bc.x_min())*i/100.0;
      for ( int j = 0; j <= 100; ++j ) {
        real_type yyy = bc.y_min() + (bc.y_max()-bc.y_min())*j/100.0;
        file_bc << bc.Dyy(xxx,yyy) << '\t';
        file_bq << bq.Dyy(xxx,yyy) << '\t';
        file_bl << bl.Dyy(xxx,yyy) << '\t';
        file_ak << ak.Dyy(xxx,yyy) << '\t';
      }
      file_bc << '\n';
      file_bq << '\n';
      file_bl << '\n';
      file_ak << '\n';
    }

    file_bc.close();
    file_bq.close();
    file_bl.close();
    file_ak.close();
  }

  {
    ofstream file_bl("out/bilinear_Dxy.txt");
    ofstream file_bc("out/bicubic_Dxy.txt");
    ofstream file_bq("out/biquintic_Dxy.txt");
    ofstream file_ak("out/akima2d_Dxy.txt");

    for ( int i = 0; i <= 100; ++i ) {
      real_type xxx = bc.x_min() + (bc.x_max()-bc.x_min())*i/100.0;
      for ( int j = 0; j <= 100; ++j ) {
        real_type yyy = bc.y_min() + (bc.y_max()-bc.y_min())*j/100.0;
        file_bc << bc.Dxy(xxx,yyy) << '\t';
        file_bq << bq.Dxy(xxx,yyy) << '\t';
        file_bl << bl.Dxy(xxx,yyy) << '\t';
        file_ak << ak.Dxy(xxx,yyy) << '\t';
      }
      file_bc << '\n';
      file_bq << '\n';
      file_bl << '\n';
      file_ak << '\n';
    }

    file_bc.close();
    file_bq.close();
    file_bl.close();
    file_ak.close();
  }

  cout << "\nALL DONE!\n\n";
}
