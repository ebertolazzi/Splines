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
#include <cmath>
#include <iomanip>
/**
 * 
 */

namespace Splines {

  using namespace std; // load standard namspace
  
  // Statement Function definitions
  inline
  real_type
  Extrapolate2(
    real_type X1,
    real_type X2,
    real_type Z0,
    real_type Z1
  ) {
    return (Z1-Z0)*X2/X1 + Z0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  inline
  real_type
  Extrapolate3(
    real_type X1,
    real_type X2,
    real_type X3,
    real_type Z0,
    real_type Z1,
    real_type Z2
  ) {
    return ( (Z2-Z0) * (X3-X1)/X2 - (Z1-Z0) * (X3-X2)/X1 ) * (X3/(X2-X1)) + Z0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  estimate(
    real_type   z0,
    real_type   z1,
    real_type   z2,
    real_type   z3,
    real_type   x1,
    real_type   x2,
    real_type   x3,
    real_type & c1,
    real_type & c2,
    real_type & c3,
    real_type & SX,
    real_type & SXX,
    real_type & b0,
    real_type & b1,
    real_type   RF[2],
    real_type   RI[2]
  ) {

    // Primary estimate of partial derivative zx
    // as the coefficient of the bicubic polynomial.
    c1 = x2*x3/((x1-x2)*(x1-x3));
    c2 = x3*x1/((x2-x3)*(x2-x1));
    c3 = x1*x2/((x3-x1)*(x3-x2));
    real_type primary_estimate = c1*(z1-z0)/x1 + c2*(z2-z0)/x2 + c3*(z3-z0)/x3;

    // Volatility factor and distance factor in the x direction
    // for the primary estimate of zx.
              SX  = x1 + x2 + x3;
    real_type SZ  = z0 + z1 + z2 + z3;
              SXX = x1*x1 + x2*x2 + x3*x3;
    real_type SXZ = x1*z1 + z2*z2 + x3*z3;
    real_type DNM = 4.0*SXX - SX*SX;
              b0  = (SXX*SZ-SX*SXZ)/DNM;
              b1  = (4.0*SXZ-SX*SZ)/DNM;
    real_type dz0 = z0 - b0;
    real_type dz1 = z1 - (b0+b1*x1);
    real_type dz2 = z2 - (b0+b1*x2);
    real_type dz3 = z3 - (b0+b1*x3);
    real_type volatility_factor = dz0*dz0 + dz1*dz1 + dz2*dz2 + dz3*dz3;
    // epsi value used to decide whether or not
    // the volatility factor is essentially zero.
    real_type epsi = (z0*z0+z1*z1+z2*z2+z3*z3)*1.0E-12;
    // Accumulates the weighted primary estimates and their weights.
    if ( volatility_factor > epsi ) { // finite weight.
      real_type WT = 1.0/ (volatility_factor*SXX);
      RF[1] += WT*primary_estimate;
      RF[0] += WT;
    } else { // infinite weight.
      RI[1] += primary_estimate;
      RI[0] += 1.0;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // 0 1 2 3 (4) 5 6 7 8
  static
  void
  AkimaSmooth(
    real_type const X[9], integer imin, integer imax,
    real_type const Y[9], integer jmin, integer jmax,
    real_type const Z[9][9],
    real_type & DX,
    real_type & DY,
    real_type & DXY
  ) {

    integer stencil[3][4] = { { -3, -2, -1, 1 },
                              { -2, -1,  1, 2 },
                              { -1,  1,  2, 3 } };

    real_type by0[4], by1[4], CY1A[4], CY2A[4], CY3A[4], SYA[4], SYYA[4];

    real_type X0  = X[4];
    real_type Y0  = Y[4];
    real_type Z00 = Z[4][4];

    // Initial setting
    real_type DXF[2]  = {0,0};
    real_type DXI[2]  = {0,0};
    real_type DYF[2]  = {0,0};
    real_type DYI[2]  = {0,0};
    real_type DXYF[2] = {0,0};
    real_type DXYI[2] = {0,0};

    for ( integer k = 0; k < 4; ++k ) {
      int j1 = 4 + stencil[0][k];
      int j2 = 4 + stencil[1][k];
      int j3 = 4 + stencil[2][k];
      if ( j1 >= int(jmin) && j3 <= int(jmax) )
        estimate( Z00, Z[4][j1], Z[4][j2], Z[4][j3],
                  Y[j1] - Y0, Y[j2] - Y0, Y[j3] - Y0,
                  CY1A[k], CY2A[k], CY3A[k], SYA[k], SYYA[k],
                  by0[k], by1[k], DYF, DYI );
    }

    for ( integer kx = 0; kx < 4; ++kx ) {
      int i1 = 4 + stencil[0][kx];
      int i2 = 4 + stencil[1][kx];
      int i3 = 4 + stencil[2][kx];

      if ( i1 < int(imin) || i3 > int(imax) ) continue;

      real_type X1  = X[i1] - X0;
      real_type X2  = X[i2] - X0;
      real_type X3  = X[i3] - X0;
      real_type Z10 = Z[i1][4];
      real_type Z20 = Z[i2][4];
      real_type Z30 = Z[i3][4];
      real_type CX1, CX2, CX3, SX, SXX, B00X, B10;
      estimate( Z00, Z10, Z20, Z30,
                X1, X2, X3,
                CX1, CX2, CX3,
                SX, SXX,
                B00X, B10, DXF, DXI );

      for ( int ky = 0; ky < 4; ++ky ) {
        int j1 = 4 + stencil[0][ky];
        int j2 = 4 + stencil[1][ky];
        int j3 = 4 + stencil[2][ky];
        if ( j1 < int(jmin) || j3 > int(jmax) ) continue;

        real_type Y1   = Y[j1] - Y0;
        real_type Y2   = Y[j2] - Y0;
        real_type Y3   = Y[j3] - Y0;
        real_type CY1  = CY1A[ky];
        real_type CY2  = CY2A[ky];
        real_type CY3  = CY3A[ky];
        real_type SY   = SYA[ky];
        real_type SYY  = SYYA[ky];
        real_type B00Y = by0[ky];
        real_type B01  = by1[ky];

        real_type Z01 = Z[4][j1],  Z02 = Z[4][j2],  Z03 = Z[4][j3];
        real_type Z11 = Z[i1][j1], Z12 = Z[i1][j2], Z13 = Z[i1][j3];
        real_type Z21 = Z[i2][j1], Z22 = Z[i2][j2], Z23 = Z[i2][j3];
        real_type Z31 = Z[i3][j1], Z32 = Z[i3][j2], Z33 = Z[i3][j3];

        // Primary estimate of partial derivative zxy
        // as the coefficient of the bicubic polynomial.
        real_type DZXY11 = (Z11-Z10-Z01+Z00)/(X1*Y1);
        real_type DZXY12 = (Z12-Z10-Z02+Z00)/(X1*Y2);
        real_type DZXY13 = (Z13-Z10-Z03+Z00)/(X1*Y3);
        real_type DZXY21 = (Z21-Z20-Z01+Z00)/(X2*Y1);
        real_type DZXY22 = (Z22-Z20-Z02+Z00)/(X2*Y2);
        real_type DZXY23 = (Z23-Z20-Z03+Z00)/(X2*Y3);
        real_type DZXY31 = (Z31-Z30-Z01+Z00)/(X3*Y1);
        real_type DZXY32 = (Z32-Z30-Z02+Z00)/(X3*Y2);
        real_type DZXY33 = (Z33-Z30-Z03+Z00)/(X3*Y3);
        real_type PEZXY  = CX1* (CY1*DZXY11+CY2*DZXY12+CY3*DZXY13) +
                           CX2* (CY1*DZXY21+CY2*DZXY22+CY3*DZXY23) +
                           CX3* (CY1*DZXY31+CY2*DZXY32+CY3*DZXY33);
        // Volatility factor and distance factor
        // in the x and y directions for the primary estimate of zxy.
        real_type B00   = (B00X+B00Y)/2.0;
        real_type SXY   = SX*SY;
        real_type SXXY  = SXX*SY;
        real_type SXYY  = SX*SYY;
        real_type SXXYY = SXX*SYY;
        real_type SXYZ  = X1 * (Y1*Z11+Y2*Z12+Y3*Z13) +
                          X2 * (Y1*Z21+Y2*Z22+Y3*Z23) +
                          X3 * (Y1*Z31+Y2*Z32+Y3*Z33);
        real_type B11  = (SXYZ-B00*SXY-B10*SXXY-B01*SXYY)/SXXYY;
        real_type DZ00 = Z00 - B00;
        real_type DZ01 = Z01 - (B00+B01*Y1);
        real_type DZ02 = Z02 - (B00+B01*Y2);
        real_type DZ03 = Z03 - (B00+B01*Y3);
        real_type DZ10 = Z10 - (B00+B10*X1);
        real_type DZ11 = Z11 - (B00+B01*Y1+X1*(B10+B11*Y1));
        real_type DZ12 = Z12 - (B00+B01*Y2+X1*(B10+B11*Y2));
        real_type DZ13 = Z13 - (B00+B01*Y3+X1*(B10+B11*Y3));
        real_type DZ20 = Z20 - (B00+B10*X2);
        real_type DZ21 = Z21 - (B00+B01*Y1+X2*(B10+B11*Y1));
        real_type DZ22 = Z22 - (B00+B01*Y2+X2*(B10+B11*Y2));
        real_type DZ23 = Z23 - (B00+B01*Y3+X2*(B10+B11*Y3));
        real_type DZ30 = Z30 - (B00+B10*X3);
        real_type DZ31 = Z31 - (B00+B01*Y1+X3*(B10+B11*Y1));
        real_type DZ32 = Z32 - (B00+B01*Y2+X3*(B10+B11*Y2));
        real_type DZ33 = Z33 - (B00+B01*Y3+X3*(B10+B11*Y3));
        real_type volatility_factor = DZ00*DZ00 + DZ01*DZ01 +
                                      DZ02*DZ02 + DZ03*DZ03 +
                                      DZ10*DZ10 + DZ11*DZ11 +
                                      DZ12*DZ12 + DZ13*DZ13 +
                                      DZ20*DZ20 + DZ21*DZ21 +
                                      DZ22*DZ22 + DZ23*DZ23 +
                                      DZ30*DZ30 + DZ31*DZ31 +
                                      DZ32*DZ32 + DZ33*DZ33;
        real_type epsi = (Z00*Z00+Z01*Z01+Z02*Z02+Z03*Z03+
                          Z10*Z10+Z11*Z11+Z12*Z12+Z13*Z13+
                          Z20*Z20+Z21*Z21+Z22*Z22+Z23*Z23+
                          Z30*Z30+Z31*Z31+Z32*Z32+Z33*Z33)*1.0E-12;
        // Accumulates the weighted primary estimates of zxy and their weights.
        if ( volatility_factor > epsi ) { // finite weight.
          real_type WT = 1 / (volatility_factor*SXX*SYY);
          DXYF[1] += WT*PEZXY;
          DXYF[0] += WT;
        } else { // infinite weight.
          DXYI[1] += PEZXY;
          DXYI[0] += 1.0;
        }
      }
    }
    DX  = DXI[0]  < 0.5 ? DXF[1]/DXF[0]   : DXI[1]/DXI[0];
    DY  = DYI[0]  < 0.5 ? DYF[1]/DYF[0]   : DYI[1]/DYI[0];
    DXY = DXYI[0] < 0.5 ? DXYF[1]/DXYF[0] : DXYI[1]/DXYI[0];
  }
    
  /*
   * This subroutine estimates three partial derivatives, zx, zy, and
   * zxy, of a bivariate function, z(x,y), on a rectangular grid in
   * the x-y plane.  It is based on the revised Akima method that has
   * the accuracy of a bicubic polynomial.
   */
  void
  Akima2Dspline::makeSpline() {
    this->DX.resize(Z.size());
    this->DY.resize(Z.size());
    this->DXY.resize(Z.size());
    // calcolo derivate
    size_t nx = size_t(X.size());
    size_t ny = size_t(Y.size());
    
    std::fill(DX.begin(),DX.end(),0);
    std::fill(DY.begin(),DY.end(),0);
    std::fill(DXY.begin(),DXY.end(),0);
    
    real_type x_loc[9], y_loc[9], z_loc[9][9];

    for ( size_t i0 = 0; i0 < nx; ++i0 ) {
      size_t imin = 4  > i0   ? 4-i0      : 0;
      size_t imax = nx < 5+i0 ? 3+(nx-i0) : 8;

      for ( size_t i = imin; i <= imax; ++i ) x_loc[i] = X[i+i0-4]-X[i0];

      for ( size_t j0 = 0; j0 < ny; ++j0 ) {
        size_t jmin = 4 > j0    ? 4-j0      : 0;
        size_t jmax = ny < 5+j0 ? 3+(ny-j0) : 8;

        for ( size_t j = jmin; j <= jmax; ++j ) y_loc[j] = Y[j+j0-4]-Y[j0];

        for ( size_t i = imin; i <= imax; ++i )
          for ( size_t j = jmin; j <= jmax; ++j )
            z_loc[i][j] = Z[size_t(ipos_C(integer(i+i0-4),
                                          integer(j+j0-4),
                                          integer(ny)))];

        // if not enough points, extrapolate
        size_t iadd = 0, jadd = 0;
        if ( imax < 3+imin ) {
          x_loc[imin-1] = 2*x_loc[imin] - x_loc[imax];
          x_loc[imax+1] = 2*x_loc[imax] - x_loc[imin];
          iadd = 1;
          if ( imax == 1+imin ) {
            real_type x0 = x_loc[imin];
            real_type x1 = x_loc[imax];
            for ( size_t j = jmin; j <= jmax; ++j ) {
              real_type z0 = z_loc[imin][j];
              real_type z1 = z_loc[imax][j];
              z_loc[imin-1][j] = Extrapolate2( x0-x1, x_loc[imin-1]-x1, z1, z0 );
              z_loc[imax+1][j] = Extrapolate2( x1-x0, x_loc[imax+1]-x0, z0, z1 );
            }
          } else {
            real_type x0 = x_loc[imin];
            real_type x1 = x_loc[imin+1];
            real_type x2 = x_loc[imax];
            for ( size_t j = jmin; j <= jmax; ++j ) {
              real_type z0 = z_loc[imin][j];
              real_type z1 = z_loc[imin+1][j];
              real_type z2 = z_loc[imax][j];
              z_loc[imin-1][j] = Extrapolate3( x1-x2, x0-x2, x_loc[imin-1]-x1, z2, z1, z0 );
              z_loc[imax+1][j] = Extrapolate3( x1-x0, x2-x0, x_loc[imax+1]-x0, z0, z1, z2 );
            }
          }
        }
        if ( jmax < 3+jmin ) {
          y_loc[jmin-1] = 2*y_loc[jmin] - y_loc[jmax];
          y_loc[jmax+1] = 2*y_loc[jmax] - y_loc[jmin];
          jadd = 1;
          if ( jmax-jmin == 1 ) {
            real_type y0 = y_loc[jmin];
            real_type y1 = y_loc[jmax];
            for ( size_t i = imin-iadd; i <= imax+iadd; ++i ) {
              real_type z0 = z_loc[i][jmin];
              real_type z1 = z_loc[i][jmax];
              z_loc[i][jmin-1] = Extrapolate2( y0-y1, y_loc[jmin-1]-y1, z1, z0 );
              z_loc[i][jmax+1] = Extrapolate2( y1-y0, y_loc[jmax+1]-y0, z0, z1 );
            }
          } else {
            real_type y0 = y_loc[jmin];
            real_type y1 = y_loc[jmin+1];
            real_type y2 = y_loc[jmax];
            for ( size_t i = imin-iadd; i <= imax+iadd; ++i ) {
              real_type z0 = z_loc[i][jmin];
              real_type z1 = z_loc[i][jmin+1];
              real_type z2 = z_loc[i][jmax];
              z_loc[i][imin-1] = Extrapolate3( y1-y2, y0-y2, y_loc[jmin-1]-y1, z2, z1, z0 );
              z_loc[i][imax+1] = Extrapolate3( y1-y0, y2-y0, y_loc[jmax+1]-y0, z0, z1, z2 );
            }
          }
        }

        size_t i0j0 = size_t(ipos_C(integer(i0),integer(j0)));

        AkimaSmooth( x_loc, integer(imin-iadd), integer(imax+iadd),
                     y_loc, integer(jmin-jadd), integer(jmax+jadd),
                     z_loc, DX[i0j0], DY[i0j0], DXY[i0j0] );
      }
    }
  }

  void
  Akima2Dspline::writeToStream( ostream_type & s ) const {
    integer ny = integer(Y.size());
    s << "Nx = " << X.size() << " Ny = " << Y.size() << '\n';
    for ( integer i = 1; i < integer(X.size()); ++i ) {
      for ( integer j = 1; j < integer(Y.size()); ++j ) {
        size_t i00 = size_t( ipos_C(i-1,j-1,ny) );
        size_t i10 = size_t( ipos_C(i,j-1,ny) );
        size_t i01 = size_t( ipos_C(i-1,j,ny) );
        size_t i11 = size_t( ipos_C(i,j,ny) );
        s << "patch (" << i << "," << j
          << ")\n DX = " << setw(10) << left << X[size_t(i)]-X[size_t(i-1)]
          <<    " DY = " << setw(10) << left << Y[size_t(j)]-Y[size_t(j-1)]
          << "\n Z00  = " << setw(10) << left << Z[i00]
          <<   " Z01  = " << setw(10) << left << Z[i01]
          <<   " Z10  = " << setw(10) << left << Z[i10]
          <<   " Z11  = " << setw(10) << left << Z[i11]
          << "\n Dx00 = " << setw(10) << left << DX[i00]
          <<   " Dx01 = " << setw(10) << left << DX[i01]
          <<   " Dx10 = " << setw(10) << left << DX[i10]
          <<   " Dx11 = " << setw(10) << left << DX[i11]
          << "\n Dy00 = " << setw(10) << left << DY[i00]
          <<   " Dy01 = " << setw(10) << left << DY[i01]
          <<   " Dy10 = " << setw(10) << left << DY[i10]
          <<   " Dy11 = " << setw(10) << left << DY[i11]
          << '\n';
      }
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  char const *
  Akima2Dspline::type_name() const
  { return "Akima2D"; }

}
