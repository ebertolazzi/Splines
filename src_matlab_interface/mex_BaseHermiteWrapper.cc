/****************************************************************************\
  Copyright (c) Enrico Bertolazzi 2016
  All Rights Reserved.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the file license.txt for more details.
\****************************************************************************/

#include "Splines.hh"
#include "mex_utils.hh"

#define MEX_ERROR_MESSAGE \
"=====================================================================================\n" \
"BaseHermite:  Compute Hermite Bases\n" \
"\n" \
"USAGE:\n" \
"    [ h0, h1, dh0, dh1 ] = BaseHermite( 'base', t );\n" \
"    [ h0, h1, dh0, dh1 ] = BaseHermite( 'base_D', t );\n" \
"    [ h0, h1, dh0, dh1 ] = BaseHermite( 'base_DD', t );\n" \
"    [ h0, h1, dh0, dh1 ] = BaseHermite( 'base_DDD', t );\n" \
"    [ P0, P1, P2, P3 ]   = BaseHermite( 'hermite_to_bezier', P0, P1, T0, T1 );\n" \
"    [ P0, P1, T0, T1 ]   = BaseHermite( 'bezier_to_hermite', P0, P1, P2, P3 );\n" \
"    [ D1, sqrtD1 ]       = BaseHermite( 'L2_first_derivative' );\n" \
"    [ D2, sqrtD2 ]       = BaseHermite( 'L2_second_derivative' );\n" \
"    [ D3, sqrtD3 ]       = BaseHermite( 'L2_third_derivative' );\n" \
"    [ h0, h1, dh0, dh1, ddh0, ddh1 ] = BaseHermite( 'base5', t );\n" \
"    [ h0, h1, dh0, dh1, ddh0, ddh1 ] = BaseHermite( 'base5_D', t );\n" \
"    [ h0, h1, dh0, dh1, ddh0, ddh1 ] = BaseHermite( 'base5_DD', t );\n" \
"    [ h0, h1, dh0, dh1, ddh0, ddh1 ] = BaseHermite( 'base5_DDD', t );\n" \
"    [ h0, h1, dh0, dh1, ddh0, ddh1 ] = BaseHermite( 'base5_DDDD', t );\n" \
"    [ h0, h1, dh0, dh1, ddh0, ddh1 ] = BaseHermite( 'base5_DDDDD', t );\n" \
"\n" \

namespace Splines {

  using namespace std;

  /*\
   *                      _____                 _   _
   *  _ __ ___   _____  _|  ___|   _ _ __   ___| |_(_) ___  _ __
   * | '_ ` _ \ / _ \ \/ / |_ | | | | '_ \ / __| __| |/ _ \| '_ \
   * | | | | | |  __/>  <|  _|| |_| | | | | (__| |_| | (_) | | | |
   * |_| |_| |_|\___/_/\_\_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|
   *
  \*/

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_base( int nlhs, mxArray       *plhs[],
           int nrhs, mxArray const *prhs[] ) {
    #define CMD "BaseHermite('base',t [,H]): "
    MEX_ASSERT( nrhs == 2 || nrhs == 3, CMD "expected 2 or 3  inputs, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 1 || nlhs == 4, CMD "expected 1 or 4  outputs, nlhs = " << nlhs );
    mwSize nr;
    double H = 1;
    if ( nrhs == 3 ) H = getScalarValue( arg_in_2, CMD "argument H" );
    double const * t = getVectorPointer( arg_in_1, nr, CMD "argument t");
    if ( nlhs == 1 ) {
      double * out = createMatrixValue( arg_out_0, nr, 4);
      for ( mwSize i = 0; i < nr; ++i ) Hermite3( t[i], H, out ); out += 4;
    } else if ( nlhs == 4 ) {
      double * out1 = createMatrixValue( arg_out_0, nr, 1);
      double * out2 = createMatrixValue( arg_out_1, nr, 1);
      double * out3 = createMatrixValue( arg_out_2, nr, 1);
      double * out4 = createMatrixValue( arg_out_3, nr, 1);
      double tmp[4];
      for ( mwSize i = 0; i < nr; ++i ) {
        Hermite3( t[i], H, tmp );
        *out1++ = tmp[0];
        *out2++ = tmp[1];
        *out3++ = tmp[2];
        *out4++ = tmp[3];
      }
    }
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_base_D( int nlhs, mxArray       *plhs[],
             int nrhs, mxArray const *prhs[] ) {
    #define CMD "BaseHermite('base_D',t [,H]): "
    MEX_ASSERT( nrhs == 2 || nrhs == 3, CMD "expected 2 or 3  inputs, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 1 || nlhs == 4, CMD "expected 1 or 4  outputs, nlhs = " << nlhs );
    mwSize nr;
    double H = 1;
    if ( nrhs == 3 ) H = getScalarValue( arg_in_2, CMD "argument H" );
    double const * t = getVectorPointer( arg_in_1, nr, CMD "argument t");
    if ( nlhs == 1 ) {
      double * out = createMatrixValue( arg_out_0, nr, 4);
      for ( mwSize i = 0; i < nr; ++i ) Hermite3_D( t[i], H, out ); out += 4;
    } else if ( nlhs == 4 ) {
      double * out1 = createMatrixValue( arg_out_0, nr, 1);
      double * out2 = createMatrixValue( arg_out_1, nr, 1);
      double * out3 = createMatrixValue( arg_out_2, nr, 1);
      double * out4 = createMatrixValue( arg_out_3, nr, 1);
      double tmp[4];
      for ( mwSize i = 0; i < nr; ++i ) {
        Hermite3_D( t[i], H, tmp );
        *out1++ = tmp[0];
        *out2++ = tmp[1];
        *out3++ = tmp[2];
        *out4++ = tmp[3];
      }
    }
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_base_DD( int nlhs, mxArray       *plhs[],
              int nrhs, mxArray const *prhs[] ) {
    #define CMD "BaseHermite('base_DD',t [,H]): "
    MEX_ASSERT( nrhs == 2 || nrhs == 3, CMD "expected 2 or 3  inputs, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 1 || nlhs == 4, CMD "expected 1 or 4  outputs, nlhs = " << nlhs );
    mwSize nr;
    double H = 1;
    if ( nrhs == 3 ) H = getScalarValue( arg_in_2, CMD "argument H" );
    double const * t = getVectorPointer( arg_in_1, nr, CMD "argument t");
    if ( nlhs == 1 ) {
      double * out = createMatrixValue( arg_out_0, nr, 4);
      for ( mwSize i = 0; i < nr; ++i ) Hermite3_DD( t[i], H, out ); out += 4;
    } else if ( nlhs == 4 ) {
      double * out1 = createMatrixValue( arg_out_0, nr, 1);
      double * out2 = createMatrixValue( arg_out_1, nr, 1);
      double * out3 = createMatrixValue( arg_out_2, nr, 1);
      double * out4 = createMatrixValue( arg_out_3, nr, 1);
      double tmp[4];
      for ( mwSize i = 0; i < nr; ++i ) {
        Hermite3_DD( t[i], H, tmp );
        *out1++ = tmp[0];
        *out2++ = tmp[1];
        *out3++ = tmp[2];
        *out4++ = tmp[3];
      }
    }
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_base_DDD( int nlhs, mxArray       *plhs[],
               int nrhs, mxArray const *prhs[] ) {
    #define CMD "BaseHermite('base_DDD',t [,H]): "
    MEX_ASSERT( nrhs == 2 || nrhs == 3, CMD "expected 2 or 3  inputs, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 1 || nlhs == 4, CMD "expected 1 or 4  outputs, nlhs = " << nlhs );
    mwSize nr;
    double H = 1;
    if ( nrhs == 3 ) H = getScalarValue( arg_in_2, CMD "argument H" );
    double const * t = getVectorPointer( arg_in_1, nr, CMD "argument t");
    if ( nlhs == 1 ) {
      double * out = createMatrixValue( arg_out_0, nr, 4);
      for ( mwSize i = 0; i < nr; ++i ) Hermite3_DDD( t[i], H, out ); out += 4;
    } else if ( nlhs == 4 ) {
      double * out1 = createMatrixValue( arg_out_0, nr, 1);
      double * out2 = createMatrixValue( arg_out_1, nr, 1);
      double * out3 = createMatrixValue( arg_out_2, nr, 1);
      double * out4 = createMatrixValue( arg_out_3, nr, 1);
      double tmp[4];
      for ( mwSize i = 0; i < nr; ++i ) {
        Hermite3_DDD( t[i], H, tmp );
        *out1++ = tmp[0];
        *out2++ = tmp[1];
        *out3++ = tmp[2];
        *out4++ = tmp[3];
      }
    }
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_hermite_to_bezier( int nlhs, mxArray       *plhs[],
                        int nrhs, mxArray const *prhs[] ) {
    #define CMD "BaseHermite('hermite_to_bezier',P0,P1,T0,T1): "
    MEX_ASSERT( nrhs == 4, CMD "expected 4 inputs, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 4, CMD "expected 4 output, nlhs = " << nlhs );
    #undef CMD
    mwSize n1, n2, n3, n4;
    double const * P0 = getVectorPointer( arg_in_1, n1, CMD "P0" );
    double const * P1 = getVectorPointer( arg_in_2, n2, CMD "P1" );
    double const * T0 = getVectorPointer( arg_in_3, n3, CMD "T0" );
    double const * T1 = getVectorPointer( arg_in_4, n4, CMD "T1" );

    MEX_ASSERT(
      n1 == n2 && n2 == n3 && n3 == n4,
      CMD "bad dimensions |P0| = " << n1 << " |P1| = " << n2 <<
      " |T0| = " << n3 <<" |T1| = " << n4
    );
    double * oP0 = createMatrixValue( arg_out_0, n1, 1);
    double * oP1 = createMatrixValue( arg_out_1, n1, 1);
    double * oP3 = createMatrixValue( arg_out_2, n1, 1);
    double * oP4 = createMatrixValue( arg_out_3, n1, 1);
    for ( mwSize i = 0; i < n1; ++i ) {
      oP0[i] = P0[i];
      oP1[i] = P0[i] + T0[i]/3 ;
      oP1[i] = P1[i] - T1[i]/3 ;
      oP3[i] = P1[i] ;
    }
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_bezier_to_hermite( int nlhs, mxArray       *plhs[],
                        int nrhs, mxArray const *prhs[] ) {
    #define CMD "BaseHermite('bezier_to_hermite',P0,P1,P2,P3): "
    MEX_ASSERT( nrhs == 4, CMD "expected 4 inputs, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 4, CMD "expected 4 output, nlhs = " << nlhs );
    #undef CMD
    mwSize n1, n2, n3, n4;
    double const * P0 = getVectorPointer( arg_in_1, n1, CMD "P0" );
    double const * P1 = getVectorPointer( arg_in_2, n2, CMD "P1" );
    double const * P2 = getVectorPointer( arg_in_3, n3, CMD "P2" );
    double const * P3 = getVectorPointer( arg_in_4, n4, CMD "P3" );

    MEX_ASSERT(
      n1 == n2 && n2 == n3 && n3 == n4,
      CMD "bad dimensions |P0| = " << n1 << " |P1| = " << n2 <<
      " |P2| = " << n3 <<" |P3| = " << n4
    );
    double * oP0 = createMatrixValue( arg_out_0, n1, 1);
    double * oP1 = createMatrixValue( arg_out_1, n1, 1);
    double * oT0 = createMatrixValue( arg_out_2, n1, 1);
    double * oT1 = createMatrixValue( arg_out_3, n1, 1);
    for ( mwSize i = 0; i < n1; ++i ) {
      oP0[i] = P0[i];
      oP1[i] = P3[i];
      oT0[i] = 3*(P1[i] - P0[i]);
      oT1[i] = 3*(P3[i] - P2[i]);
    }
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_L2_first_derivative( int nlhs, mxArray       *plhs[],
                          int nrhs, mxArray const *prhs[] ) {
    #define CMD "BaseHermite('L2_first_derivative'): "
    MEX_ASSERT( nrhs == 1, CMD "expected 1 input, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 2, CMD "expected 2 outputs, nlhs = " << nlhs );
    #undef CMD
    double * D1     = createMatrixValue( arg_out_0, 4, 4 );
    double * sqrtD1 = createMatrixValue( arg_out_0, 3, 4 );
    *D1++ = 36/30.0;
    *D1++ = -36/30.0;
    *D1++ = 3/30.0;
    *D1++ = 3/30.0;

    *D1++ = 3/30.0;
    *D1++ = 36/30.0;
    *D1++ = -3/30.0;
    *D1++ = -3/30.0;

    *D1++ = 3/30.0;
    *D1++ = -3/30.0;
    *D1++ = 4/30.0;
    *D1++ = -1/30.0;

    *D1++ = 3/30.0;
    *D1++ = -3/30.0;
    *D1++ = -1/30.0;
    *D1++ = 4/30.0;

    *sqrtD1++ = sqrt(30.0)/5.0;
    *sqrtD1++ = 0;
    *sqrtD1++ = 0;

    *sqrtD1++ = -sqrt(30.0)/5.0;
    *sqrtD1++ = 0;
    *sqrtD1++ = 0;

    *sqrtD1++ = sqrt(30.0)/60.0;
    *sqrtD1++ = sqrt(2.0)/4.0;
    *sqrtD1++ = 0;

    *sqrtD1++ = sqrt(30.0)/60.0;
    *sqrtD1++ = -sqrt(2.0)/12.0;
    *sqrtD1++ = 1.0/3.0;

  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_L2_second_derivative( int nlhs, mxArray       *plhs[],
                           int nrhs, mxArray const *prhs[] ) {
    #define CMD "BaseHermite('L2_second_derivative','): "
    MEX_ASSERT( nrhs == 1, CMD "expected 1 input, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 2, CMD "expected 2 outputs, nlhs = " << nlhs );
    #undef CMD
    double * D2     = createMatrixValue( arg_out_0, 4, 4 );
    double * sqrtD2 = createMatrixValue( arg_out_0, 2, 4 );
    *D2++ = 12;
    *D2++ = -12;
    *D2++ = 6;
    *D2++ = 6;

    *D2++ = -12;
    *D2++ = 12;
    *D2++ = -6;
    *D2++ = -6;

    *D2++ = 6;
    *D2++ = -6;
    *D2++ = 4;
    *D2++ = 2;

    *D2++ = 6;
    *D2++ = -6;
    *D2++ = 2;
    *D2++ = 4;

    *sqrtD2++ = 2*sqrt(3.0);
    *sqrtD2++ = 0;

    *sqrtD2++ = -2*sqrt(3.0);
    *sqrtD2++ = 0;

    *sqrtD2++ = sqrt(3.0);
    *sqrtD2++ = 1;

    *sqrtD2++ = sqrt(3.0);
    *sqrtD2++ = -1;
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_L2_third_derivative( int nlhs, mxArray       *plhs[],
                          int nrhs, mxArray const *prhs[] ) {
    #define CMD "BaseHermite('L2_third_derivative','): "
    MEX_ASSERT( nrhs == 1, CMD "expected 1 input, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 2, CMD "expected 2 outputs, nlhs = " << nlhs );
    #undef CMD
    double * D3     = createMatrixValue( arg_out_0, 4, 4 );
    double * sqrtD3 = createMatrixValue( arg_out_0, 1, 4 );
    *D3++ = 144;
    *D3++ = -144;
    *D3++ = 72;
    *D3++ = 72;

    *D3++ = -144;
    *D3++ = 144;
    *D3++ = -72;
    *D3++ = -72;

    *D3++ = 72;
    *D3++ = -72;
    *D3++ = 36;
    *D3++ = 36;

    *D3++ = 72;
    *D3++ = -72;
    *D3++ = 36;
    *D3++ = 36;

    *sqrtD3++ = 12;
    *sqrtD3++ = -12;
    *sqrtD3++ = 6;
    *sqrtD3++ = 6;
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_base5( int nlhs, mxArray       *plhs[],
            int nrhs, mxArray const *prhs[] ) {
    #define CMD "BaseHermite('base5',t [,H]): "
    MEX_ASSERT( nrhs == 2 || nrhs == 3, CMD "expected 2 or 3  inputs, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 1 || nlhs == 4, CMD "expected 1 or 6  outputs, nlhs = " << nlhs );
    mwSize nr;
    double H = 1;
    if ( nrhs == 3 ) H = getScalarValue( arg_in_2, CMD "argument H" );
    double const * t = getVectorPointer( arg_in_1, nr, CMD "argument t");
    if ( nlhs == 1 ) {
      double * out = createMatrixValue( arg_out_0, nr, 6);
      for ( mwSize i = 0; i < nr; ++i ) Hermite5( t[i], H, out ); out += 6;
    } else if ( nlhs == 4 ) {
      double * out1 = createMatrixValue( arg_out_0, nr, 1);
      double * out2 = createMatrixValue( arg_out_1, nr, 1);
      double * out3 = createMatrixValue( arg_out_2, nr, 1);
      double * out4 = createMatrixValue( arg_out_3, nr, 1);
      double * out5 = createMatrixValue( arg_out_4, nr, 1);
      double * out6 = createMatrixValue( arg_out_5, nr, 1);
      double tmp[6];
      for ( mwSize i = 0; i < nr; ++i ) {
        Hermite5( t[i], H, tmp );
        *out1++ = tmp[0];
        *out2++ = tmp[1];
        *out3++ = tmp[2];
        *out4++ = tmp[3];
        *out5++ = tmp[4];
        *out6++ = tmp[5];
      }
    }
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_base5_D( int nlhs, mxArray       *plhs[],
              int nrhs, mxArray const *prhs[] ) {
    #define CMD "BaseHermite('base5_D',t [,H]): "
    MEX_ASSERT( nrhs == 2 || nrhs == 3, CMD "expected 2 or 3  inputs, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 1 || nlhs == 4, CMD "expected 1 or 6  outputs, nlhs = " << nlhs );
    mwSize nr;
    double H = 1;
    if ( nrhs == 3 ) H = getScalarValue( arg_in_2, CMD "argument H" );
    double const * t = getVectorPointer( arg_in_1, nr, CMD "argument t");
    if ( nlhs == 1 ) {
      double * out = createMatrixValue( arg_out_0, nr, 6);
      for ( mwSize i = 0; i < nr; ++i ) Hermite5_D( t[i], H, out ); out += 6;
    } else if ( nlhs == 4 ) {
      double * out1 = createMatrixValue( arg_out_0, nr, 1);
      double * out2 = createMatrixValue( arg_out_1, nr, 1);
      double * out3 = createMatrixValue( arg_out_2, nr, 1);
      double * out4 = createMatrixValue( arg_out_3, nr, 1);
      double * out5 = createMatrixValue( arg_out_4, nr, 1);
      double * out6 = createMatrixValue( arg_out_5, nr, 1);
      double tmp[6];
      for ( mwSize i = 0; i < nr; ++i ) {
        Hermite5_D( t[i], H, tmp );
        *out1++ = tmp[0];
        *out2++ = tmp[1];
        *out3++ = tmp[2];
        *out4++ = tmp[3];
        *out5++ = tmp[4];
        *out6++ = tmp[5];
      }
    }
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_base5_DD( int nlhs, mxArray       *plhs[],
               int nrhs, mxArray const *prhs[] ) {
    #define CMD "BaseHermite('base5_DD',t [,H]): "
    MEX_ASSERT( nrhs == 2 || nrhs == 3, CMD "expected 2 or 3  inputs, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 1 || nlhs == 4, CMD "expected 1 or 6  outputs, nlhs = " << nlhs );
    mwSize nr;
    double H = 1;
    if ( nrhs == 3 ) H = getScalarValue( arg_in_2, CMD "argument H" );
    double const * t = getVectorPointer( arg_in_1, nr, CMD "argument t");
    if ( nlhs == 1 ) {
      double * out = createMatrixValue( arg_out_0, nr, 6);
      for ( mwSize i = 0; i < nr; ++i ) Hermite5_DD( t[i], H, out ); out += 6;
    } else if ( nlhs == 4 ) {
      double * out1 = createMatrixValue( arg_out_0, nr, 1);
      double * out2 = createMatrixValue( arg_out_1, nr, 1);
      double * out3 = createMatrixValue( arg_out_2, nr, 1);
      double * out4 = createMatrixValue( arg_out_3, nr, 1);
      double * out5 = createMatrixValue( arg_out_4, nr, 1);
      double * out6 = createMatrixValue( arg_out_5, nr, 1);
      double tmp[6];
      for ( mwSize i = 0; i < nr; ++i ) {
        Hermite5_DD( t[i], H, tmp );
        *out1++ = tmp[0];
        *out2++ = tmp[1];
        *out3++ = tmp[2];
        *out4++ = tmp[3];
        *out5++ = tmp[4];
        *out6++ = tmp[5];
      }
    }
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_base5_DDD( int nlhs, mxArray       *plhs[],
                int nrhs, mxArray const *prhs[] ) {
    #define CMD "BaseHermite('base5_DDD',t [,H]): "
    MEX_ASSERT( nrhs == 2 || nrhs == 3, CMD "expected 2 or 3  inputs, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 1 || nlhs == 4, CMD "expected 1 or 6  outputs, nlhs = " << nlhs );
    mwSize nr;
    double H = 1;
    if ( nrhs == 3 ) H = getScalarValue( arg_in_2, CMD "argument H" );
    double const * t = getVectorPointer( arg_in_1, nr, CMD "argument t");
    if ( nlhs == 1 ) {
      double * out = createMatrixValue( arg_out_0, nr, 6);
      for ( mwSize i = 0; i < nr; ++i ) Hermite5_DDD( t[i], H, out ); out += 6;
    } else if ( nlhs == 4 ) {
      double * out1 = createMatrixValue( arg_out_0, nr, 1);
      double * out2 = createMatrixValue( arg_out_1, nr, 1);
      double * out3 = createMatrixValue( arg_out_2, nr, 1);
      double * out4 = createMatrixValue( arg_out_3, nr, 1);
      double * out5 = createMatrixValue( arg_out_4, nr, 1);
      double * out6 = createMatrixValue( arg_out_5, nr, 1);
      double tmp[6];
      for ( mwSize i = 0; i < nr; ++i ) {
        Hermite5_DDD( t[i], H, tmp );
        *out1++ = tmp[0];
        *out2++ = tmp[1];
        *out3++ = tmp[2];
        *out4++ = tmp[3];
        *out5++ = tmp[4];
        *out6++ = tmp[5];
      }
    }
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_base5_DDDD( int nlhs, mxArray       *plhs[],
                 int nrhs, mxArray const *prhs[] ) {
    #define CMD "BaseHermite('base5_DDDD',t [,H]): "
    MEX_ASSERT( nrhs == 2 || nrhs == 3, CMD "expected 2 or 3  inputs, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 1 || nlhs == 4, CMD "expected 1 or 6  outputs, nlhs = " << nlhs );
    mwSize nr;
    double H = 1;
    if ( nrhs == 3 ) H = getScalarValue( arg_in_2, CMD "argument H" );
    double const * t = getVectorPointer( arg_in_1, nr, CMD "argument t");
    if ( nlhs == 1 ) {
      double * out = createMatrixValue( arg_out_0, nr, 6);
      for ( mwSize i = 0; i < nr; ++i ) Hermite5_DDDD( t[i], H, out ); out += 6;
    } else if ( nlhs == 4 ) {
      double * out1 = createMatrixValue( arg_out_0, nr, 1);
      double * out2 = createMatrixValue( arg_out_1, nr, 1);
      double * out3 = createMatrixValue( arg_out_2, nr, 1);
      double * out4 = createMatrixValue( arg_out_3, nr, 1);
      double * out5 = createMatrixValue( arg_out_4, nr, 1);
      double * out6 = createMatrixValue( arg_out_5, nr, 1);
      double tmp[6];
      for ( mwSize i = 0; i < nr; ++i ) {
        Hermite5_DDDD( t[i], H, tmp );
        *out1++ = tmp[0];
        *out2++ = tmp[1];
        *out3++ = tmp[2];
        *out4++ = tmp[3];
        *out5++ = tmp[4];
        *out6++ = tmp[5];
      }
    }
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_base5_DDDDD( int nlhs, mxArray       *plhs[],
                  int nrhs, mxArray const *prhs[] ) {
    #define CMD "BaseHermite('base5_DDDDD',t [,H]): "
    MEX_ASSERT( nrhs == 2 || nrhs == 3, CMD "expected 2 or 3  inputs, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 1 || nlhs == 4, CMD "expected 1 or 6  outputs, nlhs = " << nlhs );
    mwSize nr;
    double H = 1;
    if ( nrhs == 3 ) H = getScalarValue( arg_in_2, CMD "argument H" );
    double const * t = getVectorPointer( arg_in_1, nr, CMD "argument t");
    if ( nlhs == 1 ) {
      double * out = createMatrixValue( arg_out_0, nr, 6);
      for ( mwSize i = 0; i < nr; ++i ) Hermite5_DDDDD( t[i], H, out ); out += 6;
    } else if ( nlhs == 4 ) {
      double * out1 = createMatrixValue( arg_out_0, nr, 1);
      double * out2 = createMatrixValue( arg_out_1, nr, 1);
      double * out3 = createMatrixValue( arg_out_2, nr, 1);
      double * out4 = createMatrixValue( arg_out_3, nr, 1);
      double * out5 = createMatrixValue( arg_out_4, nr, 1);
      double * out6 = createMatrixValue( arg_out_5, nr, 1);
      double tmp[6];
      for ( mwSize i = 0; i < nr; ++i ) {
        Hermite5_DDDDD( t[i], H, tmp );
        *out1++ = tmp[0];
        *out2++ = tmp[1];
        *out3++ = tmp[2];
        *out4++ = tmp[3];
        *out5++ = tmp[4];
        *out6++ = tmp[5];
      }
    }
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  typedef enum {
    CMD_BASE,
    CMD_BASE_D,
    CMD_BASE_DD,
    CMD_BASE_DDD,
    CMD_HERMITE_TO_BEZIER,
    CMD_BEZIER_TO_HERMITE,
    CMD_L2_1,
    CMD_L2_2,
    CMD_L2_3,
    CMD_BASE5,
    CMD_BASE5_D,
    CMD_BASE5_DD,
    CMD_BASE5_DDD,
    CMD_BASE5_DDDD,
    CMD_BASE5_DDDDD
  } CMD_LIST;

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static map<string,unsigned> cmd_to_idx = {
    {"base",CMD_BASE},
    {"base_D",CMD_BASE_D},
    {"base_DD",CMD_BASE_DD},
    {"base_DDD",CMD_BASE_DDD},
    {"hermite_to_bezier",CMD_HERMITE_TO_BEZIER},
    {"bezier_to_hermite",CMD_BEZIER_TO_HERMITE},
    {"L2_first_derivative",CMD_L2_1},
    {"L2_second_derivative",CMD_L2_2},
    {"L2_third_derivative",CMD_L2_3},
    {"base5",CMD_BASE5},
    {"base5_D",CMD_BASE5_D},
    {"base5_DD",CMD_BASE5_DD},
    {"base5_DDD",CMD_BASE5_DDD},
    {"base5_DDDD",CMD_BASE5_DDDD},
    {"base5_DDDDD",CMD_BASE5_DDDDD},
    CMD_MAP_LIST
  };

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  extern "C"
  void
  mexFunction( int nlhs, mxArray       *plhs[],
               int nrhs, mxArray const *prhs[] ) {
    // the first argument must be a string
    if ( nrhs == 0 ) {
      mexErrMsgTxt(MEX_ERROR_MESSAGE);
      return;
    }

    try {

      MEX_ASSERT( mxIsChar(arg_in_0), "First argument must be a string" );
      string cmd = mxArrayToString(arg_in_0);

      switch ( cmd_to_idx.at(cmd) ) {
      case CMD_BASE:
        do_base( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_BASE_D:
        do_base_D( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_BASE_DD:
        do_base_DD( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_BASE_DDD:
        do_base_DDD( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_HERMITE_TO_BEZIER:
        do_hermite_to_bezier( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_BEZIER_TO_HERMITE:
        do_bezier_to_hermite( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_L2_1:
        do_L2_first_derivative( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_L2_2:
        do_L2_second_derivative( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_L2_3:
        do_L2_third_derivative( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_BASE5:
        do_base5( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_BASE5:
        do_base5( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_BASE5_D:
        do_base5_D( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_BASE5_DD:
        do_base5_DD( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_BASE5_DDD:
        do_base5_DDD( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_BASE5_DDDD:
        do_base5_DDDD( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_BASE5_DDDDD:
        do_base5_DDDDD( nlhs, plhs, nrhs, prhs );
        break;
      }

    } catch ( exception const & e ) {
      mexErrMsgTxt(e.what());
    } catch (...) {
      mexErrMsgTxt("BaseHermite failed\n");
    }

  }
}
