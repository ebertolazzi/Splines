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
"    [ h0, h1, dh0, dh1 ] = BaseHermite( 'base', t [,H] );\n" \
"    [ h0, h1, dh0, dh1 ] = BaseHermite( 'base_D', t [,H]);\n" \
"    [ h0, h1, dh0, dh1 ] = BaseHermite( 'base_DD', t [,H]);\n" \
"    [ h0, h1, dh0, dh1 ] = BaseHermite( 'base_DDD', t [,H]);\n" \
"    P    = BaseHermite( 'eval', t, P0, P1, T0, T1 [,H] );\n" \
"    dP   = BaseHermite( 'eval_D', t, P0, P1, T0, T1 [,H] );\n" \
"    ddP  = BaseHermite( 'eval_DD', t, P0, P1, T0, T1 [,H] );\n" \
"    dddP = BaseHermite( 'eval_DDD', t, P0, P1, T0, T1 [,H] );\n" \
"    [ h0, h1, dh0, dh1, ddh0, ddh1 ] = BaseHermite( 'base5', t );\n" \
"    [ h0, h1, dh0, dh1, ddh0, ddh1 ] = BaseHermite( 'base5_D', t );\n" \
"    [ h0, h1, dh0, dh1, ddh0, ddh1 ] = BaseHermite( 'base5_DD', t );\n" \
"    [ h0, h1, dh0, dh1, ddh0, ddh1 ] = BaseHermite( 'base5_DDD', t );\n" \
"    [ h0, h1, dh0, dh1, ddh0, ddh1 ] = BaseHermite( 'base5_DDDD', t );\n" \
"    [ h0, h1, dh0, dh1, ddh0, ddh1 ] = BaseHermite( 'base5_DDDDD', t );\n" \
"    P      = BaseHermite( 'eval5', t, P0, P1, T0, T1, J0, J1 [,H] );\n" \
"    dP     = BaseHermite( 'eval5_D', t, P0, P1, T0, T1, J0, J1 [,H] );\n" \
"    ddP    = BaseHermite( 'eval5_DD', t, P0, P1, T0, T1, J0, J1 [,H] );\n" \
"    dddP   = BaseHermite( 'eval5_DDD', t, P0, P1, T0, T1, J0, J1 [,H] );\n" \
"    ddddP  = BaseHermite( 'eval5_DDDD', t, P0, P1, T0, T1, J0, J1 [,H] );\n" \
"    dddddP = BaseHermite( 'eval5_DDDDD', t, P0, P1, T0, T1, J0, J1 [,H] );\n" \
"    [ P0, P1, P2, P3 ]   = BaseHermite( 'hermite_to_bezier', P0, P1, T0, T1 );\n" \
"    [ P0, P1, T0, T1 ]   = BaseHermite( 'bezier_to_hermite', P0, P1, P2, P3 );\n" \
"    sqrtD1 = BaseHermite( 'L2_first_derivative' );\n" \
"    sqrtD2 = BaseHermite( 'L2_second_derivative' );\n" \
"    sqrtD3 = BaseHermite( 'L2_third_derivative' );\n" \
"    [ P0, P1, T0, T1 ]   = BaseHermite( 'cut', a, b, P0, P1, T0, T1, [,H] );\n" \
"    len                  = BaseHermite( 'approximate_length', P0, P1, T0, T1, [,H] );\n" \
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
    double *out1, *out2, *out3, *out4;
    if ( nlhs == 1 ) {
      out1 = createMatrixValue( arg_out_0, nr, 4);
      out2 = out1 + nr;
      out3 = out2 + nr;
      out4 = out3 + nr;
    } else if ( nlhs == 4 ) {
      out1 = createMatrixValue( arg_out_0, nr, 1);
      out2 = createMatrixValue( arg_out_1, nr, 1);
      out3 = createMatrixValue( arg_out_2, nr, 1);
      out4 = createMatrixValue( arg_out_3, nr, 1);
    }
    double tmp[4];
    for ( mwSize i = 0; i < nr; ++i ) {
      Hermite3( t[i], H, tmp );
      *out1++ = tmp[0];
      *out2++ = tmp[1];
      *out3++ = tmp[2];
      *out4++ = tmp[3];
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
    double *out1, *out2, *out3, *out4;
    if ( nlhs == 1 ) {
      out1 = createMatrixValue( arg_out_0, nr, 4);
      out2 = out1 + nr;
      out3 = out2 + nr;
      out4 = out3 + nr;
    } else if ( nlhs == 4 ) {
      out1 = createMatrixValue( arg_out_0, nr, 1);
      out2 = createMatrixValue( arg_out_1, nr, 1);
      out3 = createMatrixValue( arg_out_2, nr, 1);
      out4 = createMatrixValue( arg_out_3, nr, 1);
    }
    double tmp[4];
    for ( mwSize i = 0; i < nr; ++i ) {
      Hermite3_D( t[i], H, tmp );
      *out1++ = tmp[0];
      *out2++ = tmp[1];
      *out3++ = tmp[2];
      *out4++ = tmp[3];
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
    double *out1, *out2, *out3, *out4;
    if ( nlhs == 1 ) {
      out1 = createMatrixValue( arg_out_0, nr, 4);
      out2 = out1 + nr;
      out3 = out2 + nr;
      out4 = out3 + nr;
    } else if ( nlhs == 4 ) {
      out1 = createMatrixValue( arg_out_0, nr, 1);
      out2 = createMatrixValue( arg_out_1, nr, 1);
      out3 = createMatrixValue( arg_out_2, nr, 1);
      out4 = createMatrixValue( arg_out_3, nr, 1);
    }
    double tmp[4];
    for ( mwSize i = 0; i < nr; ++i ) {
      Hermite3_DD( t[i], H, tmp );
      *out1++ = tmp[0];
      *out2++ = tmp[1];
      *out3++ = tmp[2];
      *out4++ = tmp[3];
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
    double *out1, *out2, *out3, *out4;
    if ( nlhs == 1 ) {
      out1 = createMatrixValue( arg_out_0, nr, 4);
      out2 = out1 + nr;
      out3 = out2 + nr;
      out4 = out3 + nr;
    } else if ( nlhs == 4 ) {
      out1 = createMatrixValue( arg_out_0, nr, 1);
      out2 = createMatrixValue( arg_out_1, nr, 1);
      out3 = createMatrixValue( arg_out_2, nr, 1);
      out4 = createMatrixValue( arg_out_3, nr, 1);
    }
    double tmp[4];
    for ( mwSize i = 0; i < nr; ++i ) {
      Hermite3_DDD( t[i], H, tmp );
      *out1++ = tmp[0];
      *out2++ = tmp[1];
      *out3++ = tmp[2];
      *out4++ = tmp[3];
    }
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_eval( int nlhs, mxArray       *plhs[],
           int nrhs, mxArray const *prhs[] ) {
    #define CMD "BaseHermite( 'eval', t, P0, P1, T0, T1 [,H] ): "
    MEX_ASSERT( nrhs == 6 || nrhs == 7, CMD "expected 6 or 7 inputs, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
    mwSize nr, nr2, nr3, nr4, nr5;
    double H = 1;
    double const * t  = getVectorPointer( arg_in_1, nr,  CMD "argument t");
    double const * P0 = getVectorPointer( arg_in_2, nr2, CMD "argument t");
    double const * P1 = getVectorPointer( arg_in_3, nr3, CMD "argument t");
    double const * T0 = getVectorPointer( arg_in_4, nr4, CMD "argument t");
    double const * T1 = getVectorPointer( arg_in_5, nr5, CMD "argument t");
    MEX_ASSERT(
      nr2 == nr3 && nr3 == nr4 && nr4 == nr5,
      CMD " bad dimension |P0| = " << nr2 << " |P1| = " << nr3 <<
      " |T0| = " << nr4 << " |T1| = " << nr5
    );
    if ( nrhs == 7 ) H = getScalarValue( arg_in_6, CMD "argument H" );
    double * out = createMatrixValue( arg_out_0, nr2, nr );
    double tmp[4];
    for ( mwSize i = 0; i < nr; ++i ) {
      Hermite3( t[i], H, tmp );
      for ( mwSize j = 0; j < nr2; ++j )
        *out++ = tmp[0] * P0[j] + tmp[1] * P1[j] + tmp[2] * T0[j] + tmp[3] * T1[j] ;
    }
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_eval_D( int nlhs, mxArray       *plhs[],
             int nrhs, mxArray const *prhs[] ) {
    #define CMD "BaseHermite( 'eval_D', t, P0, P1, T0, T1 [,H] ): "
    MEX_ASSERT( nrhs == 6 || nrhs == 7, CMD "expected 6 or 7 inputs, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
    mwSize nr, nr2, nr3, nr4, nr5;
    double H = 1;
    double const * t  = getVectorPointer( arg_in_1, nr,  CMD "argument t");
    double const * P0 = getVectorPointer( arg_in_2, nr2, CMD "argument t");
    double const * P1 = getVectorPointer( arg_in_3, nr3, CMD "argument t");
    double const * T0 = getVectorPointer( arg_in_4, nr4, CMD "argument t");
    double const * T1 = getVectorPointer( arg_in_5, nr5, CMD "argument t");
    MEX_ASSERT(
      nr2 == nr3 && nr3 == nr4 && nr4 == nr5,
      CMD " bad dimension |P0| = " << nr2 << " |P1| = " << nr3 <<
      " |T0| = " << nr4 << " |T1| = " << nr5
    );
    if ( nrhs == 7 ) H = getScalarValue( arg_in_6, CMD "argument H" );
    double * out = createMatrixValue( arg_out_0, nr2, nr );
    double tmp[4];
    for ( mwSize i = 0; i < nr; ++i ) {
      Hermite3_D( t[i], H, tmp );
      for ( mwSize j = 0; j < nr2; ++j )
        *out++ = tmp[0] * P0[j] + tmp[1] * P1[j] + tmp[2] * T0[j] + tmp[3] * T1[j] ;
    }
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_eval_DD( int nlhs, mxArray       *plhs[],
              int nrhs, mxArray const *prhs[] ) {
    #define CMD "BaseHermite( 'eval_DD', t, P0, P1, T0, T1 [,H] ): "
    MEX_ASSERT( nrhs == 6 || nrhs == 7, CMD "expected 6 or 7 inputs, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
    mwSize nr, nr2, nr3, nr4, nr5;
    double H = 1;
    double const * t  = getVectorPointer( arg_in_1, nr,  CMD "argument t");
    double const * P0 = getVectorPointer( arg_in_2, nr2, CMD "argument t");
    double const * P1 = getVectorPointer( arg_in_3, nr3, CMD "argument t");
    double const * T0 = getVectorPointer( arg_in_4, nr4, CMD "argument t");
    double const * T1 = getVectorPointer( arg_in_5, nr5, CMD "argument t");
    MEX_ASSERT(
      nr2 == nr3 && nr3 == nr4 && nr4 == nr5,
      CMD " bad dimension |P0| = " << nr2 << " |P1| = " << nr3 <<
      " |T0| = " << nr4 << " |T1| = " << nr5
    );
    if ( nrhs == 7 ) H = getScalarValue( arg_in_6, CMD "argument H" );
    double * out = createMatrixValue( arg_out_0, nr2, nr );
    double tmp[4];
    for ( mwSize i = 0; i < nr; ++i ) {
      Hermite3_DD( t[i], H, tmp );
      for ( mwSize j = 0; j < nr2; ++j )
        *out++ = tmp[0] * P0[j] + tmp[1] * P1[j] + tmp[2] * T0[j] + tmp[3] * T1[j] ;
    }
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_eval_DDD( int nlhs, mxArray       *plhs[],
               int nrhs, mxArray const *prhs[] ) {
    #define CMD "BaseHermite( 'eval_DDD', t, P0, P1, T0, T1 [,H] ): "
    MEX_ASSERT( nrhs == 6 || nrhs == 7, CMD "expected 6 or 7 inputs, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
    mwSize nr, nr2, nr3, nr4, nr5;
    double H = 1;
    double const * t  = getVectorPointer( arg_in_1, nr,  CMD "argument t");
    double const * P0 = getVectorPointer( arg_in_2, nr2, CMD "argument t");
    double const * P1 = getVectorPointer( arg_in_3, nr3, CMD "argument t");
    double const * T0 = getVectorPointer( arg_in_4, nr4, CMD "argument t");
    double const * T1 = getVectorPointer( arg_in_5, nr5, CMD "argument t");
    MEX_ASSERT(
      nr2 == nr3 && nr3 == nr4 && nr4 == nr5,
      CMD " bad dimension |P0| = " << nr2 << " |P1| = " << nr3 <<
      " |T0| = " << nr4 << " |T1| = " << nr5
    );
    if ( nrhs == 7 ) H = getScalarValue( arg_in_6, CMD "argument H" );
    double * out = createMatrixValue( arg_out_0, nr2, nr );
    double tmp[4];
    for ( mwSize i = 0; i < nr; ++i ) {
      Hermite3_DDD( t[i], H, tmp );
      for ( mwSize j = 0; j < nr2; ++j )
        *out++ = tmp[0] * P0[j] + tmp[1] * P1[j] + tmp[2] * T0[j] + tmp[3] * T1[j] ;
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
    double * oP2 = createMatrixValue( arg_out_2, n1, 1);
    double * oP3 = createMatrixValue( arg_out_3, n1, 1);
    for ( mwSize i = 0; i < n1; ++i ) {
      oP0[i] = P0[i];
      oP1[i] = P0[i] + T0[i]/3 ;
      oP2[i] = P1[i] - T1[i]/3 ;
      oP3[i] = P1[i] ;
    }
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_bezier_to_hermite( int nlhs, mxArray       *plhs[],
                        int nrhs, mxArray const *prhs[] ) {
    #define CMD "BaseHermite('bezier_to_hermite',P0,P1,P2,P3): "
    MEX_ASSERT( nrhs == 4, CMD "expected 4 inputs, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 4, CMD "expected 4 output, nlhs = " << nlhs );
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
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_L2_first_derivative( int nlhs, mxArray       *plhs[],
                          int nrhs, mxArray const *prhs[] ) {
    #define CMD "BaseHermite('L2_first_derivative'): "
    MEX_ASSERT( nrhs == 1, CMD "expected 1 input, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
    #undef CMD
    double * sqrtD1 = createMatrixValue( arg_out_0, 3, 4 );

    double t1 = sqrt(5.0)/5.0;
    double t2 = sqrt(3.0)/6.0;

    *sqrtD1++ =-1;
    *sqrtD1++ = 0;
    *sqrtD1++ = t1;

    *sqrtD1++ = 1;
    *sqrtD1++ = 0;
    *sqrtD1++ = -t1;

    *sqrtD1++ = 0;
    *sqrtD1++ = -t2;
    *sqrtD1++ = t1/2;

    *sqrtD1++ = 0;
    *sqrtD1++ = t2;
    *sqrtD1++ = t1/2;

  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_L2_second_derivative( int nlhs, mxArray       *plhs[],
                           int nrhs, mxArray const *prhs[] ) {
    #define CMD "BaseHermite('L2_second_derivative','): "
    MEX_ASSERT( nrhs == 1, CMD "expected 1 input, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
    #undef CMD
    double * sqrtD2 = createMatrixValue( arg_out_0, 2, 4 );

    double t1 = sqrt(3.0);

    *sqrtD2++ = 0;
    *sqrtD2++ = 2*t1;

    *sqrtD2++ = 0;
    *sqrtD2++ = -2*t1;

    *sqrtD2++ = -1;
    *sqrtD2++ = t1;

    *sqrtD2++ = 1;
    *sqrtD2++ = t1;
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_L2_third_derivative( int nlhs, mxArray       *plhs[],
                          int nrhs, mxArray const *prhs[] ) {
    #define CMD "BaseHermite('L2_third_derivative','): "
    MEX_ASSERT( nrhs == 1, CMD "expected 1 input, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
    #undef CMD
    double * sqrtD3 = createMatrixValue( arg_out_0, 1, 4 );

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
    MEX_ASSERT( nlhs == 1 || nlhs == 6, CMD "expected 1 or 6  outputs, nlhs = " << nlhs );
    mwSize nr;
    double H = 1;
    if ( nrhs == 3 ) H = getScalarValue( arg_in_2, CMD "argument H" );
    double const * t = getVectorPointer( arg_in_1, nr, CMD "argument t");
    double *out1, *out2, *out3, *out4, *out5, *out6;
    if ( nlhs == 1 ) {
      out1 = createMatrixValue( arg_out_0, nr, 6);
      out2 = out1 + nr;
      out3 = out2 + nr;
      out4 = out3 + nr;
      out5 = out4 + nr;
      out6 = out5 + nr;
    } else if ( nlhs == 6 ) {
      out1 = createMatrixValue( arg_out_0, nr, 1);
      out2 = createMatrixValue( arg_out_1, nr, 1);
      out3 = createMatrixValue( arg_out_2, nr, 1);
      out4 = createMatrixValue( arg_out_3, nr, 1);
      out5 = createMatrixValue( arg_out_4, nr, 1);
      out6 = createMatrixValue( arg_out_5, nr, 1);
    }
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
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_base5_D( int nlhs, mxArray       *plhs[],
              int nrhs, mxArray const *prhs[] ) {
    #define CMD "BaseHermite('base5_D',t [,H]): "
    MEX_ASSERT( nrhs == 2 || nrhs == 3, CMD "expected 2 or 3  inputs, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 1 || nlhs == 6, CMD "expected 1 or 6  outputs, nlhs = " << nlhs );
    mwSize nr;
    double H = 1;
    if ( nrhs == 3 ) H = getScalarValue( arg_in_2, CMD "argument H" );
    double const * t = getVectorPointer( arg_in_1, nr, CMD "argument t");
    double *out1, *out2, *out3, *out4, *out5, *out6;
    if ( nlhs == 1 ) {
      out1 = createMatrixValue( arg_out_0, nr, 6);
      out2 = out1 + nr;
      out3 = out2 + nr;
      out4 = out3 + nr;
      out5 = out4 + nr;
      out6 = out5 + nr;
    } else if ( nlhs == 6 ) {
      out1 = createMatrixValue( arg_out_0, nr, 1);
      out2 = createMatrixValue( arg_out_1, nr, 1);
      out3 = createMatrixValue( arg_out_2, nr, 1);
      out4 = createMatrixValue( arg_out_3, nr, 1);
      out5 = createMatrixValue( arg_out_4, nr, 1);
      out6 = createMatrixValue( arg_out_5, nr, 1);
    }
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
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_base5_DD( int nlhs, mxArray       *plhs[],
               int nrhs, mxArray const *prhs[] ) {
    #define CMD "BaseHermite('base5_DD',t [,H]): "
    MEX_ASSERT( nrhs == 2 || nrhs == 3, CMD "expected 2 or 3  inputs, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 1 || nlhs == 6, CMD "expected 1 or 6  outputs, nlhs = " << nlhs );
    mwSize nr;
    double H = 1;
    if ( nrhs == 3 ) H = getScalarValue( arg_in_2, CMD "argument H" );
    double const * t = getVectorPointer( arg_in_1, nr, CMD "argument t");
    double *out1, *out2, *out3, *out4, *out5, *out6;
    if ( nlhs == 1 ) {
      out1 = createMatrixValue( arg_out_0, nr, 6);
      out2 = out1 + nr;
      out3 = out2 + nr;
      out4 = out3 + nr;
      out5 = out4 + nr;
      out6 = out5 + nr;
    } else if ( nlhs == 6 ) {
      out1 = createMatrixValue( arg_out_0, nr, 1);
      out2 = createMatrixValue( arg_out_1, nr, 1);
      out3 = createMatrixValue( arg_out_2, nr, 1);
      out4 = createMatrixValue( arg_out_3, nr, 1);
      out5 = createMatrixValue( arg_out_4, nr, 1);
      out6 = createMatrixValue( arg_out_5, nr, 1);
    }
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
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_base5_DDD( int nlhs, mxArray       *plhs[],
                int nrhs, mxArray const *prhs[] ) {
    #define CMD "BaseHermite('base5_DDD',t [,H]): "
    MEX_ASSERT( nrhs == 2 || nrhs == 3, CMD "expected 2 or 3  inputs, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 1 || nlhs == 6, CMD "expected 1 or 6  outputs, nlhs = " << nlhs );
    mwSize nr;
    double H = 1;
    if ( nrhs == 3 ) H = getScalarValue( arg_in_2, CMD "argument H" );
    double const * t = getVectorPointer( arg_in_1, nr, CMD "argument t");
    double *out1, *out2, *out3, *out4, *out5, *out6;
    if ( nlhs == 1 ) {
      out1 = createMatrixValue( arg_out_0, nr, 6);
      out2 = out1 + nr;
      out3 = out2 + nr;
      out4 = out3 + nr;
      out5 = out4 + nr;
      out6 = out5 + nr;
    } else if ( nlhs == 6 ) {
      out1 = createMatrixValue( arg_out_0, nr, 1);
      out2 = createMatrixValue( arg_out_1, nr, 1);
      out3 = createMatrixValue( arg_out_2, nr, 1);
      out4 = createMatrixValue( arg_out_3, nr, 1);
      out5 = createMatrixValue( arg_out_4, nr, 1);
      out6 = createMatrixValue( arg_out_5, nr, 1);
    }
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
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_base5_DDDD( int nlhs, mxArray       *plhs[],
                 int nrhs, mxArray const *prhs[] ) {
    #define CMD "BaseHermite_DDDD('base5',t [,H]): "
    MEX_ASSERT( nrhs == 2 || nrhs == 3, CMD "expected 2 or 3  inputs, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 1 || nlhs == 6, CMD "expected 1 or 6  outputs, nlhs = " << nlhs );
    mwSize nr;
    double H = 1;
    if ( nrhs == 3 ) H = getScalarValue( arg_in_2, CMD "argument H" );
    double const * t = getVectorPointer( arg_in_1, nr, CMD "argument t");
    double *out1, *out2, *out3, *out4, *out5, *out6;
    if ( nlhs == 1 ) {
      out1 = createMatrixValue( arg_out_0, nr, 6);
      out2 = out1 + nr;
      out3 = out2 + nr;
      out4 = out3 + nr;
      out5 = out4 + nr;
      out6 = out5 + nr;
    } else if ( nlhs == 6 ) {
      out1 = createMatrixValue( arg_out_0, nr, 1);
      out2 = createMatrixValue( arg_out_1, nr, 1);
      out3 = createMatrixValue( arg_out_2, nr, 1);
      out4 = createMatrixValue( arg_out_3, nr, 1);
      out5 = createMatrixValue( arg_out_4, nr, 1);
      out6 = createMatrixValue( arg_out_5, nr, 1);
    }
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
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_base5_DDDDD( int nlhs, mxArray       *plhs[],
                  int nrhs, mxArray const *prhs[] ) {
    #define CMD "BaseHermite_DDDDD('base5',t [,H]): "
    MEX_ASSERT( nrhs == 2 || nrhs == 3, CMD "expected 2 or 3  inputs, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 1 || nlhs == 6, CMD "expected 1 or 6  outputs, nlhs = " << nlhs );
    mwSize nr;
    double H = 1;
    if ( nrhs == 3 ) H = getScalarValue( arg_in_2, CMD "argument H" );
    double const * t = getVectorPointer( arg_in_1, nr, CMD "argument t");
    double *out1, *out2, *out3, *out4, *out5, *out6;
    if ( nlhs == 1 ) {
      out1 = createMatrixValue( arg_out_0, nr, 6);
      out2 = out1 + nr;
      out3 = out2 + nr;
      out4 = out3 + nr;
      out5 = out4 + nr;
      out6 = out5 + nr;
    } else if ( nlhs == 6 ) {
      out1 = createMatrixValue( arg_out_0, nr, 1);
      out2 = createMatrixValue( arg_out_1, nr, 1);
      out3 = createMatrixValue( arg_out_2, nr, 1);
      out4 = createMatrixValue( arg_out_3, nr, 1);
      out5 = createMatrixValue( arg_out_4, nr, 1);
      out6 = createMatrixValue( arg_out_5, nr, 1);
    }
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
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_eval5( int nlhs, mxArray       *plhs[],
            int nrhs, mxArray const *prhs[] ) {
    #define CMD "BaseHermite( 'eval5', t, P0, P1, T0, T1, J0, J1 [,H] ): "
    MEX_ASSERT( nrhs == 8 || nrhs == 9, CMD "expected 8 or 9 inputs, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
    mwSize nr, nr2, nr3, nr4, nr5, nr6, nr7;
    double H = 1;
    double const * t  = getVectorPointer( arg_in_1, nr,  CMD "argument t");
    double const * P0 = getVectorPointer( arg_in_2, nr2, CMD "argument t");
    double const * P1 = getVectorPointer( arg_in_3, nr3, CMD "argument t");
    double const * T0 = getVectorPointer( arg_in_4, nr4, CMD "argument t");
    double const * T1 = getVectorPointer( arg_in_5, nr5, CMD "argument t");
    double const * J0 = getVectorPointer( arg_in_6, nr6, CMD "argument t");
    double const * J1 = getVectorPointer( arg_in_7, nr7, CMD "argument t");
    MEX_ASSERT(
      nr2 == nr3 && nr3 == nr4 && nr4 == nr5 && nr5 == nr6 && nr6 == nr7,
      CMD " bad dimension |P0| = " << nr2 << " |P1| = " << nr3 <<
      " |T0| = " << nr4 << " |T1| = " << nr5 <<
      " |J0| = " << nr6 << " |J1| = " << nr7
    );
    if ( nrhs == 9 ) H = getScalarValue( arg_in_8, CMD "argument H" );
    double * out = createMatrixValue( arg_out_0, nr2, nr );
    double tmp[6];
    for ( mwSize i = 0; i < nr; ++i ) {
      Hermite5( t[i], H, tmp );
      for ( mwSize j = 0; j < nr2; ++j )
        *out++ = tmp[0] * P0[j] + tmp[1] * P1[j] +
                 tmp[2] * T0[j] + tmp[3] * T1[j] +
                 tmp[4] * J0[j] + tmp[5] * J1[j] ;
    }
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_eval5_D( int nlhs, mxArray       *plhs[],
              int nrhs, mxArray const *prhs[] ) {
    #define CMD "BaseHermite( 'eval5_D', t, P0, P1, T0, T1, J0, J1 [,H] ): "
    MEX_ASSERT( nrhs == 8 || nrhs == 9, CMD "expected 8 or 9 inputs, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
    mwSize nr, nr2, nr3, nr4, nr5, nr6, nr7;
    double H = 1;
    double const * t  = getVectorPointer( arg_in_1, nr,  CMD "argument t");
    double const * P0 = getVectorPointer( arg_in_2, nr2, CMD "argument t");
    double const * P1 = getVectorPointer( arg_in_3, nr3, CMD "argument t");
    double const * T0 = getVectorPointer( arg_in_4, nr4, CMD "argument t");
    double const * T1 = getVectorPointer( arg_in_5, nr5, CMD "argument t");
    double const * J0 = getVectorPointer( arg_in_6, nr6, CMD "argument t");
    double const * J1 = getVectorPointer( arg_in_7, nr7, CMD "argument t");
    MEX_ASSERT(
      nr2 == nr3 && nr3 == nr4 && nr4 == nr5 && nr5 == nr6 && nr6 == nr7,
      CMD " bad dimension |P0| = " << nr2 << " |P1| = " << nr3 <<
      " |T0| = " << nr4 << " |T1| = " << nr5 <<
      " |J0| = " << nr6 << " |J1| = " << nr7
    );
    if ( nrhs == 9 ) H = getScalarValue( arg_in_8, CMD "argument H" );
    double * out = createMatrixValue( arg_out_0, nr2, nr );
    double tmp[6];
    for ( mwSize i = 0; i < nr; ++i ) {
      Hermite5_D( t[i], H, tmp );
      for ( mwSize j = 0; j < nr2; ++j )
        *out++ = tmp[0] * P0[j] + tmp[1] * P1[j] +
                 tmp[2] * T0[j] + tmp[3] * T1[j] +
                 tmp[4] * J0[j] + tmp[5] * J1[j] ;
    }
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_eval5_DD( int nlhs, mxArray       *plhs[],
               int nrhs, mxArray const *prhs[] ) {
    #define CMD "BaseHermite( 'eval5_DD', t, P0, P1, T0, T1, J0, J1 [,H] ): "
    MEX_ASSERT( nrhs == 8 || nrhs == 9, CMD "expected 8 or 9 inputs, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
    mwSize nr, nr2, nr3, nr4, nr5, nr6, nr7;
    double H = 1;
    double const * t  = getVectorPointer( arg_in_1, nr,  CMD "argument t");
    double const * P0 = getVectorPointer( arg_in_2, nr2, CMD "argument t");
    double const * P1 = getVectorPointer( arg_in_3, nr3, CMD "argument t");
    double const * T0 = getVectorPointer( arg_in_4, nr4, CMD "argument t");
    double const * T1 = getVectorPointer( arg_in_5, nr5, CMD "argument t");
    double const * J0 = getVectorPointer( arg_in_6, nr6, CMD "argument t");
    double const * J1 = getVectorPointer( arg_in_7, nr7, CMD "argument t");
    MEX_ASSERT(
      nr2 == nr3 && nr3 == nr4 && nr4 == nr5 && nr5 == nr6 && nr6 == nr7,
      CMD " bad dimension |P0| = " << nr2 << " |P1| = " << nr3 <<
      " |T0| = " << nr4 << " |T1| = " << nr5 <<
      " |J0| = " << nr6 << " |J1| = " << nr7
    );
    if ( nrhs == 9 ) H = getScalarValue( arg_in_8, CMD "argument H" );
    double * out = createMatrixValue( arg_out_0, nr2, nr );
    double tmp[6];
    for ( mwSize i = 0; i < nr; ++i ) {
      Hermite5_DD( t[i], H, tmp );
      for ( mwSize j = 0; j < nr2; ++j )
        *out++ = tmp[0] * P0[j] + tmp[1] * P1[j] +
                 tmp[2] * T0[j] + tmp[3] * T1[j] +
                 tmp[4] * J0[j] + tmp[5] * J1[j] ;
    }
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_eval5_DDD( int nlhs, mxArray       *plhs[],
                int nrhs, mxArray const *prhs[] ) {
    #define CMD "BaseHermite( 'eval5_DDD', t, P0, P1, T0, T1, J0, J1 [,H] ): "
    MEX_ASSERT( nrhs == 8 || nrhs == 9, CMD "expected 8 or 9 inputs, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
    mwSize nr, nr2, nr3, nr4, nr5, nr6, nr7;
    double H = 1;
    double const * t  = getVectorPointer( arg_in_1, nr,  CMD "argument t");
    double const * P0 = getVectorPointer( arg_in_2, nr2, CMD "argument t");
    double const * P1 = getVectorPointer( arg_in_3, nr3, CMD "argument t");
    double const * T0 = getVectorPointer( arg_in_4, nr4, CMD "argument t");
    double const * T1 = getVectorPointer( arg_in_5, nr5, CMD "argument t");
    double const * J0 = getVectorPointer( arg_in_6, nr6, CMD "argument t");
    double const * J1 = getVectorPointer( arg_in_7, nr7, CMD "argument t");
    MEX_ASSERT(
      nr2 == nr3 && nr3 == nr4 && nr4 == nr5 && nr5 == nr6 && nr6 == nr7,
      CMD " bad dimension |P0| = " << nr2 << " |P1| = " << nr3 <<
      " |T0| = " << nr4 << " |T1| = " << nr5 <<
      " |J0| = " << nr6 << " |J1| = " << nr7
    );
    if ( nrhs == 9 ) H = getScalarValue( arg_in_8, CMD "argument H" );
    double * out = createMatrixValue( arg_out_0, nr2, nr );
    double tmp[6];
    for ( mwSize i = 0; i < nr; ++i ) {
      Hermite5_DDD( t[i], H, tmp );
      for ( mwSize j = 0; j < nr2; ++j )
        *out++ = tmp[0] * P0[j] + tmp[1] * P1[j] +
                 tmp[2] * T0[j] + tmp[3] * T1[j] +
                 tmp[4] * J0[j] + tmp[5] * J1[j] ;
    }
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_eval5_DDDD( int nlhs, mxArray       *plhs[],
                 int nrhs, mxArray const *prhs[] ) {
    #define CMD "BaseHermite( 'eval5_DDDD', t, P0, P1, T0, T1, J0, J1 [,H] ): "
    MEX_ASSERT( nrhs == 8 || nrhs == 9, CMD "expected 8 or 9 inputs, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
    mwSize nr, nr2, nr3, nr4, nr5, nr6, nr7;
    double H = 1;
    double const * t  = getVectorPointer( arg_in_1, nr,  CMD "argument t");
    double const * P0 = getVectorPointer( arg_in_2, nr2, CMD "argument t");
    double const * P1 = getVectorPointer( arg_in_3, nr3, CMD "argument t");
    double const * T0 = getVectorPointer( arg_in_4, nr4, CMD "argument t");
    double const * T1 = getVectorPointer( arg_in_5, nr5, CMD "argument t");
    double const * J0 = getVectorPointer( arg_in_6, nr6, CMD "argument t");
    double const * J1 = getVectorPointer( arg_in_7, nr7, CMD "argument t");
    MEX_ASSERT(
      nr2 == nr3 && nr3 == nr4 && nr4 == nr5 && nr5 == nr6 && nr6 == nr7,
      CMD " bad dimension |P0| = " << nr2 << " |P1| = " << nr3 <<
      " |T0| = " << nr4 << " |T1| = " << nr5 <<
      " |J0| = " << nr6 << " |J1| = " << nr7
    );
    if ( nrhs == 9 ) H = getScalarValue( arg_in_8, CMD "argument H" );
    double * out = createMatrixValue( arg_out_0, nr2, nr );
    double tmp[6];
    for ( mwSize i = 0; i < nr; ++i ) {
      Hermite5_DDDD( t[i], H, tmp );
      for ( mwSize j = 0; j < nr2; ++j )
        *out++ = tmp[0] * P0[j] + tmp[1] * P1[j] +
                 tmp[2] * T0[j] + tmp[3] * T1[j] +
                 tmp[4] * J0[j] + tmp[5] * J1[j] ;
    }
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_eval5_DDDDD( int nlhs, mxArray       *plhs[],
                  int nrhs, mxArray const *prhs[] ) {
    #define CMD "BaseHermite( 'eval5_DDDDD', t, P0, P1, T0, T1, J0, J1 [,H] ): "
    MEX_ASSERT( nrhs == 8 || nrhs == 9, CMD "expected 8 or 9 inputs, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
    mwSize nr, nr2, nr3, nr4, nr5, nr6, nr7;
    double H = 1;
    double const * t  = getVectorPointer( arg_in_1, nr,  CMD "argument t");
    double const * P0 = getVectorPointer( arg_in_2, nr2, CMD "argument t");
    double const * P1 = getVectorPointer( arg_in_3, nr3, CMD "argument t");
    double const * T0 = getVectorPointer( arg_in_4, nr4, CMD "argument t");
    double const * T1 = getVectorPointer( arg_in_5, nr5, CMD "argument t");
    double const * J0 = getVectorPointer( arg_in_6, nr6, CMD "argument t");
    double const * J1 = getVectorPointer( arg_in_7, nr7, CMD "argument t");
    MEX_ASSERT(
      nr2 == nr3 && nr3 == nr4 && nr4 == nr5 && nr5 == nr6 && nr6 == nr7,
      CMD " bad dimension |P0| = " << nr2 << " |P1| = " << nr3 <<
      " |T0| = " << nr4 << " |T1| = " << nr5 <<
      " |J0| = " << nr6 << " |J1| = " << nr7
    );
    if ( nrhs == 9 ) H = getScalarValue( arg_in_8, CMD "argument H" );
    double * out = createMatrixValue( arg_out_0, nr2, nr );
    double tmp[6];
    for ( mwSize i = 0; i < nr; ++i ) {
      Hermite5_DDDDD( t[i], H, tmp );
      for ( mwSize j = 0; j < nr2; ++j )
        *out++ = tmp[0] * P0[j] + tmp[1] * P1[j] +
                 tmp[2] * T0[j] + tmp[3] * T1[j] +
                 tmp[4] * J0[j] + tmp[5] * J1[j] ;
    }
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_cut( int nlhs, mxArray       *plhs[],
          int nrhs, mxArray const *prhs[] ) {
    #define CMD "BaseHermite( 'cut', a, b, P0, P1, T0, T1, [,H] ): "
    MEX_ASSERT( nrhs == 7 || nrhs == 8, CMD "expected 7 or 8 inputs, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 4, CMD "expected 4 output, nlhs = " << nlhs );
    mwSize nr3, nr4, nr5, nr6;
    double H = 1;
    double const   a  = getScalarValue( arg_in_1, CMD "argument a");
    double const   b  = getScalarValue( arg_in_2, CMD "argument b");
    double const * P0 = getVectorPointer( arg_in_3, nr3, CMD "argument t");
    double const * P1 = getVectorPointer( arg_in_4, nr4, CMD "argument t");
    double const * T0 = getVectorPointer( arg_in_5, nr5, CMD "argument t");
    double const * T1 = getVectorPointer( arg_in_6, nr6, CMD "argument t");
    MEX_ASSERT(
      nr3 == nr4 && nr4 == nr5 && nr5 == nr6,
      CMD " bad dimension |P0| = " << nr3 << " |P1| = " << nr4 <<
      " |T0| = " << nr5 << " |T1| = " << nr6
    );
    if ( nrhs == 8 ) H = getScalarValue( arg_in_7, CMD "argument H" );
    double * oP0 = createMatrixValue( arg_out_0, nr3, 1 );
    double * oP1 = createMatrixValue( arg_out_1, nr3, 1 );
    double * oT0 = createMatrixValue( arg_out_2, nr3, 1 );
    double * oT1 = createMatrixValue( arg_out_3, nr3, 1 );
    double A[4], DA[4], B[4], DB[4];
    Hermite3( a, H, A ); Hermite3_D( a, H, DA );
    Hermite3( b, H, B ); Hermite3_D( b, H, DB );
    for ( mwSize j = 0; j < nr3; ++j ) {
      oP0[j] = A[0]  * P0[j] + A[1]  * P1[j] + A[2]  * T0[j] + A[3]  * T1[j];
      oT0[j] = DA[0] * P0[j] + DA[1] * P1[j] + DA[2] * T0[j] + DA[3] * T1[j];
      oP1[j] = B[0]  * P0[j] + B[1]  * P1[j] + B[2]  * T0[j] + B[3]  * T1[j];
      oT1[j] = DB[0] * P0[j] + DB[1] * P1[j] + DB[2] * T0[j] + DB[3] * T1[j];
    }
    for ( mwSize j = 0; j < nr3; ++j ) {
      oT0[j] *= b-a;
      oT1[j] *= b-a;
    }
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_alen( int nlhs, mxArray       *plhs[],
           int nrhs, mxArray const *prhs[] ) {
    #define CMD "BaseHermite( 'approximate_length', P0, P1, T0, T1, [,H] ): "
    MEX_ASSERT( nrhs == 5 || nrhs == 6, CMD "expected 5 or 6 inputs, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
    mwSize nr1, nr2, nr3, nr4;
    double H = 1;
    double const * P0 = getVectorPointer( arg_in_1, nr1, CMD "argument P0");
    double const * P1 = getVectorPointer( arg_in_2, nr2, CMD "argument P1");
    double const * T0 = getVectorPointer( arg_in_3, nr3, CMD "argument T0");
    double const * T1 = getVectorPointer( arg_in_4, nr4, CMD "argument T1");
    MEX_ASSERT(
      nr1 == nr2 && nr2 == nr3 && nr3 == nr4,
      CMD " bad dimension |P0| = " << nr1 << " |P1| = " << nr2 <<
      " |T0| = " << nr3 << " |T1| = " << nr4
    );
    if ( nrhs == 6 ) H = getScalarValue( arg_in_5, CMD "argument H" );

    vector<double> P(nr1);

    double A[4];
    Hermite3( 0, H, A );
    for ( mwSize j = 0; j < nr1; ++j )
      P[j] = A[0]  * P0[j] + A[1]  * P1[j] + A[2]  * T0[j] + A[3]  * T1[j];

    double len = 0;
    for ( mwSize j = 1; j <= 100; ++j ) {
      double t = (j*H)/100;
      Hermite3( t, H, A );
      double dst2 = 0;
      for ( mwSize j = 0; j < nr1; ++j ) {
        double Pj    = A[0]  * P0[j] + A[1]  * P1[j] + A[2]  * T0[j] + A[3]  * T1[j];
        double delta = P[j]-Pj;
        dst2 += delta*delta;
        P[j] = Pj;
      }
      len += sqrt( dst2 );
    }
    setScalarValue( arg_out_0, len );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  typedef enum {
    CMD_BASE,
    CMD_BASE_D,
    CMD_BASE_DD,
    CMD_BASE_DDD,
    CMD_EVAL,
    CMD_EVAL_D,
    CMD_EVAL_DD,
    CMD_EVAL_DDD,
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
    CMD_BASE5_DDDDD,
    CMD_EVAL5,
    CMD_EVAL5_D,
    CMD_EVAL5_DD,
    CMD_EVAL5_DDD,
    CMD_EVAL5_DDDD,
    CMD_EVAL5_DDDDD,
    CMD_CUT,
    CMD_ALEN,
  } CMD_LIST;

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static map<string,unsigned> cmd_to_idx = {
    {"base",CMD_BASE},
    {"base_D",CMD_BASE_D},
    {"base_DD",CMD_BASE_DD},
    {"base_DDD",CMD_BASE_DDD},
    {"eval",CMD_EVAL},
    {"eval_D",CMD_EVAL_D},
    {"eval_DD",CMD_EVAL_DD},
    {"eval_DDD",CMD_EVAL_DDD},
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
    {"eval5",CMD_EVAL5},
    {"eval5_D",CMD_EVAL5_D},
    {"eval5_DD",CMD_EVAL5_DD},
    {"eval5_DDD",CMD_EVAL5_DDD},
    {"eval5_DDDD",CMD_EVAL5_DDDD},
    {"eval5_DDDDD",CMD_EVAL5_DDDDD},
    {"cut",CMD_CUT},
    {"approximate_length",CMD_ALEN},
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
      case CMD_EVAL:
        do_eval( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_EVAL_D:
        do_eval_D( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_EVAL_DD:
        do_eval_DD( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_EVAL_DDD:
        do_eval_DDD( nlhs, plhs, nrhs, prhs );
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
      case CMD_EVAL5:
        do_eval5( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_EVAL5_D:
        do_eval5_D( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_EVAL5_DD:
        do_eval5_DD( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_EVAL5_DDD:
        do_eval5_DDD( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_EVAL5_DDDD:
        do_eval5_DDDD( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_EVAL5_DDDDD:
        do_eval5_DDDDD( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_CUT:
        do_cut( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_ALEN:
        do_alen( nlhs, plhs, nrhs, prhs );
        break;
      }

    } catch ( exception const & e ) {
      mexErrMsgTxt(e.what());
    } catch (...) {
      mexErrMsgTxt("BaseHermite failed\n");
    }

  }
}
