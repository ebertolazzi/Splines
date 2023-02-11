/****************************************************************************\
  Copyright (c) Enrico Bertolazzi 2016
  All Rights Reserved.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the file license.txt for more details.
\****************************************************************************/

#include "Splines.hh"
#include "Utils_mex.hh"

#ifdef __clang__
  #pragma clang diagnostic ignored "-Wexit-time-destructors"
#endif

#include <unordered_map>

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
    #define MEX_ERROR_MESSAGE_1 "res = BaseHermite('base',t [,H])"
    #define CMD MEX_ERROR_MESSAGE_1

    UTILS_MEX_ASSERT( nrhs == 2 || nrhs == 3, CMD ": expected 2 or 3  inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1 || nlhs == 4, CMD ": expected 1 or 4  outputs, nlhs = {}\n", nlhs );
    mwSize nr;
    double H = 1;
    if ( nrhs == 3 ) H = Utils::mex_get_scalar_value( arg_in_2, CMD ": argument H" );
    double const * t = Utils::mex_vector_pointer( arg_in_1, nr, CMD ": argument t");
    double *out1 = nullptr;
    double *out2 = nullptr;
    double *out3 = nullptr;
    double *out4 = nullptr;
    if ( nlhs == 1 ) {
      out1 = Utils::mex_create_matrix_value( arg_out_0, nr, 4);
      out2 = out1 + nr;
      out3 = out2 + nr;
      out4 = out3 + nr;
    } else if ( nlhs == 4 ) {
      out1 = Utils::mex_create_matrix_value( arg_out_0, nr, 1);
      out2 = Utils::mex_create_matrix_value( arg_out_1, nr, 1);
      out3 = Utils::mex_create_matrix_value( arg_out_2, nr, 1);
      out4 = Utils::mex_create_matrix_value( arg_out_3, nr, 1);
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
    #define MEX_ERROR_MESSAGE_2 "res = BaseHermite('base_D',t [,H])"
    #define CMD MEX_ERROR_MESSAGE_2

    UTILS_MEX_ASSERT( nrhs == 2 || nrhs == 3, CMD ": expected 2 or 3  inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1 || nlhs == 4, CMD ": expected 1 or 4  outputs, nlhs = {}\n", nlhs );
    mwSize nr;
    double H = 1;
    if ( nrhs == 3 ) H = Utils::mex_get_scalar_value( arg_in_2, CMD ": argument H" );
    double const * t = Utils::mex_vector_pointer( arg_in_1, nr, CMD ": argument t");
    double *out1 = nullptr;
    double *out2 = nullptr;
    double *out3 = nullptr;
    double *out4 = nullptr;
    if ( nlhs == 1 ) {
      out1 = Utils::mex_create_matrix_value( arg_out_0, nr, 4);
      out2 = out1 + nr;
      out3 = out2 + nr;
      out4 = out3 + nr;
    } else if ( nlhs == 4 ) {
      out1 = Utils::mex_create_matrix_value( arg_out_0, nr, 1);
      out2 = Utils::mex_create_matrix_value( arg_out_1, nr, 1);
      out3 = Utils::mex_create_matrix_value( arg_out_2, nr, 1);
      out4 = Utils::mex_create_matrix_value( arg_out_3, nr, 1);
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
    #define MEX_ERROR_MESSAGE_3 "res = BaseHermite('base_DD',t [,H])"
    #define CMD MEX_ERROR_MESSAGE_3

    UTILS_MEX_ASSERT( nrhs == 2 || nrhs == 3, CMD ": expected 2 or 3  inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1 || nlhs == 4, CMD ": expected 1 or 4  outputs, nlhs = {}\n", nlhs );
    mwSize nr;
    double H = 1;
    if ( nrhs == 3 ) H = Utils::mex_get_scalar_value( arg_in_2, CMD ": argument H" );
    double const * t = Utils::mex_vector_pointer( arg_in_1, nr, CMD ": argument t");
    double *out1 = nullptr;
    double *out2 = nullptr;
    double *out3 = nullptr;
    double *out4 = nullptr;
    if ( nlhs == 1 ) {
      out1 = Utils::mex_create_matrix_value( arg_out_0, nr, 4);
      out2 = out1 + nr;
      out3 = out2 + nr;
      out4 = out3 + nr;
    } else if ( nlhs == 4 ) {
      out1 = Utils::mex_create_matrix_value( arg_out_0, nr, 1);
      out2 = Utils::mex_create_matrix_value( arg_out_1, nr, 1);
      out3 = Utils::mex_create_matrix_value( arg_out_2, nr, 1);
      out4 = Utils::mex_create_matrix_value( arg_out_3, nr, 1);
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
    #define MEX_ERROR_MESSAGE_4 "res = BaseHermite('base_DDD',t [,H])"
    #define CMD MEX_ERROR_MESSAGE_4

    UTILS_MEX_ASSERT( nrhs == 2 || nrhs == 3, CMD ": expected 2 or 3 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1 || nlhs == 4, CMD ": expected 1 or 4 outputs, nlhs = {}\n", nlhs );
    mwSize nr;
    double H = 1;
    if ( nrhs == 3 ) H = Utils::mex_get_scalar_value( arg_in_2, CMD ": argument H" );
    double const * t = Utils::mex_vector_pointer( arg_in_1, nr, CMD ": argument t");
    double *out1 = nullptr;
    double *out2 = nullptr;
    double *out3 = nullptr;
    double *out4 = nullptr;
    if ( nlhs == 1 ) {
      out1 = Utils::mex_create_matrix_value( arg_out_0, nr, 4);
      out2 = out1 + nr;
      out3 = out2 + nr;
      out4 = out3 + nr;
    } else if ( nlhs == 4 ) {
      out1 = Utils::mex_create_matrix_value( arg_out_0, nr, 1);
      out2 = Utils::mex_create_matrix_value( arg_out_1, nr, 1);
      out3 = Utils::mex_create_matrix_value( arg_out_2, nr, 1);
      out4 = Utils::mex_create_matrix_value( arg_out_3, nr, 1);
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
    #define MEX_ERROR_MESSAGE_5 "res = BaseHermite( 'eval', t, P0, P1, T0, T1 [,H] )"
    #define CMD MEX_ERROR_MESSAGE_5

    UTILS_MEX_ASSERT( nrhs == 6 || nrhs == 7, CMD ": expected 6 or 7 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    mwSize nr, nr2, nr3, nr4, nr5;
    double H = 1;
    double const * t  = Utils::mex_vector_pointer( arg_in_1, nr,  CMD ": argument t");
    double const * P0 = Utils::mex_vector_pointer( arg_in_2, nr2, CMD ": argument t");
    double const * P1 = Utils::mex_vector_pointer( arg_in_3, nr3, CMD ": argument t");
    double const * T0 = Utils::mex_vector_pointer( arg_in_4, nr4, CMD ": argument t");
    double const * T1 = Utils::mex_vector_pointer( arg_in_5, nr5, CMD ": argument t");
    UTILS_MEX_ASSERT(
      nr2 == nr3 && nr3 == nr4 && nr4 == nr5,
      CMD " bad dimension |P0| = {} |P1| = {} |T0| = {} |T1| = {}\n",
      nr2, nr3, nr4, nr5
    );
    if ( nrhs == 7 ) H = Utils::mex_get_scalar_value( arg_in_6, CMD ": argument H" );
    double * out = Utils::mex_create_matrix_value( arg_out_0, nr2, nr );
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
    #define MEX_ERROR_MESSAGE_6 "res = BaseHermite( 'eval_D', t, P0, P1, T0, T1 [,H] )"
    #define CMD MEX_ERROR_MESSAGE_6

    UTILS_MEX_ASSERT( nrhs == 6 || nrhs == 7, CMD ": expected 6 or 7 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    mwSize nr, nr2, nr3, nr4, nr5;
    double H = 1;
    double const * t  = Utils::mex_vector_pointer( arg_in_1, nr,  CMD ": argument t");
    double const * P0 = Utils::mex_vector_pointer( arg_in_2, nr2, CMD ": argument t");
    double const * P1 = Utils::mex_vector_pointer( arg_in_3, nr3, CMD ": argument t");
    double const * T0 = Utils::mex_vector_pointer( arg_in_4, nr4, CMD ": argument t");
    double const * T1 = Utils::mex_vector_pointer( arg_in_5, nr5, CMD ": argument t");
    UTILS_MEX_ASSERT(
      nr2 == nr3 && nr3 == nr4 && nr4 == nr5,
      CMD " bad dimension |P0| = {} |P1| = {} |T0| = {} |T1| = {}\n",
      nr2, nr3, nr4, nr5
    );
    if ( nrhs == 7 ) H = Utils::mex_get_scalar_value( arg_in_6, CMD ": argument H" );
    double * out = Utils::mex_create_matrix_value( arg_out_0, nr2, nr );
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
    #define MEX_ERROR_MESSAGE_7 "res = BaseHermite( 'eval_DD', t, P0, P1, T0, T1 [,H] )"
    #define CMD MEX_ERROR_MESSAGE_7

    UTILS_MEX_ASSERT( nrhs == 6 || nrhs == 7, CMD ": expected 6 or 7 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    mwSize nr, nr2, nr3, nr4, nr5;
    double H = 1;
    double const * t  = Utils::mex_vector_pointer( arg_in_1, nr,  CMD ": argument t");
    double const * P0 = Utils::mex_vector_pointer( arg_in_2, nr2, CMD ": argument P0");
    double const * P1 = Utils::mex_vector_pointer( arg_in_3, nr3, CMD ": argument P1");
    double const * T0 = Utils::mex_vector_pointer( arg_in_4, nr4, CMD ": argument T0");
    double const * T1 = Utils::mex_vector_pointer( arg_in_5, nr5, CMD ": argument t1");
    UTILS_MEX_ASSERT(
      nr2 == nr3 && nr3 == nr4 && nr4 == nr5,
      CMD " bad dimension |P0| = {} |P1| = {} |T0| = {} |T1| = {}\n",
      nr2, nr3, nr4, nr5
    );
    if ( nrhs == 7 ) H = Utils::mex_get_scalar_value( arg_in_6, CMD ": argument H" );
    double * out = Utils::mex_create_matrix_value( arg_out_0, nr2, nr );
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
    #define MEX_ERROR_MESSAGE_8 "res = BaseHermite( 'eval_DDD', t, P0, P1, T0, T1 [,H] )"
    #define CMD MEX_ERROR_MESSAGE_8

    UTILS_MEX_ASSERT( nrhs == 6 || nrhs == 7, CMD ": expected 6 or 7 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    mwSize nr, nr2, nr3, nr4, nr5;
    double H = 1;
    double const * t  = Utils::mex_vector_pointer( arg_in_1, nr,  CMD ": argument t");
    double const * P0 = Utils::mex_vector_pointer( arg_in_2, nr2, CMD ": argument P0");
    double const * P1 = Utils::mex_vector_pointer( arg_in_3, nr3, CMD ": argument P1");
    double const * T0 = Utils::mex_vector_pointer( arg_in_4, nr4, CMD ": argument T0");
    double const * T1 = Utils::mex_vector_pointer( arg_in_5, nr5, CMD ": argument T1");
    UTILS_MEX_ASSERT(
      nr2 == nr3 && nr3 == nr4 && nr4 == nr5,
      CMD " bad dimension |P0| = {} |P1| = {} |T0| = {} |T1| = {}\n",
      nr2, nr3, nr4, nr5
    );
    if ( nrhs == 7 ) H = Utils::mex_get_scalar_value( arg_in_6, CMD ": argument H" );
    double * out = Utils::mex_create_matrix_value( arg_out_0, nr2, nr );
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
    #define MEX_ERROR_MESSAGE_9 "[P0,P1,P2,P3] = BaseHermite('hermite_to_bezier',P0,P1,T0,T1)"
    #define CMD MEX_ERROR_MESSAGE_9

    UTILS_MEX_ASSERT( nrhs == 4, CMD ": expected 4 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 4, CMD ": expected 4 output, nlhs = {}\n", nlhs );
    mwSize n1, n2, n3, n4;
    double const * P0 = Utils::mex_vector_pointer( arg_in_1, n1, CMD ": P0" );
    double const * P1 = Utils::mex_vector_pointer( arg_in_2, n2, CMD ": P1" );
    double const * T0 = Utils::mex_vector_pointer( arg_in_3, n3, CMD ": T0" );
    double const * T1 = Utils::mex_vector_pointer( arg_in_4, n4, CMD ": T1" );

    UTILS_MEX_ASSERT(
      n1 == n2 && n2 == n3 && n3 == n4,
      CMD ": bad dimensions |P0| = {} |P1| = {} |T0| = {} |T1| = {}\n",
      n1, n2, n3, n4
    );
    double * oP0 = Utils::mex_create_matrix_value( arg_out_0, n1, 1);
    double * oP1 = Utils::mex_create_matrix_value( arg_out_1, n1, 1);
    double * oP2 = Utils::mex_create_matrix_value( arg_out_2, n1, 1);
    double * oP3 = Utils::mex_create_matrix_value( arg_out_3, n1, 1);
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
    #define MEX_ERROR_MESSAGE_10 "[P0,P1,T0,T1] = BaseHermite('bezier_to_hermite',P0,P1,P2,P3)"
    #define CMD MEX_ERROR_MESSAGE_10

    UTILS_MEX_ASSERT( nrhs == 4, CMD ": expected 4 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 4, CMD ": expected 4 output, nlhs = {}\n", nlhs );
    mwSize n1, n2, n3, n4;
    double const * P0 = Utils::mex_vector_pointer( arg_in_1, n1, CMD ": P0" );
    double const * P1 = Utils::mex_vector_pointer( arg_in_2, n2, CMD ": P1" );
    double const * P2 = Utils::mex_vector_pointer( arg_in_3, n3, CMD ": P2" );
    double const * P3 = Utils::mex_vector_pointer( arg_in_4, n4, CMD ": P3" );

    UTILS_MEX_ASSERT(
      n1 == n2 && n2 == n3 && n3 == n4,
      CMD ": bad dimensions |P0| = {} |P1| = {} |P2| = {} |P3| = {}\n",
      n1, n2, n3, n4
    );
    double * oP0 = Utils::mex_create_matrix_value( arg_out_0, n1, 1);
    double * oP1 = Utils::mex_create_matrix_value( arg_out_1, n1, 1);
    double * oT0 = Utils::mex_create_matrix_value( arg_out_2, n1, 1);
    double * oT1 = Utils::mex_create_matrix_value( arg_out_3, n1, 1);
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
                          int nrhs, mxArray const *[] ) {
    #define MEX_ERROR_MESSAGE_11 "res = BaseHermite('L2_first_derivative')"
    #define CMD MEX_ERROR_MESSAGE_11

    UTILS_MEX_ASSERT( nrhs == 1, CMD ": expected 1 input, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );

    double * sqrtD1 = Utils::mex_create_matrix_value( arg_out_0, 3, 4 );

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

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_L2_second_derivative( int nlhs, mxArray       *plhs[],
                           int nrhs, mxArray const *[] ) {
    #define MEX_ERROR_MESSAGE_12 "BaseHermite('L2_second_derivative',')"
    #define CMD MEX_ERROR_MESSAGE_12

    UTILS_MEX_ASSERT( nrhs == 1, CMD ": expected 1 input, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );

    double * sqrtD2 = Utils::mex_create_matrix_value( arg_out_0, 2, 4 );

    double t1 = sqrt(3.0);

    *sqrtD2++ = 0;
    *sqrtD2++ = 2*t1;

    *sqrtD2++ = 0;
    *sqrtD2++ = -2*t1;

    *sqrtD2++ = -1;
    *sqrtD2++ = t1;

    *sqrtD2++ = 1;
    *sqrtD2++ = t1;

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_L2_third_derivative( int nlhs, mxArray       *plhs[],
                          int nrhs, mxArray const *[] ) {
    #define MEX_ERROR_MESSAGE_13 "res = BaseHermite('L2_third_derivative',')"
    #define CMD MEX_ERROR_MESSAGE_13

    UTILS_MEX_ASSERT( nrhs == 1, CMD ": expected 1 input, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    double * sqrtD3 = Utils::mex_create_matrix_value( arg_out_0, 1, 4 );

    *sqrtD3++ = 12;
    *sqrtD3++ = -12;
    *sqrtD3++ = 6;
    *sqrtD3++ = 6;
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_base5( int nlhs, mxArray       *plhs[],
            int nrhs, mxArray const *prhs[] ) {
    #define MEX_ERROR_MESSAGE_14 "[b0,b1,b2,b3,b4,b5] = BaseHermite('base5',t [,H])"
    #define CMD MEX_ERROR_MESSAGE_14

    UTILS_MEX_ASSERT( nrhs == 2 || nrhs == 3, CMD ": expected 2 or 3  inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1 || nlhs == 6, CMD ": expected 1 or 6  outputs, nlhs = {}\n", nlhs );
    mwSize nr;
    double H = 1;
    if ( nrhs == 3 ) H = Utils::mex_get_scalar_value( arg_in_2, CMD "argument H" );
    double const * t = Utils::mex_vector_pointer( arg_in_1, nr, CMD "argument t");
    double *out1 = nullptr;
    double *out2 = nullptr;
    double *out3 = nullptr;
    double *out4 = nullptr;
    double *out5 = nullptr;
    double *out6 = nullptr;
    if ( nlhs == 1 ) {
      out1 = Utils::mex_create_matrix_value( arg_out_0, nr, 6);
      out2 = out1 + nr;
      out3 = out2 + nr;
      out4 = out3 + nr;
      out5 = out4 + nr;
      out6 = out5 + nr;
    } else if ( nlhs == 6 ) {
      out1 = Utils::mex_create_matrix_value( arg_out_0, nr, 1);
      out2 = Utils::mex_create_matrix_value( arg_out_1, nr, 1);
      out3 = Utils::mex_create_matrix_value( arg_out_2, nr, 1);
      out4 = Utils::mex_create_matrix_value( arg_out_3, nr, 1);
      out5 = Utils::mex_create_matrix_value( arg_out_4, nr, 1);
      out6 = Utils::mex_create_matrix_value( arg_out_5, nr, 1);
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
    #define MEX_ERROR_MESSAGE_15 "[b0,b1,b2,b3,b4,b5] = BaseHermite('base5_D',t [,H])"
    #define CMD MEX_ERROR_MESSAGE_15

    UTILS_MEX_ASSERT( nrhs == 2 || nrhs == 3, CMD ": expected 2 or 3  inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1 || nlhs == 6, CMD ": expected 1 or 6  outputs, nlhs = {}\n", nlhs );
    mwSize nr;
    double H = 1;
    if ( nrhs == 3 ) H = Utils::mex_get_scalar_value( arg_in_2, CMD ": argument H" );
    double const * t = Utils::mex_vector_pointer( arg_in_1, nr, CMD ": argument t");
    double *out1 = nullptr;
    double *out2 = nullptr;
    double *out3 = nullptr;
    double *out4 = nullptr;
    double *out5 = nullptr;
    double *out6 = nullptr;
    if ( nlhs == 1 ) {
      out1 = Utils::mex_create_matrix_value( arg_out_0, nr, 6);
      out2 = out1 + nr;
      out3 = out2 + nr;
      out4 = out3 + nr;
      out5 = out4 + nr;
      out6 = out5 + nr;
    } else if ( nlhs == 6 ) {
      out1 = Utils::mex_create_matrix_value( arg_out_0, nr, 1);
      out2 = Utils::mex_create_matrix_value( arg_out_1, nr, 1);
      out3 = Utils::mex_create_matrix_value( arg_out_2, nr, 1);
      out4 = Utils::mex_create_matrix_value( arg_out_3, nr, 1);
      out5 = Utils::mex_create_matrix_value( arg_out_4, nr, 1);
      out6 = Utils::mex_create_matrix_value( arg_out_5, nr, 1);
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
    #define MEX_ERROR_MESSAGE_16 "[b0,b1,b2,b3,b4,b5] = BaseHermite('base5_DD',t [,H])"
    #define CMD MEX_ERROR_MESSAGE_16

    UTILS_MEX_ASSERT( nrhs == 2 || nrhs == 3, CMD ": expected 2 or 3  inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1 || nlhs == 6, CMD ": expected 1 or 6  outputs, nlhs = {}\n", nlhs );
    mwSize nr;
    double H = 1;
    if ( nrhs == 3 ) H = Utils::mex_get_scalar_value( arg_in_2, CMD ": argument H" );
    double const * t = Utils::mex_vector_pointer( arg_in_1, nr, CMD ": argument t");
    double *out1 = nullptr;
    double *out2 = nullptr;
    double *out3 = nullptr;
    double *out4 = nullptr;
    double *out5 = nullptr;
    double *out6 = nullptr;
    if ( nlhs == 1 ) {
      out1 = Utils::mex_create_matrix_value( arg_out_0, nr, 6);
      out2 = out1 + nr;
      out3 = out2 + nr;
      out4 = out3 + nr;
      out5 = out4 + nr;
      out6 = out5 + nr;
    } else if ( nlhs == 6 ) {
      out1 = Utils::mex_create_matrix_value( arg_out_0, nr, 1);
      out2 = Utils::mex_create_matrix_value( arg_out_1, nr, 1);
      out3 = Utils::mex_create_matrix_value( arg_out_2, nr, 1);
      out4 = Utils::mex_create_matrix_value( arg_out_3, nr, 1);
      out5 = Utils::mex_create_matrix_value( arg_out_4, nr, 1);
      out6 = Utils::mex_create_matrix_value( arg_out_5, nr, 1);
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
    #define MEX_ERROR_MESSAGE_17 "[b0,b1,b2,b3,b4,b5] = BaseHermite('base5_DDD',t [,H])"
    #define CMD MEX_ERROR_MESSAGE_17

    UTILS_MEX_ASSERT( nrhs == 2 || nrhs == 3, CMD ": expected 2 or 3  inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1 || nlhs == 6, CMD ": expected 1 or 6  outputs, nlhs = {}\n", nlhs );
    mwSize nr;
    double H = 1;
    if ( nrhs == 3 ) H = Utils::mex_get_scalar_value( arg_in_2, CMD ": argument H" );
    double const * t = Utils::mex_vector_pointer( arg_in_1, nr, CMD ": argument t");
    double *out1 = nullptr;
    double *out2 = nullptr;
    double *out3 = nullptr;
    double *out4 = nullptr;
    double *out5 = nullptr;
    double *out6 = nullptr;
    if ( nlhs == 1 ) {
      out1 = Utils::mex_create_matrix_value( arg_out_0, nr, 6);
      out2 = out1 + nr;
      out3 = out2 + nr;
      out4 = out3 + nr;
      out5 = out4 + nr;
      out6 = out5 + nr;
    } else if ( nlhs == 6 ) {
      out1 = Utils::mex_create_matrix_value( arg_out_0, nr, 1);
      out2 = Utils::mex_create_matrix_value( arg_out_1, nr, 1);
      out3 = Utils::mex_create_matrix_value( arg_out_2, nr, 1);
      out4 = Utils::mex_create_matrix_value( arg_out_3, nr, 1);
      out5 = Utils::mex_create_matrix_value( arg_out_4, nr, 1);
      out6 = Utils::mex_create_matrix_value( arg_out_5, nr, 1);
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
    #define MEX_ERROR_MESSAGE_18 "[b0,b1,b2,b3,b4,b5] = BaseHermite_DDDD('base5',t [,H])"
    #define CMD MEX_ERROR_MESSAGE_18

    UTILS_MEX_ASSERT( nrhs == 2 || nrhs == 3, CMD ": expected 2 or 3  inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1 || nlhs == 6, CMD ": expected 1 or 6  outputs, nlhs = {}\n", nlhs );
    mwSize nr;
    double H = 1;
    if ( nrhs == 3 ) H = Utils::mex_get_scalar_value( arg_in_2, CMD ": argument H" );
    double const * t = Utils::mex_vector_pointer( arg_in_1, nr, CMD ": argument t");
    double *out1 = nullptr;
    double *out2 = nullptr;
    double *out3 = nullptr;
    double *out4 = nullptr;
    double *out5 = nullptr;
    double *out6 = nullptr;
    if ( nlhs == 1 ) {
      out1 = Utils::mex_create_matrix_value( arg_out_0, nr, 6);
      out2 = out1 + nr;
      out3 = out2 + nr;
      out4 = out3 + nr;
      out5 = out4 + nr;
      out6 = out5 + nr;
    } else if ( nlhs == 6 ) {
      out1 = Utils::mex_create_matrix_value( arg_out_0, nr, 1);
      out2 = Utils::mex_create_matrix_value( arg_out_1, nr, 1);
      out3 = Utils::mex_create_matrix_value( arg_out_2, nr, 1);
      out4 = Utils::mex_create_matrix_value( arg_out_3, nr, 1);
      out5 = Utils::mex_create_matrix_value( arg_out_4, nr, 1);
      out6 = Utils::mex_create_matrix_value( arg_out_5, nr, 1);
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
    #define MEX_ERROR_MESSAGE_19 "[b0,b1,b2,b3,b4,b5] = BaseHermite_DDDDD('base5',t [,H])"
    #define CMD MEX_ERROR_MESSAGE_19

    UTILS_MEX_ASSERT( nrhs == 2 || nrhs == 3, CMD ": expected 2 or 3  inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1 || nlhs == 6, CMD ": expected 1 or 6  outputs, nlhs = {}\n", nlhs );
    mwSize nr;
    double H = 1;
    if ( nrhs == 3 ) H = Utils::mex_get_scalar_value( arg_in_2, CMD ": argument H" );
    double const * t = Utils::mex_vector_pointer( arg_in_1, nr, CMD ": argument t");
    double *out1 = nullptr;
    double *out2 = nullptr;
    double *out3 = nullptr;
    double *out4 = nullptr;
    double *out5 = nullptr;
    double *out6 = nullptr;
    if ( nlhs == 1 ) {
      out1 = Utils::mex_create_matrix_value( arg_out_0, nr, 6);
      out2 = out1 + nr;
      out3 = out2 + nr;
      out4 = out3 + nr;
      out5 = out4 + nr;
      out6 = out5 + nr;
    } else if ( nlhs == 6 ) {
      out1 = Utils::mex_create_matrix_value( arg_out_0, nr, 1);
      out2 = Utils::mex_create_matrix_value( arg_out_1, nr, 1);
      out3 = Utils::mex_create_matrix_value( arg_out_2, nr, 1);
      out4 = Utils::mex_create_matrix_value( arg_out_3, nr, 1);
      out5 = Utils::mex_create_matrix_value( arg_out_4, nr, 1);
      out6 = Utils::mex_create_matrix_value( arg_out_5, nr, 1);
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
    #define MEX_ERROR_MESSAGE_20 "res = BaseHermite( 'eval5', t, P0, P1, T0, T1, J0, J1 [,H] )"
    #define CMD MEX_ERROR_MESSAGE_20

    UTILS_MEX_ASSERT( nrhs == 8 || nrhs == 9, CMD ": expected 8 or 9 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    mwSize nr, nr2, nr3, nr4, nr5, nr6, nr7;
    double H = 1;
    double const * t  = Utils::mex_vector_pointer( arg_in_1, nr,  CMD ": argument t");
    double const * P0 = Utils::mex_vector_pointer( arg_in_2, nr2, CMD ": argument P0");
    double const * P1 = Utils::mex_vector_pointer( arg_in_3, nr3, CMD ": argument P1");
    double const * T0 = Utils::mex_vector_pointer( arg_in_4, nr4, CMD ": argument T0");
    double const * T1 = Utils::mex_vector_pointer( arg_in_5, nr5, CMD ": argument T1");
    double const * J0 = Utils::mex_vector_pointer( arg_in_6, nr6, CMD ": argument J0");
    double const * J1 = Utils::mex_vector_pointer( arg_in_7, nr7, CMD ": argument J1");
    UTILS_MEX_ASSERT(
      nr2 == nr3 && nr3 == nr4 && nr4 == nr5 && nr5 == nr6 && nr6 == nr7,
      CMD ": bad dimensions\n"
      "|P0| = {} |P1| = {} |T0| = {} |T1| = {} |J0| = {} |J1| = {}",
      nr2, nr3, nr4, nr5, nr6, nr7
    );
    if ( nrhs == 9 ) H = Utils::mex_get_scalar_value( arg_in_8, CMD ": argument H" );
    double * out = Utils::mex_create_matrix_value( arg_out_0, nr2, nr );
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
    #define MEX_ERROR_MESSAGE_21 "res = BaseHermite( 'eval5_D', t, P0, P1, T0, T1, J0, J1 [,H] )"
    #define CMD MEX_ERROR_MESSAGE_21

    UTILS_MEX_ASSERT( nrhs == 8 || nrhs == 9, CMD ": expected 8 or 9 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    mwSize nr, nr2, nr3, nr4, nr5, nr6, nr7;
    double H = 1;
    double const * t  = Utils::mex_vector_pointer( arg_in_1, nr,  CMD ": argument t");
    double const * P0 = Utils::mex_vector_pointer( arg_in_2, nr2, CMD ": argument P0");
    double const * P1 = Utils::mex_vector_pointer( arg_in_3, nr3, CMD ": argument P1");
    double const * T0 = Utils::mex_vector_pointer( arg_in_4, nr4, CMD ": argument T0");
    double const * T1 = Utils::mex_vector_pointer( arg_in_5, nr5, CMD ": argument T1");
    double const * J0 = Utils::mex_vector_pointer( arg_in_6, nr6, CMD ": argument J0");
    double const * J1 = Utils::mex_vector_pointer( arg_in_7, nr7, CMD ": argument J1");
    UTILS_MEX_ASSERT(
      nr2 == nr3 && nr3 == nr4 && nr4 == nr5 && nr5 == nr6 && nr6 == nr7,
      CMD ": bad dimension |P0|={} |P1|={} |T0|={} |T1|={} |J0|={} |J1|={}\n",
      nr2, nr3, nr4, nr5, nr6, nr7
    );
    if ( nrhs == 9 ) H = Utils::mex_get_scalar_value( arg_in_8, CMD ": argument H" );
    double * out = Utils::mex_create_matrix_value( arg_out_0, nr2, nr );
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
    #define MEX_ERROR_MESSAGE_22 "res = BaseHermite( 'eval5_DD', t, P0, P1, T0, T1, J0, J1 [,H] )"
    #define CMD MEX_ERROR_MESSAGE_22

    UTILS_MEX_ASSERT( nrhs == 8 || nrhs == 9, CMD ": expected 8 or 9 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    mwSize nr, nr2, nr3, nr4, nr5, nr6, nr7;
    double H = 1;
    double const * t  = Utils::mex_vector_pointer( arg_in_1, nr,  CMD ": argument t");
    double const * P0 = Utils::mex_vector_pointer( arg_in_2, nr2, CMD ": argument P0");
    double const * P1 = Utils::mex_vector_pointer( arg_in_3, nr3, CMD ": argument P1");
    double const * T0 = Utils::mex_vector_pointer( arg_in_4, nr4, CMD ": argument T0");
    double const * T1 = Utils::mex_vector_pointer( arg_in_5, nr5, CMD ": argument T1");
    double const * J0 = Utils::mex_vector_pointer( arg_in_6, nr6, CMD ": argument J0");
    double const * J1 = Utils::mex_vector_pointer( arg_in_7, nr7, CMD ": argument J1");
    UTILS_MEX_ASSERT(
      nr2 == nr3 && nr3 == nr4 && nr4 == nr5 && nr5 == nr6 && nr6 == nr7,
      CMD ": bad dimension |P0|={} |P1|={} |T0|={} |T1|={} |J0|={} |J1|={}\n",
      nr2, nr3, nr4, nr5, nr6, nr7
    );
    if ( nrhs == 9 ) H = Utils::mex_get_scalar_value( arg_in_8, CMD ": argument H" );
    double * out = Utils::mex_create_matrix_value( arg_out_0, nr2, nr );
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
    #define MEX_ERROR_MESSAGE_23 "res = BaseHermite( 'eval5_DDD', t, P0, P1, T0, T1, J0, J1 [,H] )"
    #define CMD MEX_ERROR_MESSAGE_23

    UTILS_MEX_ASSERT( nrhs == 8 || nrhs == 9, CMD ": expected 8 or 9 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    mwSize nr, nr2, nr3, nr4, nr5, nr6, nr7;
    double H = 1;
    double const * t  = Utils::mex_vector_pointer( arg_in_1, nr,  CMD ": argument t");
    double const * P0 = Utils::mex_vector_pointer( arg_in_2, nr2, CMD ": argument P0");
    double const * P1 = Utils::mex_vector_pointer( arg_in_3, nr3, CMD ": argument P1");
    double const * T0 = Utils::mex_vector_pointer( arg_in_4, nr4, CMD ": argument T0");
    double const * T1 = Utils::mex_vector_pointer( arg_in_5, nr5, CMD ": argument T1");
    double const * J0 = Utils::mex_vector_pointer( arg_in_6, nr6, CMD ": argument J0");
    double const * J1 = Utils::mex_vector_pointer( arg_in_7, nr7, CMD ": argument J1");
    UTILS_MEX_ASSERT(
      nr2 == nr3 && nr3 == nr4 && nr4 == nr5 && nr5 == nr6 && nr6 == nr7,
      CMD ": bad dimension |P0|={} |P1|={} |T0|={} |T1|={} |J0|={} |J1|={}\n",
      nr2, nr3, nr4, nr5, nr6, nr7
    );
    if ( nrhs == 9 ) H = Utils::mex_get_scalar_value( arg_in_8, CMD ": argument H" );
    double * out = Utils::mex_create_matrix_value( arg_out_0, nr2, nr );
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
    #define MEX_ERROR_MESSAGE_24 "res = BaseHermite( 'eval5_DDDD', t, P0, P1, T0, T1, J0, J1 [,H] )"
    #define CMD MEX_ERROR_MESSAGE_24

    UTILS_MEX_ASSERT( nrhs == 8 || nrhs == 9, CMD ": expected 8 or 9 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    mwSize nr, nr2, nr3, nr4, nr5, nr6, nr7;
    double H = 1;
    double const * t  = Utils::mex_vector_pointer( arg_in_1, nr,  CMD ": argument t");
    double const * P0 = Utils::mex_vector_pointer( arg_in_2, nr2, CMD ": argument P0");
    double const * P1 = Utils::mex_vector_pointer( arg_in_3, nr3, CMD ": argument P1");
    double const * T0 = Utils::mex_vector_pointer( arg_in_4, nr4, CMD ": argument T0");
    double const * T1 = Utils::mex_vector_pointer( arg_in_5, nr5, CMD ": argument T1");
    double const * J0 = Utils::mex_vector_pointer( arg_in_6, nr6, CMD ": argument J0");
    double const * J1 = Utils::mex_vector_pointer( arg_in_7, nr7, CMD ": argument J1");
    UTILS_MEX_ASSERT(
      nr2 == nr3 && nr3 == nr4 && nr4 == nr5 && nr5 == nr6 && nr6 == nr7,
      CMD ": bad dimension |P0|={} |P1|={} |T0|={} |T1|={} |J0|={} |J1|={}\n",
      nr2, nr3, nr4, nr5, nr6, nr7
    );
    if ( nrhs == 9 ) H = Utils::mex_get_scalar_value( arg_in_8, CMD ": argument H" );
    double * out = Utils::mex_create_matrix_value( arg_out_0, nr2, nr );
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
    #define MEX_ERROR_MESSAGE_25 "res = BaseHermite( 'eval5_DDDDD', t, P0, P1, T0, T1, J0, J1 [,H] )"
    #define CMD MEX_ERROR_MESSAGE_25

    UTILS_MEX_ASSERT( nrhs == 8 || nrhs == 9, CMD ": expected 8 or 9 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    mwSize nr, nr2, nr3, nr4, nr5, nr6, nr7;
    double H = 1;
    double const * t  = Utils::mex_vector_pointer( arg_in_1, nr,  CMD ": argument t");
    double const * P0 = Utils::mex_vector_pointer( arg_in_2, nr2, CMD ": argument P0");
    double const * P1 = Utils::mex_vector_pointer( arg_in_3, nr3, CMD ": argument P1");
    double const * T0 = Utils::mex_vector_pointer( arg_in_4, nr4, CMD ": argument T0");
    double const * T1 = Utils::mex_vector_pointer( arg_in_5, nr5, CMD ": argument T1");
    double const * J0 = Utils::mex_vector_pointer( arg_in_6, nr6, CMD ": argument J0");
    double const * J1 = Utils::mex_vector_pointer( arg_in_7, nr7, CMD ": argument J1");
    UTILS_MEX_ASSERT(
      nr2 == nr3 && nr3 == nr4 && nr4 == nr5 && nr5 == nr6 && nr6 == nr7,
      CMD ": bad dimension |P0|={} |P1|={} |T0|={} |T1|={} |J0|={} |J1|={}\n",
      nr2, nr3, nr4, nr5, nr6, nr7
    );
    if ( nrhs == 9 ) H = Utils::mex_get_scalar_value( arg_in_8, CMD ": argument H" );
    double * out = Utils::mex_create_matrix_value( arg_out_0, nr2, nr );
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
    #define MEX_ERROR_MESSAGE_26 "[P0,P1,T0,T1] = BaseHermite( 'cut', a, b, P0, P1, T0, T1, [,H] )"
    #define CMD MEX_ERROR_MESSAGE_26

    UTILS_MEX_ASSERT( nrhs == 7 || nrhs == 8, CMD ": expected 7 or 8 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 4, CMD ": expected 4 output, nlhs = {}\n", nlhs );
    mwSize nr3, nr4, nr5, nr6;
    double H = 1;
    double const   a  = Utils::mex_get_scalar_value( arg_in_1, CMD ": argument a");
    double const   b  = Utils::mex_get_scalar_value( arg_in_2, CMD ": argument b");
    double const * P0 = Utils::mex_vector_pointer( arg_in_3, nr3, CMD ": argument P0");
    double const * P1 = Utils::mex_vector_pointer( arg_in_4, nr4, CMD ": argument P1");
    double const * T0 = Utils::mex_vector_pointer( arg_in_5, nr5, CMD ": argument T0");
    double const * T1 = Utils::mex_vector_pointer( arg_in_6, nr6, CMD ": argument T1");
    UTILS_MEX_ASSERT(
      nr3 == nr4 && nr4 == nr5 && nr5 == nr6,
      CMD ": bad dimension |P0| = {} |P1| = {} |T0| = {} |T1| = {}\n",
      nr3, nr4, nr5, nr6
    );
    if ( nrhs == 8 ) H = Utils::mex_get_scalar_value( arg_in_7, CMD ": argument H" );
    double * oP0 = Utils::mex_create_matrix_value( arg_out_0, nr3, 1 );
    double * oP1 = Utils::mex_create_matrix_value( arg_out_1, nr3, 1 );
    double * oT0 = Utils::mex_create_matrix_value( arg_out_2, nr3, 1 );
    double * oT1 = Utils::mex_create_matrix_value( arg_out_3, nr3, 1 );
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
  do_approximate_length( int nlhs, mxArray       *plhs[],
                         int nrhs, mxArray const *prhs[] ) {
    #define MEX_ERROR_MESSAGE_27 "res = BaseHermite( 'approximate_length', P0, P1, T0, T1, [,H] )"
    #define CMD MEX_ERROR_MESSAGE_27

    UTILS_MEX_ASSERT( nrhs == 5 || nrhs == 6, CMD ": expected 5 or 6 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    mwSize nr1, nr2, nr3, nr4;
    double H = 1;
    double const * P0 = Utils::mex_vector_pointer( arg_in_1, nr1, CMD ": argument P0");
    double const * P1 = Utils::mex_vector_pointer( arg_in_2, nr2, CMD ": argument P1");
    double const * T0 = Utils::mex_vector_pointer( arg_in_3, nr3, CMD ": argument T0");
    double const * T1 = Utils::mex_vector_pointer( arg_in_4, nr4, CMD ": argument T1");
    UTILS_MEX_ASSERT(
      nr1 == nr2 && nr2 == nr3 && nr3 == nr4,
      CMD " bad dimension |P0| = {} |P1| = {} |T0| = {} |T1| = {}\n",
      nr1, nr2, nr3, nr4
    );
    if ( nrhs == 6 ) H = Utils::mex_get_scalar_value( arg_in_5, CMD ": argument H" );

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
    Utils::mex_set_scalar_value( arg_out_0, len );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  typedef void (*DO_CMD)( int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[] );

  static unordered_map<string,DO_CMD> cmd_to_fun = {
    {"base",do_base},
    {"base_D",do_base_D},
    {"base_DD",do_base_DD},
    {"base_DDD",do_base_DDD},
    {"eval",do_eval},
    {"eval_D",do_eval_D},
    {"eval_DD",do_eval_DD},
    {"eval_DDD",do_eval_DDD},
    {"hermite_to_bezier",do_hermite_to_bezier},
    {"bezier_to_hermite",do_bezier_to_hermite},
    {"L2_first_derivative",do_L2_first_derivative},
    {"L2_second_derivative",do_L2_second_derivative},
    {"L2_third_derivative",do_L2_third_derivative},
    {"base5",do_base5},
    {"base5_D",do_base5_D},
    {"base5_DD",do_base5_DD},
    {"base5_DDD",do_base5_DDD},
    {"base5_DDDD",do_base5_DDDD},
    {"base5_DDDDD",do_base5_DDDDD},
    {"eval5",do_eval5},
    {"eval5_D",do_eval5_D},
    {"eval5_DD",do_eval5_DD},
    {"eval5_DDD",do_eval5_DDD},
    {"eval5_DDDD",do_eval5_DDDD},
    {"eval5_DDDDD",do_eval5_DDDDD},
    {"cut",do_cut},
    {"approximate_length",do_approximate_length},
  };

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  #define MEX_ERROR_MESSAGE \
"%======================================================================%\n" \
"BaseHermiteMexWrapper:  Compute Hermite base\n" \
"\n" \
"USAGE:\n" \
"\n" \
MEX_ERROR_MESSAGE_1 "\n" \
MEX_ERROR_MESSAGE_2 "\n" \
MEX_ERROR_MESSAGE_3 "\n" \
MEX_ERROR_MESSAGE_4 "\n" \
MEX_ERROR_MESSAGE_5 "\n" \
MEX_ERROR_MESSAGE_6 "\n" \
MEX_ERROR_MESSAGE_7 "\n" \
MEX_ERROR_MESSAGE_8 "\n" \
MEX_ERROR_MESSAGE_9 "\n" \
MEX_ERROR_MESSAGE_10 "\n" \
MEX_ERROR_MESSAGE_11 "\n" \
MEX_ERROR_MESSAGE_12 "\n" \
MEX_ERROR_MESSAGE_13 "\n" \
MEX_ERROR_MESSAGE_14 "\n" \
MEX_ERROR_MESSAGE_15 "\n" \
MEX_ERROR_MESSAGE_16 "\n" \
MEX_ERROR_MESSAGE_17 "\n" \
MEX_ERROR_MESSAGE_18 "\n" \
MEX_ERROR_MESSAGE_19 "\n" \
MEX_ERROR_MESSAGE_20 "\n" \
MEX_ERROR_MESSAGE_21 "\n" \
MEX_ERROR_MESSAGE_22 "\n" \
MEX_ERROR_MESSAGE_23 "\n" \
MEX_ERROR_MESSAGE_24 "\n" \
MEX_ERROR_MESSAGE_25 "\n" \
MEX_ERROR_MESSAGE_26 "\n" \
MEX_ERROR_MESSAGE_27 "\n" \
"\n" \
"%======================================================================%\n" \
"%                                                                      %\n" \
"%  Autor: Enrico Bertolazzi                                            %\n" \
"%         Department of Industrial Engineering                         %\n" \
"%         University of Trento                                         %\n" \
"%         enrico.bertolazzi@unitn.it                                   %\n" \
"%                                                                      %\n" \
"%======================================================================%\n"

  extern "C"
  void
  mexFunction( int nlhs, mxArray       *plhs[],
               int nrhs, mxArray const *prhs[] ) {

    char cmd[256];

    // the first argument must be a string
    if ( nrhs == 0 ) { mexErrMsgTxt(MEX_ERROR_MESSAGE); return; }

    try {
      UTILS_MEX_ASSERT0( mxIsChar(arg_in_0), "First argument must be a string" );
      mxGetString( arg_in_0, cmd, 256 );
      cmd_to_fun.at(cmd)( nlhs, plhs, nrhs, prhs );
    } catch ( exception const & e ) {
      mexErrMsgTxt( fmt::format( "BaseHermite Error: {}", e.what() ).c_str() );
    } catch (...) {
      mexErrMsgTxt("BaseHermite failed\n");
    }

  }
}
