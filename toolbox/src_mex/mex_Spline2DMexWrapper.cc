/****************************************************************************\
Copyright (c) 2015, Enrico Bertolazzi
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the FreeBSD Project.
\****************************************************************************/

#include "Splines.hh"
#include "Utils_mex.hh"

#define ASSERT(COND,MSG)                          \
  if ( !(COND) ) {                                \
    std::ostringstream ost;                       \
    ost << "Spline2DMexWrapper: " << MSG << '\n'; \
    mexErrMsgTxt(ost.str().c_str());              \
  }

#ifdef __clang__
  #pragma clang diagnostic ignored "-Wexit-time-destructors"
#endif

#include <unordered_map>

namespace Splines {

  using namespace std;

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_new( int nlhs, mxArray       *plhs[],
          int nrhs, mxArray const *prhs[] ) {

    #define MEX_ERROR_MESSAGE_1 "Spline1DMexWrapper( 'new', kind )"
    #define CMD MEX_ERROR_MESSAGE_1

    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );

    UTILS_MEX_ASSERT(
      mxIsChar(arg_in_1),
      CMD ": second argument must be a string, found ``{}''\n",
      mxGetClassName(arg_in_1)
    );

    string tname = mxArrayToString(arg_in_1);

    SplineSurf * v = nullptr;

    if      ( tname == "bilinear"  ) v = new Splines::BilinearSpline();
    else if ( tname == "bicubic"   ) v = new Splines::BiCubicSpline();
    else if ( tname == "akima"     ) v = new Splines::Akima2Dspline();
    else if ( tname == "biquintic" ) v = new Splines::BiQuinticSpline();
    else {
      UTILS_MEX_ASSERT0(
        false,
        CMD ": second argument must be one of the strings:\n"
        "'linear', 'cubic', 'akima', 'bessel', 'pchip', 'quintic'"
      );
    }

    arg_out_0 = Utils::mex_convert_ptr_to_mx<SplineSurf>( v );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_delete( int nlhs, mxArray       *[],
             int nrhs, mxArray const *prhs[] ) {

    #define MEX_ERROR_MESSAGE_2 "Spline1DMexWrapper( 'delete', OBJ )"
    #define CMD MEX_ERROR_MESSAGE_2

    UTILS_MEX_ASSERT( nlhs == 0, CMD ": expected 0 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 inputs, nrhs = {}\n", nrhs );

    // Destroy the C++ object
    Utils::mex_destroy_object<SplineSurf>(arg_in_1);
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_build( int nlhs, mxArray       *[],
            int nrhs, mxArray const *prhs[] ) {

    #define MEX_ERROR_MESSAGE_3 "Spline2DMexWrapper('build',OBJ,x,y,z)"
    #define CMD MEX_ERROR_MESSAGE_3

    UTILS_MEX_ASSERT( nlhs == 0, CMD ": expected no output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 5, CMD ": expected 5 input, nrhs = {}\n", nrhs );

    mwSize nx, ny, nnx, nny;
    real_type const * x = Utils::mex_vector_pointer(
      arg_in_2, nx, CMD ": error in reading 'x'"
    );
    real_type const * y = Utils::mex_vector_pointer(
      arg_in_3, ny, CMD ": error in reading 'y'"
    );
    real_type const * z = Utils::mex_matrix_pointer(
      arg_in_4, nnx, nny, CMD ": error in reading 'z'"
    );

    UTILS_MEX_ASSERT(
      nx == nnx,
      CMD ": lenght of 'x' ({}) must be the number of row of 'z' ({})\n",
      nx, nnx
    );

    UTILS_MEX_ASSERT(
      ny == nny,
      CMD "lenght of 'y' ({}) must be the number of column of 'z' ({})\n",
      ny, nny
    );

    bool fortran_storage = true;
    bool transposed      = true;

    integer ldZ = nx;

    SplineSurf * ptr = Utils::mex_convert_mx_to_ptr<SplineSurf>( arg_in_1 );

    ptr->build ( x, 1, y, 1, z, ldZ, nx, ny, fortran_storage, transposed );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_nx( int nlhs, mxArray       *plhs[],
         int nrhs, mxArray const *prhs[] ) {

    #define MEX_ERROR_MESSAGE_4 "Spline2DMexWrapper('nx',OBJ)"
    #define CMD MEX_ERROR_MESSAGE_4

    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 input, nrhs = {}\n", nrhs );

    SplineSurf * ptr = Utils::mex_convert_mx_to_ptr<SplineSurf>( arg_in_1 );

    Utils::mex_set_scalar_int32( arg_out_0, ptr->num_point_x() );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_ny( int nlhs, mxArray       *plhs[],
         int nrhs, mxArray const *prhs[] ) {

    #define MEX_ERROR_MESSAGE_5 "Spline2DMexWrapper('nx',OBJ)"
    #define CMD MEX_ERROR_MESSAGE_5

    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 input, nrhs = {}\n", nrhs );

    SplineSurf * ptr = Utils::mex_convert_mx_to_ptr<SplineSurf>( arg_in_1 );

    Utils::mex_set_scalar_int32( arg_out_0, ptr->num_point_y() );

    #undef CMD

  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_x_min( int nlhs, mxArray       *plhs[],
            int nrhs, mxArray const *prhs[] ) {

    #define MEX_ERROR_MESSAGE_6 "Spline2DMexWrapper('x_min',OBJ)"
    #define CMD MEX_ERROR_MESSAGE_6

    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 input, nrhs = {}\n", nrhs );

    SplineSurf * ptr = Utils::mex_convert_mx_to_ptr<SplineSurf>( arg_in_1 );

    Utils::mex_set_scalar_value( arg_out_0, ptr->x_min() );

    #undef CMD

  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_x_max( int nlhs, mxArray       *plhs[],
            int nrhs, mxArray const *prhs[] ) {

    #define MEX_ERROR_MESSAGE_7 "Spline2DMexWrapper('x_max',OBJ)"
    #define CMD MEX_ERROR_MESSAGE_7

    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 input, nrhs = {}\n", nrhs );

    SplineSurf * ptr = Utils::mex_convert_mx_to_ptr<SplineSurf>( arg_in_1 );

    Utils::mex_set_scalar_value( arg_out_0, ptr->x_max() );

    #undef CMD

  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_y_min( int nlhs, mxArray       *plhs[],
            int nrhs, mxArray const *prhs[] ) {

    #define MEX_ERROR_MESSAGE_8 "Spline2DMexWrapper('y_min',OBJ)"
    #define CMD MEX_ERROR_MESSAGE_8

    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 input, nrhs = {}\n", nrhs );

    SplineSurf * ptr = Utils::mex_convert_mx_to_ptr<SplineSurf>( arg_in_1 );

    Utils::mex_set_scalar_value( arg_out_0, ptr->y_min() );

    #undef CMD

  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_y_max( int nlhs, mxArray       *plhs[],
            int nrhs, mxArray const *prhs[] ) {

    #define MEX_ERROR_MESSAGE_9 "Spline2DMexWrapper('y_max',OBJ)"
    #define CMD MEX_ERROR_MESSAGE_9

    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 input, nrhs = {}\n", nrhs );

    SplineSurf * ptr = Utils::mex_convert_mx_to_ptr<SplineSurf>( arg_in_1 );

    Utils::mex_set_scalar_value( arg_out_0, ptr->y_max() );

    #undef CMD

  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_zMin( int nlhs, mxArray       *plhs[],
           int nrhs, mxArray const *prhs[] ) {

    #define MEX_ERROR_MESSAGE_10 "Spline2DMexWrapper('zMin',OBJ)"
    #define CMD MEX_ERROR_MESSAGE_10

    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 input, nrhs = {}\n", nrhs );

    SplineSurf * ptr = Utils::mex_convert_mx_to_ptr<SplineSurf>( arg_in_1 );

    Utils::mex_set_scalar_value( arg_out_0, ptr->z_min() );

    #undef CMD

  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_zMax( int nlhs, mxArray       *plhs[],
           int nrhs, mxArray const *prhs[] ) {

    #define MEX_ERROR_MESSAGE_11 "Spline2DMexWrapper('zMax',OBJ)"
    #define CMD MEX_ERROR_MESSAGE_11

    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 input, nrhs = {}\n", nrhs );

    SplineSurf * ptr = Utils::mex_convert_mx_to_ptr<SplineSurf>( arg_in_1 );

    Utils::mex_set_scalar_value( arg_out_0, ptr->z_max() );

    #undef CMD

  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_eval( int nlhs, mxArray       *plhs[],
           int nrhs, mxArray const *prhs[] ) {

    #define MEX_ERROR_MESSAGE_12 "Spline2DMexWrapper('eval',OBJ,x,y)"
    #define CMD MEX_ERROR_MESSAGE_12

    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 4, CMD ": expected 4 input, nrhs = {}\n", nrhs );

    SplineSurf * ptr = Utils::mex_convert_mx_to_ptr<SplineSurf>( arg_in_1 );

    mwSize nx, mx, ny, my;
    real_type const * x = Utils::mex_matrix_pointer(
      arg_in_2, nx, mx, CMD ": error in reading `x`"
    );
    real_type const * y = Utils::mex_matrix_pointer(
      arg_in_3, ny, my, CMD ": error in reading `y`"
    );
    UTILS_MEX_ASSERT(
      nx == ny && mx == my,
      CMD ": size(x) = {} x {} must be equal to size(y) = {} x {}",
      nx, ny, mx, my
    );
    real_type * z = Utils::mex_create_matrix_value( arg_out_0, nx, mx );
    for ( mwSize j = 0; j < mx*nx; ++j )
      *z++ = (*ptr)( *x++, *y++ );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_eval_Dx( int nlhs, mxArray       *plhs[],
              int nrhs, mxArray const *prhs[] ) {

    #define MEX_ERROR_MESSAGE_13 "Spline2DMexWrapper('eval_Dx',OBJ,x,y)"
    #define CMD MEX_ERROR_MESSAGE_13

    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 4, CMD ": expected 4 input, nrhs = {}\n", nrhs );

    SplineSurf * ptr = Utils::mex_convert_mx_to_ptr<SplineSurf>( arg_in_1 );

    mwSize nx, mx, ny, my;
    real_type const * x = Utils::mex_matrix_pointer(
      arg_in_2, nx, mx, CMD ": error in reading `x`"
    );
    real_type const * y = Utils::mex_matrix_pointer(
      arg_in_3, ny, my, CMD ": error in reading `y`"
    );
    UTILS_MEX_ASSERT(
      nx == ny && mx == my,
      CMD ": size(x) = {} x {} must be equal to size(y) = {} x {}\n",
      nx, ny, mx, my
    );
    real_type * z = Utils::mex_create_matrix_value( arg_out_0, nx, mx );
    for ( mwSize j = 0; j < mx*nx; ++j )
      *z++ = ptr->Dx( *x++, *y++ );
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_eval_Dy( int nlhs, mxArray       *plhs[],
              int nrhs, mxArray const *prhs[] ) {

    #define MEX_ERROR_MESSAGE_14 "Spline2DMexWrapper('eval_Dy',OBJ,x,y)"
    #define CMD MEX_ERROR_MESSAGE_14

    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 4, CMD ": expected 4 input, nrhs = {}\n", nrhs );

    SplineSurf * ptr = Utils::mex_convert_mx_to_ptr<SplineSurf>( arg_in_1 );

    mwSize nx, mx, ny, my;
    real_type const * x = Utils::mex_matrix_pointer(
      arg_in_2, nx, mx, CMD ": error in reading `x`"
    );
    real_type const * y = Utils::mex_matrix_pointer(
      arg_in_3, ny, my, CMD ": error in reading `y`"
    );
    UTILS_MEX_ASSERT(
      nx == ny && mx == my,
      CMD ": size(x) = {} x {} must be equal to size(y) = {} x {}",
      nx, ny, mx, my
    );
    real_type * z = Utils::mex_create_matrix_value( arg_out_0, nx, mx );
    for ( mwSize j = 0; j < mx*nx; ++j )
      *z++ = ptr->Dy( *x++, *y++ );
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_eval_Dxx( int nlhs, mxArray       *plhs[],
               int nrhs, mxArray const *prhs[] ) {

    #define MEX_ERROR_MESSAGE_15 "Spline2DMexWrapper('eval_Dxx',OBJ,x,y)"
    #define CMD MEX_ERROR_MESSAGE_15

    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 4, CMD ": expected 4 input, nrhs = {}\n", nrhs );

    SplineSurf * ptr = Utils::mex_convert_mx_to_ptr<SplineSurf>( arg_in_1 );

    mwSize nx, mx, ny, my;
    real_type const * x = Utils::mex_matrix_pointer(
      arg_in_2, nx, mx, CMD ": error in reading `x`"
    );
    real_type const * y = Utils::mex_matrix_pointer(
      arg_in_3, ny, my, CMD ": error in reading `y`"
    );
    UTILS_MEX_ASSERT(
      nx == ny && mx == my,
      CMD ": size(x) = {} x {} must be equal to size(y) = {} x {}\n",
      nx, ny, mx, my
    );
    real_type * z = Utils::mex_create_matrix_value( arg_out_0, nx, mx );
    for ( mwSize j = 0; j < mx*nx; ++j )
      *z++ = ptr->Dxx( *x++, *y++ );
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_eval_Dxy( int nlhs, mxArray       *plhs[],
               int nrhs, mxArray const *prhs[] ) {

    #define MEX_ERROR_MESSAGE_16 "Spline2DMexWrapper('eval_Dxy',OBJ,x,y)"
    #define CMD MEX_ERROR_MESSAGE_16

    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 4, CMD ": expected 4 input, nrhs = {}\n", nrhs );

    SplineSurf * ptr = Utils::mex_convert_mx_to_ptr<SplineSurf>( arg_in_1 );

    mwSize nx, mx, ny, my;
    real_type const * x = Utils::mex_matrix_pointer(
      arg_in_2, nx, mx, CMD ": error in reading `x`"
    );
    real_type const * y = Utils::mex_matrix_pointer(
      arg_in_3, ny, my, CMD ": error in reading `y`"
    );
    UTILS_MEX_ASSERT(
      nx == ny && mx == my,
      CMD ": size(x) = {} x {} must be equal to size(y) = {} x {}\n",
      nx, ny, mx, my
    );
    real_type * z = Utils::mex_create_matrix_value( arg_out_0, nx, mx );
    for ( mwSize j = 0; j < mx*nx; ++j )
      *z++ = ptr->Dxy( *x++, *y++ );
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_eval_Dyy( int nlhs, mxArray       *plhs[],
               int nrhs, mxArray const *prhs[] ) {

    #define MEX_ERROR_MESSAGE_17 "Spline2DMexWrapper('eval_Dyy',OBJ,x,y)"
    #define CMD MEX_ERROR_MESSAGE_17

    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 4, CMD ": expected 4 input, nrhs = {}\n", nrhs );

    SplineSurf * ptr = Utils::mex_convert_mx_to_ptr<SplineSurf>( arg_in_1 );

    mwSize nx, mx, ny, my;
    real_type const * x = Utils::mex_matrix_pointer(
      arg_in_2, nx, mx, CMD ": error in reading `x`"
    );
    real_type const * y = Utils::mex_matrix_pointer(
      arg_in_3, ny, my, CMD ": error in reading `y`"
    );
    UTILS_MEX_ASSERT(
      nx == ny && mx == my,
      CMD ": size(x) = {} x {} must be equal to size(y) = {} x {}\n",
      nx, ny, mx, my
    );
    real_type * z = Utils::mex_create_matrix_value( arg_out_0, nx, mx );
    for ( mwSize j = 0; j < mx*nx; ++j )
      *z++ = ptr->Dyy( *x++, *y++ );
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_make_x_closed( int nlhs, mxArray       *[],
                    int nrhs, mxArray const *prhs[] ) {

    #define MEX_ERROR_MESSAGE_18 "Spline1DMexWrapper('make_x_closed',OBJ)"
    #define CMD MEX_ERROR_MESSAGE_18

    UTILS_MEX_ASSERT( nlhs == 0, CMD ": expected 0 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 input, nrhs = {}\n", nrhs );

    SplineSurf * ptr = Utils::mex_convert_mx_to_ptr<SplineSurf>( arg_in_1 );
    ptr->make_x_closed();

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_make_x_opened( int nlhs, mxArray       *[],
                    int nrhs, mxArray const *prhs[] ) {

    #define MEX_ERROR_MESSAGE_19 "Spline1DMexWrapper('make_x_opened',OBJ)"
    #define CMD MEX_ERROR_MESSAGE_19

    UTILS_MEX_ASSERT( nlhs == 0, CMD ": expected 0 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 input, nrhs = {}\n", nrhs );

    SplineSurf * ptr = Utils::mex_convert_mx_to_ptr<SplineSurf>( arg_in_1 );
    ptr->make_x_opened();

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_is_x_closed( int nlhs, mxArray       *plhs[],
                  int nrhs, mxArray const *prhs[] ) {

    #define MEX_ERROR_MESSAGE_20 "Spline1DMexWrapper('is_x_closed',OBJ)"
    #define CMD MEX_ERROR_MESSAGE_20

    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 input, nrhs = {}\n", nrhs );

    SplineSurf * ptr = Utils::mex_convert_mx_to_ptr<SplineSurf>( arg_in_1 );
    Utils::mex_set_scalar_bool( arg_out_0, ptr->is_x_closed() );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_make_x_bounded( int nlhs, mxArray       *[],
                     int nrhs, mxArray const *prhs[] ) {

    #define MEX_ERROR_MESSAGE_21 "Spline1DMexWrapper('make_x_bounded',OBJ)"
    #define CMD MEX_ERROR_MESSAGE_21

    UTILS_MEX_ASSERT( nlhs == 0, CMD ": expected 0 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 input, nrhs = {}\n", nrhs );

    SplineSurf * ptr = Utils::mex_convert_mx_to_ptr<SplineSurf>( arg_in_1 );
    ptr->make_y_bounded();

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_make_x_unbounded( int nlhs, mxArray       *[],
                       int nrhs, mxArray const *prhs[] ) {

    #define MEX_ERROR_MESSAGE_22 "Spline1DMexWrapper('make_x_unbounded',OBJ)"
    #define CMD MEX_ERROR_MESSAGE_22

    UTILS_MEX_ASSERT( nlhs == 0, CMD ": expected 0 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 input, nrhs = {}\n", nrhs );

    SplineSurf * ptr = Utils::mex_convert_mx_to_ptr<SplineSurf>( arg_in_1 );
    ptr->make_x_unbounded();

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_is_x_bounded( int nlhs, mxArray       *plhs[],
                   int nrhs, mxArray const *prhs[] ) {

    #define MEX_ERROR_MESSAGE_23 "Spline1DMexWrapper('is_x_bounded',OBJ)"
    #define CMD MEX_ERROR_MESSAGE_23

    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 input, nrhs = {}\n", nrhs );

    SplineSurf * ptr = Utils::mex_convert_mx_to_ptr<SplineSurf>( arg_in_1 );
    Utils::mex_set_scalar_bool( arg_out_0, ptr->is_x_bounded() );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_make_y_closed( int nlhs, mxArray       *[],
                    int nrhs, mxArray const *prhs[] ) {

    #define MEX_ERROR_MESSAGE_24 "Spline1DMexWrapper('make_y_closed',OBJ)"
    #define CMD MEX_ERROR_MESSAGE_24

    UTILS_MEX_ASSERT( nlhs == 0, CMD ": expected 0 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 input, nrhs = {}\n", nrhs );

    SplineSurf * ptr = Utils::mex_convert_mx_to_ptr<SplineSurf>( arg_in_1 );
    ptr->make_y_closed();

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_make_y_opened( int nlhs, mxArray       *[],
                    int nrhs, mxArray const *prhs[] ) {

    #define MEX_ERROR_MESSAGE_25 "Spline1DMexWrapper('make_y_opened',OBJ)"
    #define CMD MEX_ERROR_MESSAGE_25

    UTILS_MEX_ASSERT( nlhs == 0, CMD ": expected 0 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 input, nrhs = {}\n", nrhs );

    SplineSurf * ptr = Utils::mex_convert_mx_to_ptr<SplineSurf>( arg_in_1 );
    ptr->make_y_opened();

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_is_y_closed( int nlhs, mxArray       *plhs[],
                  int nrhs, mxArray const *prhs[] ) {

    #define MEX_ERROR_MESSAGE_26 "Spline1DMexWrapper('is_y_closed',OBJ)"
    #define CMD MEX_ERROR_MESSAGE_26

    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 input, nrhs = {}\n", nrhs );

    SplineSurf * ptr = Utils::mex_convert_mx_to_ptr<SplineSurf>( arg_in_1 );
    Utils::mex_set_scalar_bool( arg_out_0, ptr->is_y_closed() );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_make_y_bounded( int nlhs, mxArray       *[],
                     int nrhs, mxArray const *prhs[] ) {

    #define MEX_ERROR_MESSAGE_27 "Spline1DMexWrapper('make_y_bounded',OBJ)"
    #define CMD MEX_ERROR_MESSAGE_27

    UTILS_MEX_ASSERT( nlhs == 0, CMD ": expected 0 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 input, nrhs = {}\n", nrhs );

    SplineSurf * ptr = Utils::mex_convert_mx_to_ptr<SplineSurf>( arg_in_1 );
    ptr->make_y_bounded();

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_make_y_unbounded( int nlhs, mxArray       *[],
                       int nrhs, mxArray const *prhs[] ) {

    #define MEX_ERROR_MESSAGE_28 "Spline1DMexWrapper('make_y_unbounded',OBJ)"
    #define CMD MEX_ERROR_MESSAGE_28

    UTILS_MEX_ASSERT( nlhs == 0, CMD ": expected 0 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 input, nrhs = {}\n", nrhs );

    SplineSurf * ptr = Utils::mex_convert_mx_to_ptr<SplineSurf>( arg_in_1 );
    ptr->make_y_unbounded();

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_is_y_bounded( int nlhs, mxArray       *plhs[],
                   int nrhs, mxArray const *prhs[] ) {

    #define MEX_ERROR_MESSAGE_29 "Spline1DMexWrapper('is_y_bounded',OBJ)"
    #define CMD MEX_ERROR_MESSAGE_29

    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 input, nrhs = {}\n", nrhs );

    SplineSurf * ptr = Utils::mex_convert_mx_to_ptr<SplineSurf>( arg_in_1 );
    Utils::mex_set_scalar_bool( arg_out_0, ptr->is_y_bounded() );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  typedef void (*DO_CMD)( int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[] );

  static unordered_map<string,DO_CMD> cmd_to_fun = {
    {"new",do_new},
    {"delete",do_delete},
    {"build",do_build},
    {"nx",do_nx},
    {"ny",do_ny},
    {"x_min",do_x_min},
    {"x_max",do_x_max},
    {"y_min",do_y_min},
    {"y_max",do_y_max},
    {"zMin",do_zMin},
    {"zMax",do_zMax},
    {"eval",do_eval},
    {"eval_Dx",do_eval_Dx},
    {"eval_Dy",do_eval_Dy},
    {"eval_Dxx",do_eval_Dxx},
    {"eval_Dxy",do_eval_Dxy},
    {"eval_Dyy",do_eval_Dyy},
    {"is_x_closed",do_is_x_closed},
    {"is_y_closed",do_is_y_closed},
    {"make_x_closed",do_make_x_closed},
    {"make_y_closed",do_make_y_closed},
    {"make_x_opened",do_make_x_opened},
    {"make_y_opened",do_make_y_opened},
    {"is_x_bounded",do_is_x_bounded},
    {"is_y_bounded",do_is_y_bounded},
    {"make_x_bounded",do_make_x_bounded},
    {"make_y_bounded",do_make_y_bounded},
    {"make_x_unbounded",do_make_x_unbounded},
    {"make_y_unbounded",do_make_y_unbounded}
  };

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

#define MEX_ERROR_MESSAGE \
"%======================================================================%\n" \
"Spline2DMexWrapper:  Compute spline surface\n" \
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
MEX_ERROR_MESSAGE_28 "\n" \
MEX_ERROR_MESSAGE_29 "\n" \
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
      mexErrMsgTxt( fmt::format( "Spline2DMexWrapper Error: {}", e.what() ).c_str() );
    } catch (...) {
      mexErrMsgTxt("Spline2DMexWrapper failed\n");
    }

  }


}
