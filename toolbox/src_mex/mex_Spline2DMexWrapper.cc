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

#define MEX_ERROR_MESSAGE \
"%======================================================================%\n" \
"% Spline2DMexWrapper:  Compute spline curve                            %\n" \
"%                                                                      %\n" \
"% USAGE:                                                               %\n" \
"%   obj = Spline2DMexWrapper( 'new', kind );                           %\n" \
"%   Spline2DMexWrapper( 'delete', obj );                               %\n" \
"%   Spline2DMexWrapper( 'build', obj, X, Y, Z );                       %\n" \
"%   P    = Spline2DMexWrapper( 'eval', obj, X, Y );                    %\n" \
"%   Dx   = Spline2DMexWrapper( 'eval_Dx', obj, X, Y );                 %\n" \
"%   Dy   = Spline2DMexWrapper( 'eval_Dy', obj, X, Y );                 %\n" \
"%   Dxx  = Spline2DMexWrapper( 'eval_Dxx', obj, X, Y );                %\n" \
"%   Dxy  = Spline2DMexWrapper( 'eval_Dxy', obj, X, Y );                %\n" \
"%   Dyy  = Spline2DMexWrapper( 'eval_Dyy', obj, X, Y );                %\n" \
"%                                                                      %\n" \
"% On input:                                                            %\n" \
"%                                                                      %\n" \
"%  kind = string with the kind of spline, any of:                      %\n" \
"%         'bilinear', 'bicubic', 'akima', 'biquintic'                  %\n" \
"%  X = vector of X coordinates                                         %\n" \
"%  Y = vector of Y coordinates                                         %\n" \
"%  Z = matrix of Z coordinates                                         %\n" \
"%                                                                      %\n" \
"% On output:                                                           %\n" \
"%                                                                      %\n" \
"%  P    = vector of Z values                                           %\n" \
"%                                                                      %\n" \
"%======================================================================%\n" \
"%                                                                      %\n" \
"%  Autor: Enrico Bertolazzi                                            %\n" \
"%         Department of Industrial Engineering                         %\n" \
"%         University of Trento                                         %\n" \
"%         enrico.bertolazzi@unitn.it                                   %\n" \
"%                                                                      %\n" \
"%======================================================================%\n"


using namespace std;

namespace Splines {

  static
  void
  DATA_NEW( mxArray * & mx_id, SplineSurf * ptr ) {
    mx_id = Utils::convert_ptr_to_mat<SplineSurf>(ptr);
  }

  static
  inline
  SplineSurf *
  DATA_GET( mxArray const * & mx_id ) {
    return Utils::convert_ptr_to_mat<SplineSurf>(mx_id);
  }

  static
  void
  DATA_DELETE( mxArray const * & mx_id ) {
    Utils::destroy_object<SplineSurf>(mx_id);
  }


  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_new( int nlhs, mxArray       *plhs[],
          int nrhs, mxArray const *prhs[] ) {

    #define CMD "Spline1DMexWrapper( 'new', kind ): "
    MEX_ASSERT2( nrhs == 2, CMD "expected 2 inputs, nrhs = {}\n", nrhs );
    MEX_ASSERT2( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );

    MEX_ASSERT2(
      mxIsChar(arg_in_1),
      CMD "second argument must be a string, found ``{}''\n",
      mxGetClassName(arg_in_1)
    );

    string tname = mxArrayToString(arg_in_1);

    SplineSurf * v = nullptr;

    if      ( tname == "bilinear"  ) v = new Splines::BilinearSpline();
    else if ( tname == "bicubic"   ) v = new Splines::BiCubicSpline();
    else if ( tname == "akima"     ) v = new Splines::Akima2Dspline();
    else if ( tname == "biquintic" ) v = new Splines::BiQuinticSpline();
    else {
      MEX_ASSERT(
        false,
        "Second argument must be one of the strings:\n"
        "'linear', 'cubic', 'akima', 'bessel', 'pchip', 'quintic'"
      );
    }

    DATA_NEW( arg_out_0, v );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_delete( int nlhs, mxArray       *plhs[],
             int nrhs, mxArray const *prhs[] ) {

    #define CMD "Spline1DMexWrapper( 'delete', OBJ ): "
    MEX_ASSERT2( nrhs == 2, CMD "expected 2 inputs, nrhs = {}\n", nrhs );
    MEX_ASSERT2( nlhs == 0, CMD "expected 0 output, nlhs = {}\n", nlhs );

    // Destroy the C++ object
    DATA_DELETE(arg_in_1);
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_build( int nlhs, mxArray       *plhs[],
            int nrhs, mxArray const *prhs[] ) {

    #define CMD "Spline2DMexWrapper('build',OBJ,x,y,z): "

    MEX_ASSERT2( nlhs == 0, CMD "expected no output, nlhs = {}\n", nlhs );
    MEX_ASSERT2( nrhs == 5, CMD "expected 5 input, nrhs = {}\n", nrhs );

    mwSize nx, ny, nnx, nny;
    real_type const * x = Utils::mex_vector_pointer(
      arg_in_2, nx, CMD "error in reading 'x'"
    );
    real_type const * y = Utils::mex_vector_pointer(
      arg_in_3, ny, CMD "error in reading 'y'"
    );
    real_type const * z = Utils::mex_matrix_pointer(
      arg_in_4, nnx, nny, CMD "error in reading 'z'"
    );

    MEX_ASSERT2(
      nx == nnx,
      CMD "lenght of 'x' ({}) must be the number of row of 'z' ({})\n",
      nx, nnx
    );

    MEX_ASSERT2(
      ny == nny,
      CMD "lenght of 'y' ({}) must be the number of column of 'z' ({})\n",
      ny, nny
    );

    bool fortran_storage = true;
    bool transposed      = true;

    integer ldZ = nx;

    SplineSurf * ptr = DATA_GET( arg_in_1 );

    ptr->build ( x, 1, y, 1, z, ldZ, nx, ny, fortran_storage, transposed );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_nx( int nlhs, mxArray       *plhs[],
         int nrhs, mxArray const *prhs[] ) {

    #define CMD "Spline2DMexWrapper('nx',OBJ): "

    MEX_ASSERT2( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );
    MEX_ASSERT2( nrhs == 2, CMD "expected 2 input, nrhs = {}\n", nrhs );

    SplineSurf * ptr = DATA_GET( arg_in_1 );

    Utils::mex_set_scalar_int32( arg_out_0, ptr->numPointX() );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_ny( int nlhs, mxArray       *plhs[],
         int nrhs, mxArray const *prhs[] ) {

    #define CMD "Spline2DMexWrapper('nx',OBJ): "

    MEX_ASSERT2( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );
    MEX_ASSERT2( nrhs == 2, CMD "expected 2 input, nrhs = {}\n", nrhs );

    SplineSurf * ptr = DATA_GET( arg_in_1 );

    Utils::mex_set_scalar_int32( arg_out_0, ptr->numPointY() );

    #undef CMD

  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_xMin( int nlhs, mxArray       *plhs[],
           int nrhs, mxArray const *prhs[] ) {

    #define CMD "Spline2DMexWrapper('xMin',OBJ): "

    MEX_ASSERT2( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );
    MEX_ASSERT2( nrhs == 2, CMD "expected 2 input, nrhs = {}\n", nrhs );

    SplineSurf * ptr = DATA_GET( arg_in_1 );

    Utils::mex_set_scalar_value( arg_out_0, ptr->xMin() );

    #undef CMD

  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_xMax( int nlhs, mxArray       *plhs[],
           int nrhs, mxArray const *prhs[] ) {

    #define CMD "Spline2DMexWrapper('xMax',OBJ): "

    MEX_ASSERT2( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );
    MEX_ASSERT2( nrhs == 2, CMD "expected 2 input, nrhs = {}\n", nrhs );

    SplineSurf * ptr = DATA_GET( arg_in_1 );

    Utils::mex_set_scalar_value( arg_out_0, ptr->xMax() );

    #undef CMD

  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_yMin( int nlhs, mxArray       *plhs[],
           int nrhs, mxArray const *prhs[] ) {

    #define CMD "Spline2DMexWrapper('yMin',OBJ): "

    MEX_ASSERT2( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );
    MEX_ASSERT2( nrhs == 2, CMD "expected 2 input, nrhs = {}\n", nrhs );

    SplineSurf * ptr = DATA_GET( arg_in_1 );

    Utils::mex_set_scalar_value( arg_out_0, ptr->yMin() );

    #undef CMD

  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_yMax( int nlhs, mxArray       *plhs[],
           int nrhs, mxArray const *prhs[] ) {

    #define CMD "Spline2DMexWrapper('yMax',OBJ): "

    MEX_ASSERT2( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );
    MEX_ASSERT2( nrhs == 2, CMD "expected 2 input, nrhs = {}\n", nrhs );

    SplineSurf * ptr = DATA_GET( arg_in_1 );

    Utils::mex_set_scalar_value( arg_out_0, ptr->yMax() );

    #undef CMD

  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_zMin( int nlhs, mxArray       *plhs[],
           int nrhs, mxArray const *prhs[] ) {

    #define CMD "Spline2DMexWrapper('zMin',OBJ): "

    MEX_ASSERT2( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );
    MEX_ASSERT2( nrhs == 2, CMD "expected 2 input, nrhs = {}\n", nrhs );

    SplineSurf * ptr = DATA_GET( arg_in_1 );

    Utils::mex_set_scalar_value( arg_out_0, ptr->zMin() );

    #undef CMD

  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_zMax( int nlhs, mxArray       *plhs[],
           int nrhs, mxArray const *prhs[] ) {

    #define CMD "Spline2DMexWrapper('zMax',OBJ): "

    MEX_ASSERT2( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );
    MEX_ASSERT2( nrhs == 2, CMD "expected 2 input, nrhs = {}\n", nrhs );

    SplineSurf * ptr = DATA_GET( arg_in_1 );

    Utils::mex_set_scalar_value( arg_out_0, ptr->zMax() );

    #undef CMD

  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_eval( int nlhs, mxArray       *plhs[],
           int nrhs, mxArray const *prhs[] ) {

    #define CMD "Spline2DMexWrapper('eval',OBJ,x,y): "

    MEX_ASSERT2( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );
    MEX_ASSERT2( nrhs == 4, CMD "expected 4 input, nrhs = {}\n", nrhs );

    SplineSurf * ptr = DATA_GET( arg_in_1 );

    mwSize nx, mx, ny, my;
    real_type const * x = Utils::mex_matrix_pointer(
      arg_in_2, nx, mx, CMD "error in reading `x`"
    );
    real_type const * y = Utils::mex_matrix_pointer(
      arg_in_3, ny, my, CMD "error in reading `y`"
    );
    MEX_ASSERT(
      nx == ny && mx == my,
      CMD "size(x) = " << nx << " x " << ny <<
      " must be equal to size(y) = " << mx << " x " << my
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

    #define CMD "Spline2DMexWrapper('eval_Dx',OBJ,x,y): "

    MEX_ASSERT2( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );
    MEX_ASSERT2( nrhs == 4, CMD "expected 4 input, nrhs = {}\n", nrhs );

    SplineSurf * ptr = DATA_GET( arg_in_1 );

    mwSize nx, mx, ny, my;
    real_type const * x = Utils::mex_matrix_pointer(
      arg_in_2, nx, mx, CMD "error in reading `x`"
    );
    real_type const * y = Utils::mex_matrix_pointer(
      arg_in_3, ny, my, CMD "error in reading `y`"
    );
    MEX_ASSERT2(
      nx == ny && mx == my,
      CMD "size(x) = {} x {} must be equal to size(y) = {} x {}\n",
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

    #define CMD "Spline2DMexWrapper('eval_Dy',OBJ,x,y): "

    MEX_ASSERT2( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );
    MEX_ASSERT2( nrhs == 4, CMD "expected 4 input, nrhs = {}\n", nrhs );

    SplineSurf * ptr = DATA_GET( arg_in_1 );

    mwSize nx, mx, ny, my;
    real_type const * x = Utils::mex_matrix_pointer(
      arg_in_2, nx, mx, CMD "error in reading `x`"
    );
    real_type const * y = Utils::mex_matrix_pointer(
      arg_in_3, ny, my, CMD "error in reading `y`"
    );
    MEX_ASSERT(
      nx == ny && mx == my,
      CMD "size(x) = " << nx << " x " << ny <<
      " must be equal to size(y) = " << mx << " x " << my
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

    #define CMD "Spline2DMexWrapper('eval_Dxx',OBJ,x,y): "

    MEX_ASSERT2( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );
    MEX_ASSERT2( nrhs == 4, CMD "expected 4 input, nrhs = {}\n", nrhs );

    SplineSurf * ptr = DATA_GET( arg_in_1 );

    mwSize nx, mx, ny, my;
    real_type const * x = Utils::mex_matrix_pointer(
      arg_in_2, nx, mx, CMD "error in reading `x`"
    );
    real_type const * y = Utils::mex_matrix_pointer(
      arg_in_3, ny, my, CMD "error in reading `y`"
    );
    MEX_ASSERT2(
      nx == ny && mx == my,
      CMD "size(x) = {} x {} must be equal to size(y) = {} x {}\n",
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

    #define CMD "Spline2DMexWrapper('eval_Dxy',OBJ,x,y): "

    MEX_ASSERT2( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );
    MEX_ASSERT2( nrhs == 4, CMD "expected 4 input, nrhs = {}\n", nrhs );

    SplineSurf * ptr = DATA_GET( arg_in_1 );

    mwSize nx, mx, ny, my;
    real_type const * x = Utils::mex_matrix_pointer(
      arg_in_2, nx, mx, CMD "error in reading `x`"
    );
    real_type const * y = Utils::mex_matrix_pointer(
      arg_in_3, ny, my, CMD "error in reading `y`"
    );
    MEX_ASSERT2(
      nx == ny && mx == my,
      CMD "size(x) = {} x {} must be equal to size(y) = {} x {}\n",
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

    #define CMD "Spline2DMexWrapper('eval_Dyy',OBJ,x,y): "

    MEX_ASSERT2( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );
    MEX_ASSERT2( nrhs == 4, CMD "expected 4 input, nrhs = {}\n", nrhs );

    SplineSurf * ptr = DATA_GET( arg_in_1 );

    mwSize nx, mx, ny, my;
    real_type const * x = Utils::mex_matrix_pointer(
      arg_in_2, nx, mx, CMD "error in reading `x`"
    );
    real_type const * y = Utils::mex_matrix_pointer(
      arg_in_3, ny, my, CMD "error in reading `y`"
    );
    MEX_ASSERT2(
      nx == ny && mx == my,
      CMD "size(x) = {} x {} must be equal to size(y) = {} x {}\n",
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
  do_make_x_closed( int nlhs, mxArray       *plhs[],
                    int nrhs, mxArray const *prhs[] ) {

    #define CMD "Spline1DMexWrapper('make_x_closed',OBJ): "

    MEX_ASSERT2( nlhs == 0, CMD "expected 0 output, nlhs = {}\n", nlhs );
    MEX_ASSERT2( nrhs == 2, CMD "expected 2 input, nrhs = {}\n", nrhs );

    SplineSurf * ptr = DATA_GET( arg_in_1 );
    ptr->make_x_closed();

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_make_x_opened( int nlhs, mxArray       *plhs[],
                  int nrhs, mxArray const *prhs[] ) {

    #define CMD "Spline1DMexWrapper('make_x_opened',OBJ): "

    MEX_ASSERT2( nlhs == 0, CMD "expected 0 output, nlhs = {}\n", nlhs );
    MEX_ASSERT2( nrhs == 2, CMD "expected 2 input, nrhs = {}\n", nrhs );

    SplineSurf * ptr = DATA_GET( arg_in_1 );
    ptr->make_x_opened();

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_is_x_closed( int nlhs, mxArray       *plhs[],
                  int nrhs, mxArray const *prhs[] ) {

    #define CMD "Spline1DMexWrapper('is_x_closed',OBJ): "

    MEX_ASSERT2( nlhs == 0, CMD "expected 1 output, nlhs = {}\n", nlhs );
    MEX_ASSERT2( nrhs == 2, CMD "expected 2 input, nrhs = {}\n", nrhs );

    SplineSurf * ptr = DATA_GET( arg_in_1 );
    Utils::mex_set_scalar_bool( arg_out_0, ptr->is_x_closed() );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_make_x_bounded( int nlhs, mxArray       *plhs[],
                     int nrhs, mxArray const *prhs[] ) {

    #define CMD "Spline1DMexWrapper('make_x_bounded',OBJ): "

    MEX_ASSERT2( nlhs == 0, CMD "expected 0 output, nlhs = {}\n", nlhs );
    MEX_ASSERT2( nrhs == 2, CMD "expected 2 input, nrhs = {}\n", nrhs );

    SplineSurf * ptr = DATA_GET( arg_in_1 );
    ptr->make_y_bounded();

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_make_x_unbounded( int nlhs, mxArray       *plhs[],
                       int nrhs, mxArray const *prhs[] ) {

    #define CMD "Spline1DMexWrapper('make_x_unbounded',OBJ): "

    MEX_ASSERT2( nlhs == 0, CMD "expected 0 output, nlhs = {}\n", nlhs );
    MEX_ASSERT2( nrhs == 2, CMD "expected 2 input, nrhs = {}\n", nrhs );

    SplineSurf * ptr = DATA_GET( arg_in_1 );
    ptr->make_x_unbounded();

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_is_x_bounded( int nlhs, mxArray       *plhs[],
                   int nrhs, mxArray const *prhs[] ) {

    #define CMD "Spline1DMexWrapper('is_x_bounded',OBJ): "

    MEX_ASSERT2( nlhs == 0, CMD "expected 1 output, nlhs = {}\n", nlhs );
    MEX_ASSERT2( nrhs == 2, CMD "expected 2 input, nrhs = {}\n", nrhs );

    SplineSurf * ptr = DATA_GET( arg_in_1 );
    Utils::mex_set_scalar_bool( arg_out_0, ptr->is_x_bounded() );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_make_y_closed( int nlhs, mxArray       *plhs[],
                    int nrhs, mxArray const *prhs[] ) {

    #define CMD "Spline1DMexWrapper('make_y_closed',OBJ): "

    MEX_ASSERT2( nlhs == 0, CMD "expected 0 output, nlhs = {}\n", nlhs );
    MEX_ASSERT2( nrhs == 2, CMD "expected 2 input, nrhs = {}\n", nrhs );

    SplineSurf * ptr = DATA_GET( arg_in_1 );
    ptr->make_y_closed();

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_make_y_opened( int nlhs, mxArray       *plhs[],
                  int nrhs, mxArray const *prhs[] ) {

    #define CMD "Spline1DMexWrapper('make_y_opened',OBJ): "

    MEX_ASSERT2( nlhs == 0, CMD "expected 0 output, nlhs = {}\n", nlhs );
    MEX_ASSERT2( nrhs == 2, CMD "expected 2 input, nrhs = {}\n", nrhs );

    SplineSurf * ptr = DATA_GET( arg_in_1 );
    ptr->make_y_opened();

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_is_y_closed( int nlhs, mxArray       *plhs[],
                  int nrhs, mxArray const *prhs[] ) {

    #define CMD "Spline1DMexWrapper('is_y_closed',OBJ): "

    MEX_ASSERT2( nlhs == 0, CMD "expected 1 output, nlhs = {}\n", nlhs );
    MEX_ASSERT2( nrhs == 2, CMD "expected 2 input, nrhs = {}\n", nrhs );

    SplineSurf * ptr = DATA_GET( arg_in_1 );
    Utils::mex_set_scalar_bool( arg_out_0, ptr->is_y_closed() );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_make_y_bounded( int nlhs, mxArray       *plhs[],
                     int nrhs, mxArray const *prhs[] ) {

    #define CMD "Spline1DMexWrapper('make_y_bounded',OBJ): "

    MEX_ASSERT2( nlhs == 0, CMD "expected 0 output, nlhs = {}\n", nlhs );
    MEX_ASSERT2( nrhs == 2, CMD "expected 2 input, nrhs = {}\n", nrhs );

    SplineSurf * ptr = DATA_GET( arg_in_1 );
    ptr->make_y_bounded();

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_make_y_unbounded( int nlhs, mxArray       *plhs[],
                       int nrhs, mxArray const *prhs[] ) {

    #define CMD "Spline1DMexWrapper('make_y_unbounded',OBJ): "

    MEX_ASSERT2( nlhs == 0, CMD "expected 0 output, nlhs = {}\n", nlhs );
    MEX_ASSERT2( nrhs == 2, CMD "expected 2 input, nrhs = {}\n", nrhs );

    SplineSurf * ptr = DATA_GET( arg_in_1 );
    ptr->make_y_unbounded();

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_is_y_bounded( int nlhs, mxArray       *plhs[],
                   int nrhs, mxArray const *prhs[] ) {

    #define CMD "Spline1DMexWrapper('is_y_bounded',OBJ): "

    MEX_ASSERT2( nlhs == 0, CMD "expected 1 output, nlhs = {}\n", nlhs );
    MEX_ASSERT2( nrhs == 2, CMD "expected 2 input, nrhs = {}\n", nrhs );

    SplineSurf * ptr = DATA_GET( arg_in_1 );
    Utils::mex_set_scalar_bool( arg_out_0, ptr->is_y_bounded() );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  typedef void (*DO_CMD)( int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[] );

  static map<string,DO_CMD> cmd_to_fun = {
    {"new",do_new},
    {"delete",do_delete},
    {"build",do_build},
    {"nx",do_nx},
    {"ny",do_ny},
    {"xMin",do_xMin},
    {"xMax",do_xMax},
    {"yMin",do_yMin},
    {"yMax",do_yMax},
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
      DO_CMD pfun = cmd_to_fun.at(cmd);
      pfun( nlhs, plhs, nrhs, prhs );
    } catch ( exception const & e ) {
      mexErrMsgTxt( fmt::format( "Spline2DMexWrapper Error: {}", e.what() ).c_str() );
    } catch (...) {
      mexErrMsgTxt("Spline2DMexWrapper failed\n");
    }

  }


}
