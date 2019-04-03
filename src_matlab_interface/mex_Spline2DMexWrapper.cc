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
#include "mex_utils.hh"

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
    mx_id = convertPtr2Mat<SplineSurf>(ptr);
  }

  static
  inline
  SplineSurf *
  DATA_GET( mxArray const * & mx_id ) {
    return convertMat2Ptr<SplineSurf>(mx_id);
  }

  static
  void
  DATA_DELETE( mxArray const * & mx_id ) {
    destroyObject<SplineSurf>(mx_id);
  }


  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_new( int nlhs, mxArray       *plhs[],
          int nrhs, mxArray const *prhs[] ) {

    #define CMD "Spline1DMexWrapper( 'new', kind ): "
    MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );

    MEX_ASSERT( mxIsChar(arg_in_1),
                CMD "second argument must be a string, found ``" <<
                mxGetClassName(arg_in_1) << "''" );

    string tname = mxArrayToString(arg_in_1);

    SplineSurf * v = nullptr;

    if      ( tname == "bilinear"  ) v = new Splines::BilinearSpline();
    else if ( tname == "bicubic"   ) v = new Splines::BiCubicSpline();
    else if ( tname == "akima"     ) v = new Splines::Akima2Dspline();
    else if ( tname == "biquintic" ) v = new Splines::BiQuinticSpline();
    else {
      MEX_ASSERT(
        false,
        "Second argument must be one of the strings:\n" <<
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
    MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 0, CMD "expected 0 output, nlhs = " << nlhs );

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

    MEX_ASSERT( nlhs == 0, CMD "expected no output, nlhs = " << nlhs );
    MEX_ASSERT( nrhs == 5, CMD "expected 5 input, nrhs = " << nrhs );

    mwSize nx, ny, nnx, nny;
    real_type const * x = getVectorPointer(
      arg_in_2, nx, CMD "error in reading 'x'"
    );
    real_type const * y = getVectorPointer(
      arg_in_3, ny, CMD "error in reading 'y'"
    );
    real_type const * z = getMatrixPointer(
      arg_in_4, nnx, nny, CMD "error in reading 'z'"
    );

    MEX_ASSERT(
      nx == nnx,
      CMD "lenght of 'x' (" << nx << ") must be the number of row of 'z' (" << nnx << ")"
    );

    MEX_ASSERT(
      ny == nny,
      CMD "lenght of 'y' (" << ny << ") must be the number of column of 'z' (" << nny << ")"
    );

    bool fortran_storage = true;
    bool transposed      = false;

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

    MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
    MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );

    SplineSurf * ptr = DATA_GET( arg_in_1 );

    setScalarInt( arg_out_0, ptr->numPointX() );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_ny( int nlhs, mxArray       *plhs[],
         int nrhs, mxArray const *prhs[] ) {

    #define CMD "Spline2DMexWrapper('nx',OBJ): "

    MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
    MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );

    SplineSurf * ptr = DATA_GET( arg_in_1 );

    setScalarInt( arg_out_0, ptr->numPointY() );

    #undef CMD

  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_xMin( int nlhs, mxArray       *plhs[],
           int nrhs, mxArray const *prhs[] ) {

    #define CMD "Spline2DMexWrapper('xMin',OBJ): "

    MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
    MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );

    SplineSurf * ptr = DATA_GET( arg_in_1 );

    setScalarValue( arg_out_0, ptr->xMin() );

    #undef CMD

  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_xMax( int nlhs, mxArray       *plhs[],
           int nrhs, mxArray const *prhs[] ) {

    #define CMD "Spline2DMexWrapper('xMax',OBJ): "

    MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
    MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );

    SplineSurf * ptr = DATA_GET( arg_in_1 );

    setScalarValue( arg_out_0, ptr->xMax() );

    #undef CMD

  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_yMin( int nlhs, mxArray       *plhs[],
           int nrhs, mxArray const *prhs[] ) {

    #define CMD "Spline2DMexWrapper('yMin',OBJ): "

    MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
    MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );

    SplineSurf * ptr = DATA_GET( arg_in_1 );

    setScalarValue( arg_out_0, ptr->yMin() );

    #undef CMD

  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_yMax( int nlhs, mxArray       *plhs[],
           int nrhs, mxArray const *prhs[] ) {

    #define CMD "Spline2DMexWrapper('yMax',OBJ): "

    MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
    MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );

    SplineSurf * ptr = DATA_GET( arg_in_1 );

    setScalarValue( arg_out_0, ptr->yMax() );

    #undef CMD

  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_zMin( int nlhs, mxArray       *plhs[],
           int nrhs, mxArray const *prhs[] ) {

    #define CMD "Spline2DMexWrapper('zMin',OBJ): "

    MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
    MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );

    SplineSurf * ptr = DATA_GET( arg_in_1 );

    setScalarValue( arg_out_0, ptr->zMin() );

    #undef CMD

  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_zMax( int nlhs, mxArray       *plhs[],
           int nrhs, mxArray const *prhs[] ) {

    #define CMD "Spline2DMexWrapper('zMax',OBJ): "

    MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
    MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );

    SplineSurf * ptr = DATA_GET( arg_in_1 );

    setScalarValue( arg_out_0, ptr->zMax() );

    #undef CMD

  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_eval( int nlhs, mxArray       *plhs[],
           int nrhs, mxArray const *prhs[] ) {

    #define CMD "Spline2DMexWrapper('eval',OBJ,x,y): "

    MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
    MEX_ASSERT( nrhs == 4, CMD "expected 4 input, nrhs = " << nrhs );

    SplineSurf * ptr = DATA_GET( arg_in_1 );

    mwSize nx, mx, ny, my;
    real_type const * x = getMatrixPointer(
      arg_in_2, nx, mx, CMD "error in reading `x`"
    );
    real_type const * y = getMatrixPointer(
      arg_in_3, ny, my, CMD "error in reading `y`"
    );
    MEX_ASSERT(
      nx == ny && mx == my,
      CMD "size(x) = " << nx << " x " << ny <<
      " must be equal to size(y) = " << mx << " x " << my
    );
    real_type * z = createMatrixValue( arg_out_0, nx, mx );
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

    MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
    MEX_ASSERT( nrhs == 4, CMD "expected 4 input, nrhs = " << nrhs );

    SplineSurf * ptr = DATA_GET( arg_in_1 );

    mwSize nx, mx, ny, my;
    real_type const * x = getMatrixPointer(
      arg_in_2, nx, mx, CMD "error in reading `x`"
    );
    real_type const * y = getMatrixPointer(
      arg_in_3, ny, my, CMD "error in reading `y`"
    );
    MEX_ASSERT(
      nx == ny && mx == my,
      CMD "size(x) = " << nx << " x " << ny <<
      " must be equal to size(y) = " << mx << " x " << my
    );
    real_type * z = createMatrixValue( arg_out_0, nx, mx );
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

    MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
    MEX_ASSERT( nrhs == 4, CMD "expected 4 input, nrhs = " << nrhs );

    SplineSurf * ptr = DATA_GET( arg_in_1 );

    mwSize nx, mx, ny, my;
    real_type const * x = getMatrixPointer(
      arg_in_2, nx, mx, CMD "error in reading `x`"
    );
    real_type const * y = getMatrixPointer(
      arg_in_3, ny, my, CMD "error in reading `y`"
    );
    MEX_ASSERT(
      nx == ny && mx == my,
      CMD "size(x) = " << nx << " x " << ny <<
      " must be equal to size(y) = " << mx << " x " << my
    );
    real_type * z = createMatrixValue( arg_out_0, nx, mx );
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

    MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
    MEX_ASSERT( nrhs == 4, CMD "expected 4 input, nrhs = " << nrhs );

    SplineSurf * ptr = DATA_GET( arg_in_1 );

    mwSize nx, mx, ny, my;
    real_type const * x = getMatrixPointer(
      arg_in_2, nx, mx, CMD "error in reading `x`"
    );
    real_type const * y = getMatrixPointer(
      arg_in_3, ny, my, CMD "error in reading `y`"
    );
    MEX_ASSERT(
      nx == ny && mx == my,
      CMD "size(x) = " << nx << " x " << ny <<
      " must be equal to size(y) = " << mx << " x " << my
    );
    real_type * z = createMatrixValue( arg_out_0, nx, mx );
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

    MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
    MEX_ASSERT( nrhs == 4, CMD "expected 4 input, nrhs = " << nrhs );

    SplineSurf * ptr = DATA_GET( arg_in_1 );

    mwSize nx, mx, ny, my;
    real_type const * x = getMatrixPointer(
      arg_in_2, nx, mx, CMD "error in reading `x`"
    );
    real_type const * y = getMatrixPointer(
      arg_in_3, ny, my, CMD "error in reading `y`"
    );
    MEX_ASSERT(
      nx == ny && mx == my,
      CMD "size(x) = " << nx << " x " << ny <<
      " must be equal to size(y) = " << mx << " x " << my
    );
    real_type * z = createMatrixValue( arg_out_0, nx, mx );
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

    MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
    MEX_ASSERT( nrhs == 4, CMD "expected 4 input, nrhs = " << nrhs );

    SplineSurf * ptr = DATA_GET( arg_in_1 );

    mwSize nx, mx, ny, my;
    real_type const * x = getMatrixPointer(
      arg_in_2, nx, mx, CMD "error in reading `x`"
    );
    real_type const * y = getMatrixPointer(
      arg_in_3, ny, my, CMD "error in reading `y`"
    );
    MEX_ASSERT(
      nx == ny && mx == my,
      CMD "size(x) = " << nx << " x " << ny <<
      " must be equal to size(y) = " << mx << " x " << my
    );
    real_type * z = createMatrixValue( arg_out_0, nx, mx );
    for ( mwSize j = 0; j < mx*nx; ++j )
      *z++ = ptr->Dyy( *x++, *y++ );
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  typedef enum {
    CMD_NEW,
    CMD_DELETE,
    CMD_BUILD,
    CMD_NX,
    CMD_NY,
    CMD_XMIN,
    CMD_XMAX,
    CMD_YMIN,
    CMD_YMAX,
    CMD_ZMIN,
    CMD_ZMAX,
    CMD_EVAL,
    CMD_EVAL_DX,
    CMD_EVAL_DY,
    CMD_EVAL_DXX,
    CMD_EVAL_DXY,
    CMD_EVAL_DYY
  } CMD_LIST;

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static map<string,unsigned> cmd_to_idx = {
    {"new",CMD_NEW},
    {"delete",CMD_DELETE},
    {"build",CMD_BUILD},
    {"nx",CMD_NX},
    {"ny",CMD_NY},
    {"xMin",CMD_XMIN},
    {"xMax",CMD_XMAX},
    {"yMin",CMD_YMIN},
    {"yMax",CMD_YMAX},
    {"zMin",CMD_ZMIN},
    {"zMax",CMD_ZMAX},
    {"eval",CMD_EVAL},
    {"eval_Dx",CMD_EVAL_DX},
    {"eval_Dy",CMD_EVAL_DY},
    {"eval_Dxx",CMD_EVAL_DXX},
    {"eval_Dxy",CMD_EVAL_DXY},
    {"eval_Dyy",CMD_EVAL_DYY}
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
      case CMD_NEW:
        do_new( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_DELETE:
        do_delete( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_BUILD:
        do_build( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_NX:
        do_nx( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_NY:
        do_ny( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_XMIN:
        do_xMin( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_XMAX:
        do_xMax( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_YMIN:
        do_yMin( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_YMAX:
        do_yMax( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_ZMIN:
        do_zMin( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_ZMAX:
        do_zMax( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_EVAL:
        do_eval( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_EVAL_DX:
        do_eval_Dx( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_EVAL_DY:
        do_eval_Dy( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_EVAL_DXX:
        do_eval_Dxx( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_EVAL_DXY:
        do_eval_Dxy( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_EVAL_DYY:
        do_eval_Dyy( nlhs, plhs, nrhs, prhs );
        break;
      }

    } catch ( exception const & e ) {
      mexErrMsgTxt(e.what());
    } catch (...) {
      mexErrMsgTxt("Spline2DMexWrapper failed\n");
    }

  }


}
