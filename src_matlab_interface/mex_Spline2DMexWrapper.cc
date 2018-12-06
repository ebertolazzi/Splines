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

static
void
mexErrorMessage() {
  // Check for proper number of arguments, etc
  mexErrMsgTxt(
"%======================================================================%\n"
"% Spline2DMexWrapper:  Compute spline curve                            %\n"
"%                                                                      %\n"
"% USAGE:                                                               %\n"
"%   obj = Spline2DMexWrapper( 'new', kind );                           %\n"
"%   Spline2DMexWrapper( 'delete', obj );                               %\n"
"%   Spline2DMexWrapper( 'build', obj, X, Y, Z );                       %\n"
"%   P    = Spline2DMexWrapper( 'eval', obj, X, Y );                    %\n"
"%   Dx   = Spline2DMexWrapper( 'eval_Dx', obj, X, Y );                 %\n"
"%   Dy   = Spline2DMexWrapper( 'eval_Dy', obj, X, Y );                 %\n"
"%   Dxx  = Spline2DMexWrapper( 'eval_Dxx', obj, X, Y );                %\n"
"%   Dxy  = Spline2DMexWrapper( 'eval_Dxy', obj, X, Y );                %\n"
"%   Dyy  = Spline2DMexWrapper( 'eval_Dyy', obj, X, Y );                %\n"
"%                                                                      %\n"
"% On input:                                                            %\n"
"%                                                                      %\n"
"%  kind = string with the kind of spline, any of:                      %\n"
"%         'bilinear', 'bicubic', 'akima', 'biquintic'                  %\n"
"%  X = vector of X coordinates                                         %\n"
"%  Y = vector of Y coordinates                                         %\n"
"%  Z = matrix of Z coordinates                                         %\n"
"%                                                                      %\n"
"% On output:                                                           %\n"
"%                                                                      %\n"
"%  P    = vector of Z values                                           %\n"
"%                                                                      %\n"
"%======================================================================%\n"
"%                                                                      %\n"
"%  Autor: Enrico Bertolazzi                                            %\n"
"%         Department of Industrial Engineering                         %\n"
"%         University of Trento                                         %\n"
"%         enrico.bertolazzi@unitn.it                                   %\n"
"%                                                                      %\n"
"%======================================================================%\n" );
}

using namespace std;

namespace Splines {

  static
  SplineSurf *
  DATA_NEW( mxArray * & mx_id, SplineSurf * ptr ) {
    mx_id = convertPtr2Mat<SplineSurf>(ptr);
    return ptr;
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

  extern "C"
  void
  mexFunction( int nlhs, mxArray       *plhs[],
               int nrhs, mxArray const *prhs[] ) {

    // the first argument must be a string
    if ( nrhs == 0 ) {
      mexErrorMessage();
      return;
    }

    try {

      MEX_ASSERT( mxIsChar(arg_in_0),
                  "Spline2DMexWrapper: First argument must be a string" );
      string cmd = mxArrayToString(arg_in_0);

      bool do_new = cmd == "new";
      SplineSurf * ptr = nullptr;

      if ( !do_new ) ptr = DATA_GET(arg_in_1);

      if ( do_new ) {

        #define CMD "Spline2DMexWrapper('new',kind): "

        MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
        MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );

        MEX_ASSERT( mxIsChar(arg_in_1),
                    CMD "second argument must be a string, found ``" <<
                    mxGetClassName(arg_in_1) << "''" );

        // the first argument must be a string
        string tname = mxArrayToString(arg_in_1);

        SplineSurf * p = nullptr;

        if      ( tname == "bilinear"  ) p = new Splines::BilinearSpline();
        else if ( tname == "bicubic"   ) p = new Splines::BiCubicSpline();
        else if ( tname == "akima"     ) p = new Splines::Akima2Dspline();
        else if ( tname == "biquintic" ) p = new Splines::BiQuinticSpline();
        else {
         ASSERT( false, "Second argument must be one of the strings:\n" <<
                 "'linear', 'cubic', 'akima', 'bessel', 'pchip', 'quintic'" );
        }

        ptr = DATA_NEW( arg_out_0, p );

        #undef CMD

    } else if ( cmd == "build" ) {

        #define CMD "Spline2DMexWrapper('build',OBJ,x,y,z): "

        MEX_ASSERT( nlhs == 0, CMD "expected no output, nlhs = " << nlhs );
        MEX_ASSERT( nrhs == 5, CMD "expected 5 input, nrhs = " << nrhs );

        mwSize nx, ny, nnx, nny;
        real_type const * x = getVectorPointer( arg_in_2, nx,
                                                CMD "error in reading 'x'");
        real_type const * y = getVectorPointer( arg_in_3, ny,
                                                CMD "error in reading 'y'");
        real_type const * z = getMatrixPointer( arg_in_4, nnx, nny,
                                                CMD "error in reading 'z'");

        MEX_ASSERT( nx == nnx,
                    CMD "lenght of 'x' (" << nx <<
                    ") must be the number of row of 'z' (" << nnx << ")" );

        MEX_ASSERT( ny == nny,
                    CMD "lenght of 'y' (" << ny <<
                    ") must be the number of column of 'z' (" << nny << ")" );

        bool fortran_storage = true;
        bool transposed      = false;

        integer ldZ = nx;
        ptr -> build ( x, 1, y, 1, z, ldZ, nx, ny,
                       fortran_storage, transposed );

        #undef CMD

      } else if ( cmd == "nx" ) {

        #define CMD "Spline2DMexWrapper('nx',OBJ): "

        MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
        MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );

        setScalarInt( arg_out_0, ptr->numPointX() );

        #undef CMD

      } else if ( cmd == "ny" ) {

        #define CMD "Spline2DMexWrapper('nx',OBJ): "

        MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
        MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );

        setScalarInt( arg_out_0, ptr->numPointY() );

        #undef CMD

      } else if ( cmd == "xMin" ) {

        #define CMD "Spline2DMexWrapper('xMin',OBJ): "

        MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
        MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );

        setScalarValue( arg_out_0, ptr->xMin() );

        #undef CMD

      } else if ( cmd == "xMax" ) {

        #define CMD "Spline2DMexWrapper('xMax',OBJ): "

        MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
        MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );

        setScalarValue( arg_out_0, ptr->xMax() );

        #undef CMD

      } else if ( cmd == "yMin" ) {

        #define CMD "Spline2DMexWrapper('yMin',OBJ): "

        MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
        MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );

        setScalarValue( arg_out_0, ptr->yMin() );

        #undef CMD

      } else if ( cmd == "yMax" ) {

        #define CMD "Spline2DMexWrapper('yMax',OBJ): "

        MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
        MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );

        setScalarValue( arg_out_0, ptr->yMax() );

        #undef CMD

      } else if ( cmd == "zMin" ) {

        #define CMD "Spline2DMexWrapper('zMin',OBJ): "

        MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
        MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );

        setScalarValue( arg_out_0, ptr->zMin() );

        #undef CMD

      } else if ( cmd == "zMax" ) {

        #define CMD "Spline2DMexWrapper('zMax',OBJ): "

        MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
        MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );

        setScalarValue( arg_out_0, ptr->zMax() );

        #undef CMD

      } else if ( cmd == "eval" ) {

#define EVAL_LOOP(EVAL) \
MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs ); \
MEX_ASSERT( nrhs == 4, CMD "expected 4 input, nrhs = " << nrhs ); \
\
mwSize nx, mx, ny, my; \
real_type const * x = getMatrixPointer( arg_in_2, nx, mx, \
                                        CMD "error in reading `x`" ); \
real_type const * y = getMatrixPointer( arg_in_3, ny, my, \
                                        CMD "error in reading `y`" ); \
\
MEX_ASSERT( nx == ny && mx == my, \
            CMD "size(x) = " << nx << " x " << ny << \
            " must be equal to size(y) = " << mx << " x " << my ); \
\
real_type * z = createMatrixValue( arg_out_0, nx, mx ); \
\
for ( mwSize j = 0; j < mx; ++j ) \
  for ( mwSize i = 0; i < nx; ++i ) \
    *z++ = EVAL( *x++, *y++ );

        #define CMD "Spline2DMexWrapper('eval',OBJ,x,y): "
        EVAL_LOOP( (*ptr) );
        #undef CMD

      } else if ( cmd == "eval_Dx" ) {

        #define CMD "Spline2DMexWrapper('eval_Dx',OBJ,x,y): "
        EVAL_LOOP( ptr->Dx );
        #undef CMD

      } else if ( cmd == "eval_Dy" ) {

        #define CMD "Spline2DMexWrapper('eval_Dy',OBJ,x,y): "
        EVAL_LOOP( ptr->Dy );
        #undef CMD

      } else if ( cmd == "eval_Dxx" ) {

        #define CMD "Spline2DMexWrapper('eval_Dxx',OBJ,x,y): "
        EVAL_LOOP( ptr->Dxx );
        #undef CMD

      } else if ( cmd == "eval_Dxy" ) {

        #define CMD "Spline2DMexWrapper('eval_Dxy',OBJ,x,y): "
        EVAL_LOOP( ptr->Dxy );
        #undef CMD

      } else if ( cmd == "eval_Dyy" ) {

        #define CMD "Spline2DMexWrapper('eval_Dyy',OBJ,x,y): "
        EVAL_LOOP( ptr->Dyy );
        #undef CMD

      } else if ( cmd == "delete" ) {

        #define CMD "Spline2DMexWrapper('delete',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs, nrhs = " << nrhs );
        MEX_ASSERT(nlhs == 0, CMD "expected no output, nlhs = " << nlhs );

        // Destroy the C++ object
        DATA_DELETE(arg_in_1);

        // Warn if other commands were ignored

        #undef CMD

      } else {
        MEX_ASSERT( false,
                    "Spline2DMexWrapper('" << cmd <<
                    "',...): Unknown command" );
      }

    } catch ( std::exception const & e ) {
      mexErrMsgTxt(e.what());

    } catch (...) {
      mexErrMsgTxt("Spline2DMexWrapper failed\n");

    }
  }
}
