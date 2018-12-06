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
    ost << "Spline1DMexWrapper: " << MSG << '\n'; \
    mexErrMsgTxt(ost.str().c_str());              \
  }

static
void
mexErrorMessage() {
  // Check for proper number of arguments, etc
  mexErrMsgTxt(
"%======================================================================%\n"
"% Spline1DMexWrapper:  Compute spline curve                            %\n"
"%                                                                      %\n"
"% USAGE:                                                               %\n"
"%   obj = Spline1DMexWrapper( 'new', kind );                           %\n"
"%   Spline1DMexWrapper( 'delete', obj );                               %\n"
"%   Spline1DMexWrapper( 'build', obj, X, Y );                          %\n"
"%   P    = Spline1DMexWrapper( 'eval', obj, X );                       %\n"
"%   DP   = Spline1DMexWrapper( 'eval_D', obj, X );                     %\n"
"%   DDP  = Spline1DMexWrapper( 'eval_DD', obj, X );                    %\n"
"%   DDDP = Spline1DMexWrapper( 'eval_DDD', obj, X );                   %\n"
"%                                                                      %\n"
"% On input:                                                            %\n"
"%                                                                      %\n"
"%  kind = string with the kind of spline, any of:                      %\n"
"%         'linear', 'cubic', 'akima', 'bessel', 'pchip', 'quintic'     %\n"
"%  X = vector of X coordinates                                         %\n"
"%  Y = vector of Y coordinates                                         %\n"
"%                                                                      %\n"
"% On output:                                                           %\n"
"%                                                                      %\n"
"%  P    = vector of Y values                                           %\n"
"%  DP   = vector of dimension size(X) with derivative                  %\n"
"%  DDP  = vector of dimension size(X) with second derivative           %\n"
"%  DDDP = vector of dimension size(X) with third derivative            %\n"
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
  Spline *
  DATA_NEW( mxArray * & mx_id, Spline * ptr ) {
    mx_id = convertPtr2Mat<Spline>(ptr);
    return ptr;
  }

  static
  inline
  Spline *
  DATA_GET( mxArray const * & mx_id ) {
    return convertMat2Ptr<Spline>(mx_id);
  }

  static
  void
  DATA_DELETE( mxArray const * & mx_id ) {
    destroyObject<Spline>(mx_id);
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
                  "Spline1DMexWrapper: First argument must be a string" );
      string cmd = mxArrayToString(arg_in_0);

      bool do_new = cmd == "new";
      Spline * ptr = nullptr;

      if ( !do_new ) ptr = DATA_GET(arg_in_1);

      if ( do_new ) {

        #define CMD "Spline1DMexWrapper('new',kind): "

        MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
        MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );

        MEX_ASSERT( mxIsChar(arg_in_1),
                    CMD "second argument must be a string, found ``" <<
                    mxGetClassName(arg_in_1) << "''" );

        // the first argument must be a string
        string tname = mxArrayToString(arg_in_1);

        Spline * p = nullptr;

        if      ( tname == "linear"  ) p = new Splines::LinearSpline();
        else if ( tname == "cubic"   ) p = new Splines::CubicSpline();
        else if ( tname == "akima"   ) p = new Splines::AkimaSpline();
        else if ( tname == "bessel"  ) p = new Splines::BesselSpline();
        else if ( tname == "pchip"   ) p = new Splines::PchipSpline();
        else if ( tname == "quintic" ) p = new Splines::QuinticSpline();
        else {
         ASSERT( false, "Second argument must be one of the strings:\n" <<
                 "'linear', 'cubic', 'akima', 'bessel', 'pchip', 'quintic'" );
        }

        ptr = DATA_NEW( arg_out_0, p );

        #undef CMD

    } else if ( cmd == "build" ) {

        #define CMD "Spline1DMexWrapper('build',OBJ,x,y): "

        MEX_ASSERT( nlhs == 0, CMD "expected no output, nlhs = " << nlhs );
        MEX_ASSERT( nrhs == 4, CMD "expected 4 input, nrhs = " << nrhs );

        mwSize nx, ny;
        real_type const * x = getVectorPointer( arg_in_2, nx,
                                                CMD "error in reading 'x'");
        real_type const * y = getVectorPointer( arg_in_3, ny,
                                                CMD "error in reading 'y'");

        MEX_ASSERT( nx == ny,
                    "lenght of 'x' must be the lenght of 'y'" );

        ptr -> build( x, y, nx );

        #undef CMD

      } else if ( cmd == "degree" ) {

        #define CMD "Spline1DMexWrapper('degree',OBJ): "

        MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
        MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );

        setScalarInt( arg_out_0, ptr->order()-1 );

        #undef CMD

      } else if ( cmd == "order" ) {

        #define CMD "Spline1DMexWrapper('order',OBJ): "

        MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
        MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );

        setScalarInt( arg_out_0, ptr->order() );

        #undef CMD

      } else if ( cmd == "eval" ) {

        #define CMD "Spline1DMexWrapper('eval',OBJ,x): "

        MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
        MEX_ASSERT( nrhs == 3, CMD "expected 3 input, nrhs = " << nrhs );

        mwSize nx;
        real_type const * x = getVectorPointer( arg_in_2, nx,
                                                CMD "error in reading `x`" );
        real_type * y = createMatrixValue( arg_out_0, nx, 1 );

        for ( mwSize i = 0; i < nx; ++i ) *y++ = ptr->eval( *x++ );


        #undef CMD

      } else if ( cmd == "eval_D" ) {

        #define CMD "Spline1DMexWrapper('eval_D',OBJ,x): "

        MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
        MEX_ASSERT( nrhs == 3, CMD "expected 3 input, nrhs = " << nrhs );

        mwSize nx;
        real_type const * x = getVectorPointer( arg_in_2, nx,
                                                CMD "error in reading `x`" );
        real_type * y = createMatrixValue( arg_out_0, nx, 1 );

        for ( mwSize i = 0; i < nx; ++i ) *y++ = ptr->eval_D( *x++ );


        #undef CMD

      } else if ( cmd == "eval_DD" ) {

        #define CMD "Spline1DMexWrapper('eval_DD',OBJ,x): "

        MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
        MEX_ASSERT( nrhs == 3, CMD "expected 3 input, nrhs = " << nrhs );

        mwSize nx;
        real_type const * x = getVectorPointer( arg_in_2, nx,
                                                CMD "error in reading `x`" );
        real_type * y = createMatrixValue( arg_out_0, nx, 1 );

        for ( mwSize i = 0; i < nx; ++i ) *y++ = ptr->eval_DD( *x++ );


        #undef CMD

      } else if ( cmd == "eval_DDD" ) {

        #define CMD "Spline1DMexWrapper('eval_DDD',OBJ,x): "

        MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
        MEX_ASSERT( nrhs == 3, CMD "expected 3 input, nrhs = " << nrhs );

        mwSize nx;
        real_type const * x = getVectorPointer( arg_in_2, nx,
                                                CMD "error in reading `x`" );
        real_type * y = createMatrixValue( arg_out_0, nx, 1 );

        for ( mwSize i = 0; i < nx; ++i ) *y++ = ptr->eval_DDD( *x++ );


        #undef CMD

      } else if ( cmd == "delete" ) {

        #define CMD "Spline1DMexWrapper('delete',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs, nrhs = " << nrhs );
        MEX_ASSERT(nlhs == 0, CMD "expected no output, nlhs = " << nlhs );

        // Destroy the C++ object
        DATA_DELETE(arg_in_1);

        // Warn if other commands were ignored

        #undef CMD

      } else {
        MEX_ASSERT( false,
                    "Spline1DMexWrapper('" << cmd <<
                    "',...): Unknown command" );
      }

    } catch ( std::exception const & e ) {
      mexErrMsgTxt(e.what());

    } catch (...) {
      mexErrMsgTxt("Spline1DMexWrapper failed\n");

    }
  }
}
