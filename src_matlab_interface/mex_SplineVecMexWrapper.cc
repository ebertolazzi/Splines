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

#define ASSERT(COND,MSG)                           \
  if ( !(COND) ) {                                 \
    std::ostringstream ost;                        \
    ost << "SplineVecMexWrapper: " << MSG << '\n'; \
    mexErrMsgTxt(ost.str().c_str());               \
  }

static
void
mexErrorMessage() {
  // Check for proper number of arguments, etc
  mexErrMsgTxt(
"%======================================================================%\n"
"% SplineVecMexWrapper:  Compute spline curve                           %\n"
"%                                                                      %\n"
"% USAGE:                                                               %\n"
"%   obj = SplineVecMexWrapper( 'new' );                                %\n"
"%   SplineVecMexWrapper( 'delete', obj );                              %\n"
"%   SplineVecMexWrapper( 'setup', obj, Y );                            %\n"
"%   SplineVecMexWrapper( 'knots', obj, X );                            %\n"
"%   SplineVecMexWrapper( 'chord', obj );                               %\n"
"%   SplineVecMexWrapper( 'centripetal', obj );                         %\n"
"%   SplineVecMexWrapper( 'CatmullRom', obj );                          %\n"
"%   P    = SplineVecMexWrapper( 'eval', obj, X );                      %\n"
"%   DP   = SplineVecMexWrapper( 'eval_D', obj, X );                    %\n"
"%   DDP  = SplineVecMexWrapper( 'eval_DD', obj, X );                   %\n"
"%   DDDP = SplineVecMexWrapper( 'eval_DDD', obj, X );                  %\n"
"%                                                                      %\n"
"% On input:                                                            %\n"
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
  SplineVec *
  DATA_NEW( mxArray * & mx_id, SplineVec * ptr ) {
    mx_id = convertPtr2Mat<SplineVec>(ptr);
    return ptr;
  }

  static
  inline
  SplineVec *
  DATA_GET( mxArray const * & mx_id ) {
    return convertMat2Ptr<SplineVec>(mx_id);
  }

  static
  void
  DATA_DELETE( mxArray const * & mx_id ) {
    destroyObject<SplineVec>(mx_id);
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
                  "SplineVecMexWrapper: First argument must be a string" );
      string cmd = mxArrayToString(arg_in_0);

      bool do_new = cmd == "new";
      SplineVec * ptr = nullptr;

      if ( !do_new ) {
        MEX_ASSERT( nrhs >= 1,
                    "SplineVecMexWrapper: expected at least 2 inputs, nrhs = " << nrhs );
        ptr = DATA_GET(arg_in_1);
      }

      if ( do_new ) {

        #define CMD "SplineVecMexWrapper('new'): "

        MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
        MEX_ASSERT( nrhs == 1, CMD "expected 1 input, nrhs = " << nrhs );

        ptr = DATA_NEW( arg_out_0, new SplineVec() );

        #undef CMD

      } else if ( cmd == "setup" ) {
        #define CMD "SplineVecMexWrapper( 'setup', obj, Y ): "

        MEX_ASSERT( nlhs == 0, CMD "expected NO output, nlhs = " << nlhs );
        MEX_ASSERT( nrhs == 3, CMD "expected 3 input, nrhs = " << nrhs );

        mwSize dim, npts;
        real_type const * Y = getMatrixPointer( arg_in_2, dim, npts,
                                                CMD "error in reading 'Y'");
        ptr->setup( dim, npts, Y, dim );

        #undef CMD
      } else if ( cmd == "knots" ) {
        #define CMD "SplineVecMexWrapper( 'setup', obj, X ): "

        MEX_ASSERT( nlhs == 0, CMD "expected NO output, nlhs = " << nlhs );
        MEX_ASSERT( nrhs == 3, CMD "expected 3 input, nrhs = " << nrhs );

        mwSize npts;
        real_type const * X = getVectorPointer( arg_in_2, npts,
                                                CMD "error in reading 'X'");
        MEX_ASSERT( npts == ptr->numPoints(),
                    CMD "size(X) = " << npts <<
                    " must be = " << ptr->dimension() );
        ptr->setKnots( X );

        #undef CMD
      } else if ( cmd == "chord" ) {
        #define CMD "SplineVecMexWrapper( 'chord', obj ): "

        MEX_ASSERT( nlhs == 0, CMD "expected NO output, nlhs = " << nlhs );
        MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );

        ptr->setKnotsChordLength();

        #undef CMD
      } else if ( cmd == "centripetal" ) {
        #define CMD "SplineVecMexWrapper( 'centripetal', obj ): "

        MEX_ASSERT( nlhs == 0, CMD "expected NO output, nlhs = " << nlhs );
        MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );

        ptr->setKnotsCentripetal();

        #undef CMD
      } else if ( cmd == "CatmullRom" ) {
        #define CMD "SplineVecMexWrapper( 'CatmullRom', obj ): "

        MEX_ASSERT( nlhs == 0, CMD "expected NO output, nlhs = " << nlhs );
        MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );

        ptr->CatmullRom();

        #undef CMD

      } else if ( cmd == "eval"    ||
                  cmd == "eval_D"  ||
                  cmd == "eval_DD" ||
                  cmd == "eval_DDD") {

        #define CMD "SplineVecMexWrapper('eval[_*]',OBJ,x): "

        MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
        MEX_ASSERT( nrhs == 3, CMD "expected 3 input, nrhs = " << nrhs );

        mwSize nx;
        real_type const * x = getVectorPointer( arg_in_2, nx,
                                                CMD "error in reading `x`" );

        mwSize dim = ptr->dimension();
        real_type * Y = createMatrixValue( arg_out_0, dim, nx );

        if ( cmd == "eval" ) {
          for ( mwSize i = 0 ; i < nx ; ++i, Y += dim )
            ptr->eval( x[i], Y, 1 );
        } else if ( cmd == "eval_D" ) {
          for ( mwSize i = 0 ; i < nx ; ++i, Y += dim )
            ptr->eval_D( x[i], Y, 1 );
        } else if ( cmd == "eval_DD" ) {
          for ( mwSize i = 0 ; i < nx ; ++i, Y += dim )
            ptr->eval_DD( x[i], Y, 1 );
        } else if ( cmd == "eval_DDD") {
          for ( mwSize i = 0 ; i < nx ; ++i, Y += dim )
            ptr->eval_DDD( x[i], Y, 1 );
        }

        #undef CMD

      } else if ( cmd == "tmin" ) {

        #define CMD "SplineVecMexWrapper('tmin',OBJ): "

        MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
        MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );

        setScalarValue( arg_out_0, ptr->xMin() );

        #undef CMD

      } else if ( cmd == "tmax" ) {

        #define CMD "SplineVecMexWrapper('tmax',OBJ): "

        MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
        MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );

        setScalarValue( arg_out_0, ptr->xMax() );

        #undef CMD

      } else if ( cmd == "delete" ) {

        #define CMD "SplineVecMexWrapper('delete',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs, nrhs = " << nrhs );
        MEX_ASSERT(nlhs == 0, CMD "expected no output, nlhs = " << nlhs );

        // Destroy the C++ object
        DATA_DELETE(arg_in_1);

        // Warn if other commands were ignored

        #undef CMD

      } else {
        MEX_ASSERT( false,
                    "SplineVecMexWrapper('" << cmd <<
                    "',...): Unknown command" );
      }

    } catch ( std::exception const & e ) {
      mexErrMsgTxt(e.what());

    } catch (...) {
      mexErrMsgTxt("SplineVecMexWrapper failed\n");

    }
  }
}
