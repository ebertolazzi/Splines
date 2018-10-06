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
    ost << "SplineSetMexWrapper: " << MSG << '\n'; \
    mexErrMsgTxt(ost.str().c_str());               \
  }

static
void
mexErrorMessage() {
  // Check for proper number of arguments, etc
  mexErrMsgTxt(
"%======================================================================%\n"
"% SplineSetMexWrapper:  Compute spline curve                           %\n"
"%                                                                      %\n"
"% USAGE:                                                               %\n"
"%   obj = SplineSetMexWrapper( 'new' );                                %\n"
"%   SplineSetMexWrapper( 'delete', obj );                              %\n"
"%   SplineSetMexWrapper( 'build', obj, kinds, X, Y );                  %\n"
"%   P    = SplineSetMexWrapper( 'eval', obj, X );                      %\n"
"%   DP   = SplineSetMexWrapper( 'eval_D', obj, X );                    %\n"
"%   DDP  = SplineSetMexWrapper( 'eval_DD', obj, X );                   %\n"
"%   DDDP = SplineSetMexWrapper( 'eval_DDD', obj, X );                  %\n"
"%                                                                      %\n"
"% On input:                                                            %\n"
"%                                                                      %\n"
"%  kinds = string cell array with the kind of spline, any of:          %\n"
"%          'linear', 'cubic', 'akima', 'bessel', 'pchip', 'quintic'    %\n"
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
  SplineSet *
  DATA_NEW( mxArray * & mx_id, SplineSet * ptr ) {
    mx_id = convertPtr2Mat<SplineSet>(ptr);
    return ptr;
  }

  static
  inline
  SplineSet *
  DATA_GET( mxArray const * & mx_id ) {
    return convertMat2Ptr<SplineSet>(mx_id);
  }

  static
  void
  DATA_DELETE( mxArray const * & mx_id ) {
    destroyObject<SplineSet>(mx_id);
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
      SplineSet * ptr = nullptr;

      if ( !do_new ) {
        MEX_ASSERT( nrhs >= 1,
                    "SplineVecMexWrapper: expected at least 2 inputs, nrhs = " << nrhs );
        ptr = DATA_GET(arg_in_1);
      }

      if ( do_new ) {

        #define CMD "SplineVecMexWrapper('new'): "

        MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
        MEX_ASSERT( nrhs == 1, CMD "expected 1 input, nrhs = " << nrhs );

        ptr = DATA_NEW( arg_out_0, new SplineSet() );

        #undef CMD

      } else if ( cmd == "build" ) {

        #define CMD "SplineVecMexWrapper('build',OBJ,kinds,X,Y): "

        MEX_ASSERT( nlhs == 0, CMD "expected NO output, nlhs = " << nlhs );
        MEX_ASSERT( nrhs == 5, CMD "expected 5 input, nrhs = " << nrhs );

        mwSize n, nn, nspl;
        real_type const * x = getVectorPointer( arg_in_3, n,
                                                CMD "error in reading 'x'");
        real_type const * y = getMatrixPointer( arg_in_4, nn, nspl,
                                                CMD "error in reading 'y'");

        MEX_ASSERT( n == nn,
                    CMD "lenght of 'x' must be the numnber of rows of 'y'" );

        std::vector<SplineType> types;
        types.reserve(nspl);

        if ( mxIsChar(arg_in_2) ) {
          string tname = mxArrayToString(arg_in_2);
          Splines::SplineType st;
          if      ( tname == "linear"  ) st = LINEAR_TYPE;
          else if ( tname == "cubic"   ) st = CUBIC_TYPE;
          else if ( tname == "akima"   ) st = AKIMA_TYPE;
          else if ( tname == "bessel"  ) st = BESSEL_TYPE;
          else if ( tname == "pchip"   ) st = PCHIP_TYPE;
          else if ( tname == "quintic" ) st = QUINTIC_TYPE;
          else {
            MEX_ASSERT( false,
                        CMD "Cell array of strings must contains the strings:\n" <<
                        "'linear', 'cubic', 'akima', 'bessel', 'pchip', 'quintic'" );
          }
          for ( mwSize i = 0; i < nspl; ++i ) types.push_back(st);

        } else if ( mxIsCell(arg_in_2) ) {

          mwSize const *dims = mxGetDimensions(arg_in_2);
          MEX_ASSERT( dims[0] == nspl && dims[1] == 1,
                      CMD "Third argument expected to be cell array " << nspl <<
                      " x 1, found " << dims[0] << " x " << dims[1] );

          for ( mwSize i = 0; i < nspl; ++i ) {
            mxArray const * cell = mxGetCell(arg_in_2,i);
            MEX_ASSERT( mxIsChar(cell),
                        CMD "Third argument expected to be cell array of strings" );
            string tname = mxArrayToString(cell);
            Splines::SplineType st;
            if      ( tname == "linear"  ) st = LINEAR_TYPE;
            else if ( tname == "cubic"   ) st = CUBIC_TYPE;
            else if ( tname == "akima"   ) st = AKIMA_TYPE;
            else if ( tname == "bessel"  ) st = BESSEL_TYPE;
            else if ( tname == "pchip"   ) st = PCHIP_TYPE;
            else if ( tname == "quintic" ) st = QUINTIC_TYPE;
            else {
              MEX_ASSERT( false,
                          CMD "Cell array of strings must contains the strings:\n" <<
                          "'linear', 'cubic', 'akima', 'bessel', 'pchip', 'quintic'" );
            }
            types.push_back(st);
          }
        } else {
          MEX_ASSERT( false,
                      CMD "Third argument expected to be a string of a cell array of strings, found ``" <<
                      mxGetClassName(arg_in_2) << "''" );
        }

        std::vector<char const *>      headers;
        std::vector<real_type const *> pY;
        headers.reserve(nspl);
        pY.reserve(nspl);
        for ( mwSize i = 0; i < nspl; ++i ) {
          headers.push_back("noname");
          pY.push_back(y+i*nn);
        }
        ptr -> build( nspl, n,
                      &headers.front(), &types.front(),
                      x, &pY.front(), nullptr );

        #undef CMD

      } else if ( cmd == "eval"    ||
                  cmd == "eval_D"  ||
                  cmd == "eval_DD" ||
                  cmd == "eval_DDD") {

        #define CMD "SplineSetMexWrapper('eval[_*]',OBJ,x): "

        MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
        MEX_ASSERT( nrhs == 3, CMD "expected 3 input, nrhs = " << nrhs );

        mwSize nx;
        real_type const * x = getVectorPointer( arg_in_2, nx,
                                                CMD "error in reading `x`" );

        mwSize dim = ptr->numSplines();
        real_type * Y = createMatrixValue( arg_out_0, dim, nx );

        for ( mwSize nsp = 0 ; nsp < dim; ++nsp ) {
          Spline const * S = ptr->getSpline( nsp );
          real_type * y = Y+nsp;
          if ( cmd == "eval" ) {
            for ( mwSize i = 0 ; i < nx ; ++i, y += dim ) *y = S->eval( x[i] );
          } else if ( cmd == "eval_D" ) {
            for ( mwSize i = 0 ; i < nx ; ++i, y += dim ) *y = S->eval_D( x[i] );
          } else if ( cmd == "eval_DD" ) {
            for ( mwSize i = 0 ; i < nx ; ++i, y += dim ) *y = S->eval_DD( x[i] );
          } else if ( cmd == "eval_DDD") {
            for ( mwSize i = 0 ; i < nx ; ++i, y += dim ) *y = S->eval_DDD( x[i] );
          }
        }

        #undef CMD

      } else if ( cmd == "tmin" ) {

        #define CMD "SplineSetMexWrapper('tmin',OBJ): "

        MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
        MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );

        setScalarValue( arg_out_0, ptr->xMin() );

        #undef CMD

      } else if ( cmd == "tmax" ) {

        #define CMD "SplineSetMexWrapper('tmax',OBJ): "

        MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
        MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );

        setScalarValue( arg_out_0, ptr->xMax() );

        #undef CMD

      } else if ( cmd == "delete" ) {

        #define CMD "SplineSetMexWrapper('delete',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs, nrhs = " << nrhs );
        MEX_ASSERT(nlhs == 0, CMD "expected no output, nlhs = " << nlhs );

        // Destroy the C++ object
        DATA_DELETE(arg_in_1);

        // Warn if other commands were ignored

        #undef CMD

      } else {
        MEX_ASSERT( false,
                    "SplineSetMexWrapper('" << cmd <<
                    "',...): Unknown command" );
      }

    } catch ( std::exception const & e ) {
      mexErrMsgTxt(e.what());

    } catch (...) {
      mexErrMsgTxt("SplineSetMexWrapper failed\n");

    }
  }
}
