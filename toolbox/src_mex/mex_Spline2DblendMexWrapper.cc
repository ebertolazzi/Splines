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

#include "GenericContainer/GenericContainerInterface_matlab.hh"

#define ASSERT(COND,MSG)                          \
  if ( !(COND) ) {                                \
    std::ostringstream ost;                       \
    ost << "Spline2DMexWrapper: " << MSG << '\n'; \
    mexErrMsgTxt(ost.str().c_str());              \
  }

#ifdef __clang__
  #pragma clang diagnostic ignored "-Wexit-time-destructors"
  #pragma clang diagnostic ignored "-Wsign-conversion"
#endif

#include <unordered_map>

namespace Splines {

  using namespace std;

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_new( int nlhs, mxArray       *plhs[],
          int nrhs, mxArray const *prhs[] ) {

    #define MEX_ERROR_MESSAGE_1 "Spline2DblendMexWrapper( 'new', struct )"
    #define CMD MEX_ERROR_MESSAGE_1

    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );

    GenericContainer gc;
    mxArray_to_GenericContainer( arg_in_1, gc );
    
    Spline2Dblend * ptr{ new Splines::Spline2Dblend("2D blend") };

    std::ostringstream ss;
    gc.print(ss);

    try {
      ptr->build( gc );
    } catch ( exception const & e ) {
      mexErrMsgTxt( fmt::format( "Spline2DblendMexWrapper Error: {}", e.what() ).c_str() );
    } catch (...) {
      mexErrMsgTxt("Spline2DblendMexWrapper failed\n");
    }

    arg_out_0 = Utils::mex_convert_ptr_to_mx<Spline2Dblend>( ptr );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_delete( int nlhs, mxArray       *[],
             int nrhs, mxArray const *prhs[] ) {

    #define MEX_ERROR_MESSAGE_2 "Spline2DblendMexWrapper( 'delete', OBJ )"
    #define CMD MEX_ERROR_MESSAGE_2

    UTILS_MEX_ASSERT( nlhs == 0, CMD ": expected 0 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 inputs, nrhs = {}\n", nrhs );

    // Destroy the C++ object
    Utils::mex_destroy_object<Spline2Dblend>(arg_in_1);
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_build( int nlhs, mxArray       *[],
            int nrhs, mxArray const *prhs[] ) {

    #define MEX_ERROR_MESSAGE_3 "Spline2DMexWrapper('build',OBJ,struct)"
    #define CMD MEX_ERROR_MESSAGE_3

    UTILS_MEX_ASSERT( nlhs == 0, CMD ": expected no output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 3, CMD ": expected 3 input, nrhs = {}\n", nrhs );
    
    Spline2Dblend * ptr{ Utils::mex_convert_mx_to_ptr<Spline2Dblend>( arg_in_1 ) };
    
    GenericContainer gc;
    mxArray_to_GenericContainer( arg_in_1, gc );

    try {
      ptr->build( gc );
    } catch ( exception const & e ) {
      mexErrMsgTxt( fmt::format( "Spline2DblendMexWrapper Error: {}", e.what() ).c_str() );
    } catch (...) {
      mexErrMsgTxt("Spline2DblendMexWrapper failed\n");
    }

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_x_min( int nlhs, mxArray       *plhs[],
            int nrhs, mxArray const *prhs[] ) {

    #define MEX_ERROR_MESSAGE_4 "Spline2DMexWrapper('x_min',OBJ)"
    #define CMD MEX_ERROR_MESSAGE_4

    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 input, nrhs = {}\n", nrhs );

    Spline2Dblend * ptr{ Utils::mex_convert_mx_to_ptr<Spline2Dblend>( arg_in_1 ) };

    Utils::mex_set_scalar_value( arg_out_0, ptr->x_min() );

    #undef CMD

  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_x_max( int nlhs, mxArray       *plhs[],
            int nrhs, mxArray const *prhs[] ) {

    #define MEX_ERROR_MESSAGE_5 "Spline2DMexWrapper('x_max',OBJ)"
    #define CMD MEX_ERROR_MESSAGE_5

    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 input, nrhs = {}\n", nrhs );

    Spline2Dblend * ptr{ Utils::mex_convert_mx_to_ptr<Spline2Dblend>( arg_in_1 ) };

    Utils::mex_set_scalar_value( arg_out_0, ptr->x_max() );

    #undef CMD

  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_y_min( int nlhs, mxArray       *plhs[],
            int nrhs, mxArray const *prhs[] ) {

    #define MEX_ERROR_MESSAGE_6 "Spline2DMexWrapper('y_min',OBJ)"
    #define CMD MEX_ERROR_MESSAGE_6

    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 input, nrhs = {}\n", nrhs );

    Spline2Dblend * ptr{ Utils::mex_convert_mx_to_ptr<Spline2Dblend>( arg_in_1 ) };

    Utils::mex_set_scalar_value( arg_out_0, ptr->y_min() );

    #undef CMD

  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_y_max( int nlhs, mxArray       *plhs[],
            int nrhs, mxArray const *prhs[] ) {

    #define MEX_ERROR_MESSAGE_7 "Spline2DMexWrapper('y_max',OBJ)"
    #define CMD MEX_ERROR_MESSAGE_7

    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 input, nrhs = {}\n", nrhs );

    Spline2Dblend * ptr{ Utils::mex_convert_mx_to_ptr<Spline2Dblend>( arg_in_1 ) };

    Utils::mex_set_scalar_value( arg_out_0, ptr->y_max() );

    #undef CMD

  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  
  #define EVALUATE_LOOP( CMD, OP ) \
    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs ); \
    UTILS_MEX_ASSERT( nrhs == 5, CMD ": expected 4 input, nrhs = {}\n", nrhs );  \
    \
    Spline2Dblend * ptr{ Utils::mex_convert_mx_to_ptr<Spline2Dblend>( arg_in_1 ) }; \
    \
    mwSize nx, mx, ny, my; \
    real_type const * x{ Utils::mex_matrix_pointer( arg_in_2, nx, mx, CMD ": error in reading `x`" ) }; \
    real_type const * y{ Utils::mex_matrix_pointer( arg_in_3, ny, my, CMD ": error in reading `y`" ) }; \
    real_type const   s{ Utils::mex_get_scalar_value( arg_in_4, CMD ": error in reading `s`" ) }; \
    UTILS_MEX_ASSERT( \
      nx == ny && mx == my, \
      CMD ": size(x) = {} x {} must be equal to size(y) = {} x {}", \
      nx, ny, mx, my \
    ); \
    real_type * z{ Utils::mex_create_matrix_value( arg_out_0, nx, mx ) }; \
    mwSize nnn{ nx * mx }; \
    for ( mwSize j{0}; j < nnn; ++j ) z[j] = OP( x[j], y[j], s )

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_eval( int nlhs, mxArray       *plhs[],
           int nrhs, mxArray const *prhs[] ) {
    #define MEX_ERROR_MESSAGE_8 "Spline2DMexWrapper('eval',OBJ,x,y,s)"
    EVALUATE_LOOP( MEX_ERROR_MESSAGE_8, ptr->eval );
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_eval_Dx( int nlhs, mxArray       *plhs[],
              int nrhs, mxArray const *prhs[] ) {
    #define MEX_ERROR_MESSAGE_9 "Spline2DMexWrapper('eval_Dx',OBJ,x,y,s)"
    EVALUATE_LOOP( MEX_ERROR_MESSAGE_9, ptr->Dx );
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_eval_Dy( int nlhs, mxArray       *plhs[],
              int nrhs, mxArray const *prhs[] ) {
    #define MEX_ERROR_MESSAGE_10 "Spline2DMexWrapper('eval_Dy',OBJ,x,y,s)"
    EVALUATE_LOOP( MEX_ERROR_MESSAGE_10, ptr->Dy );
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_eval_Dxx( int nlhs, mxArray       *plhs[],
               int nrhs, mxArray const *prhs[] ) {
    #define MEX_ERROR_MESSAGE_11 "Spline2DMexWrapper('eval_Dxx',OBJ,x,y,s)"
    EVALUATE_LOOP( MEX_ERROR_MESSAGE_11, ptr->Dxx );
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_eval_Dxy( int nlhs, mxArray       *plhs[],
               int nrhs, mxArray const *prhs[] ) {
    #define MEX_ERROR_MESSAGE_12 "Spline2DMexWrapper('eval_Dxy',OBJ,x,y,s)"
    EVALUATE_LOOP( MEX_ERROR_MESSAGE_12, ptr->Dxy );
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_eval_Dyy( int nlhs, mxArray       *plhs[],
               int nrhs, mxArray const *prhs[] ) {
    #define MEX_ERROR_MESSAGE_13 "Spline2DMexWrapper('eval_Dyy',OBJ,x,y,s)"
    EVALUATE_LOOP( MEX_ERROR_MESSAGE_13, ptr->Dyy );
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  
  typedef void (*DO_CMD)( int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[] );

  static unordered_map<string,DO_CMD> cmd_to_fun = {
    {"new",do_new},
    {"delete",do_delete},
    {"build",do_build},
    {"x_min",do_x_min},
    {"x_max",do_x_max},
    {"y_min",do_y_min},
    {"y_max",do_y_max},
    {"eval",do_eval},
    {"eval_Dx",do_eval_Dx},
    {"eval_Dy",do_eval_Dy},
    {"eval_Dxx",do_eval_Dxx},
    {"eval_Dxy",do_eval_Dxy},
    {"eval_Dyy",do_eval_Dyy}
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
