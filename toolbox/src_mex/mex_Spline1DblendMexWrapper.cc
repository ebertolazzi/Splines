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

#define ASSERT(COND,MSG)                               \
  if ( !(COND) ) {                                     \
    std::ostringstream ost;                            \
    ost << "Spline1DblendMexWrapper: " << MSG << '\n'; \
    mexErrMsgTxt(ost.str().c_str());                   \
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
  do_new(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define MEX_ERROR_MESSAGE_1 "Spline1DblendMexWrapper( 'new', struct )"
    #define CMD MEX_ERROR_MESSAGE_1

    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );

    GenericContainer gc;
    mxArray_to_GenericContainer( arg_in_1, gc );
    
    Spline1Dblend * ptr{ new Splines::Spline1Dblend("1D blend") };

    std::ostringstream ss;
    gc.print(ss);

    try {
      ptr->build( gc );
    } catch ( exception const & e ) {
      mexErrMsgTxt( fmt::format( "Spline1DblendMexWrapper Error: {}", e.what() ).c_str() );
    } catch (...) {
      mexErrMsgTxt("Spline1DblendMexWrapper failed\n");
    }

    arg_out_0 = Utils::mex_convert_ptr_to_mx<Spline1Dblend>( ptr );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_delete(
    int nlhs, mxArray       *[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_2 "Spline1DblendMexWrapper( 'delete', OBJ )"
    #define CMD MEX_ERROR_MESSAGE_2

    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 0, CMD ": expected 0 output, nlhs = {}\n", nlhs );

    // Destroy the C++ object
    Utils::mex_destroy_object<Spline1Dblend>(arg_in_1);
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_build(
    int nlhs, mxArray       *[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define MEX_ERROR_MESSAGE_3 "Spline1DblendMexWrapper( 'build', OBJ, struct )"
    #define CMD MEX_ERROR_MESSAGE_3

    UTILS_MEX_ASSERT( nlhs == 0, CMD ": expected no output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 3, CMD ": expected 3 input, nrhs = {}\n", nrhs );

    Spline1Dblend * ptr{ Utils::mex_convert_mx_to_ptr<Spline1Dblend>( arg_in_1 ) };

    GenericContainer gc;
    mxArray_to_GenericContainer( arg_in_1, gc );

    try {
      ptr->build( gc );
    } catch ( exception const & e ) {
      mexErrMsgTxt( fmt::format( "Spline1DblendMexWrapper Error: {}", e.what() ).c_str() );
    } catch (...) {
      mexErrMsgTxt("Spline1DblendMexWrapper failed\n");
    }

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_degree(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define MEX_ERROR_MESSAGE_4 "degree = Spline1DblendMexWrapper('degree',OBJ)"
    #define CMD MEX_ERROR_MESSAGE_4

    UTILS_MEX_ASSERT( nlhs == 2, CMD ": expected 2 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 input, nrhs = {}\n", nrhs );
    Spline1Dblend const * ptr{ Utils::mex_convert_mx_to_ptr<Spline1Dblend>( arg_in_1 ) };
    Utils::mex_set_scalar_int32( arg_out_0, ptr->order0()-1 );
    Utils::mex_set_scalar_int32( arg_out_1, ptr->order1()-1 );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_order(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define MEX_ERROR_MESSAGE_5 "order = Spline1DblendMexWrapper('order',OBJ)"
    #define CMD MEX_ERROR_MESSAGE_5

    UTILS_MEX_ASSERT( nlhs == 2, CMD ": expected 2 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 input, nrhs = {}\n", nrhs );
    Spline1Dblend const * ptr{ Utils::mex_convert_mx_to_ptr<Spline1Dblend>( arg_in_1 ) };
    Utils::mex_set_scalar_int32( arg_out_0, ptr->order0() );
    Utils::mex_set_scalar_int32( arg_out_1, ptr->order1() );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_eval(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define MEX_ERROR_MESSAGE_6 "y = Spline1DblendMexWrapper('eval',OBJ,x,s)"
    #define CMD MEX_ERROR_MESSAGE_6

    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 4, CMD ": expected 4 input, nrhs = {}\n", nrhs );

    Spline1Dblend const * ptr{ Utils::mex_convert_mx_to_ptr<Spline1Dblend>( arg_in_1 ) };

    mwSize nx;
    real_type const * x{ Utils::mex_vector_pointer( arg_in_2, nx, CMD ": error in reading `x`" ) };
    real_type const   s{ Utils::mex_get_scalar_value( arg_in_3, CMD ": error in reading `s`" ) };
    real_type       * y{ Utils::mex_create_matrix_value( arg_out_0, nx, 1 ) };

    for ( mwSize i{0}; i < nx; ++i ) *y++ = ptr->eval( *x++, s );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_eval_D(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define MEX_ERROR_MESSAGE_7 "y'(x) = Spline1DblendMexWrapper('eval_D',OBJ,x,s)"
    #define CMD MEX_ERROR_MESSAGE_7

    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 3, CMD ": expected 3 input, nrhs = {}\n", nrhs );

    Spline1Dblend const * ptr{ Utils::mex_convert_mx_to_ptr<Spline1Dblend>( arg_in_1 ) };

    mwSize nx;
    real_type const * x{ Utils::mex_vector_pointer( arg_in_2, nx, CMD ": error in reading `x`"  ) };
    real_type const   s{ Utils::mex_get_scalar_value( arg_in_3, CMD ": error in reading `s`" ) };
    real_type       * y{ Utils::mex_create_matrix_value( arg_out_0, nx, 1 ) };

    for ( mwSize i{0}; i < nx; ++i ) *y++ = ptr->eval_D( *x++, s );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_eval_DD(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define MEX_ERROR_MESSAGE_8 "y''(x) = Spline1DblendMexWrapper('eval_DD',OBJ,x,s): "
    #define CMD MEX_ERROR_MESSAGE_8

    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 3, CMD ": expected 3 input, nrhs = {}\n", nrhs );

    Spline1Dblend const * ptr{ Utils::mex_convert_mx_to_ptr<Spline1Dblend>( arg_in_1 ) };

    mwSize nx;
    real_type const * x{ Utils::mex_vector_pointer( arg_in_2, nx, CMD ": error in reading `x`" ) };
    real_type const   s{ Utils::mex_get_scalar_value( arg_in_3, CMD ": error in reading `s`" ) };
    real_type       * y{ Utils::mex_create_matrix_value( arg_out_0, nx, 1 ) };

    for ( mwSize i{0}; i < nx; ++i ) *y++ = ptr->eval_DD( *x++, s );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_eval_DDD(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define MEX_ERROR_MESSAGE_9 "y'''(x) = Spline1DblendMexWrapper('eval_DDD',OBJ,x,s): "
    #define CMD MEX_ERROR_MESSAGE_9

    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 3, CMD ": expected 3 input, nrhs = {}\n", nrhs );

    Spline1Dblend const * ptr{ Utils::mex_convert_mx_to_ptr<Spline1Dblend>( arg_in_1 ) };

    mwSize nx;
    real_type const * x{ Utils::mex_vector_pointer( arg_in_2, nx, CMD ": error in reading `x`" ) };
    real_type const   s{ Utils::mex_get_scalar_value( arg_in_3, CMD ": error in reading `s`" ) };
    real_type       * y{ Utils::mex_create_matrix_value( arg_out_0, nx, 1 ) };

    for ( mwSize i{0}; i < nx; ++i ) *y++ = ptr->eval_DDD( *x++, s );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_eval_DDDD(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define MEX_ERROR_MESSAGE_10 "y''''(x) = Spline1DblendMexWrapper('eval_DDDD',OBJ,x,s)"
    #define CMD MEX_ERROR_MESSAGE_10

    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 3, CMD ": expected 3 input, nrhs = {}\n", nrhs );

    Spline1Dblend const * ptr{ Utils::mex_convert_mx_to_ptr<Spline1Dblend>( arg_in_1 ) };

    mwSize nx;
    real_type const * x{ Utils::mex_vector_pointer(  arg_in_2, nx, CMD ": error in reading `x`" ) };
    real_type const   s{ Utils::mex_get_scalar_value( arg_in_3, CMD ": error in reading `s`" ) };
    real_type       * y{ Utils::mex_create_matrix_value( arg_out_0, nx, 1 ) };

    for ( mwSize i{0}; i < nx; ++i ) *y++ = ptr->eval_DDDD( *x++, s );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_eval_DDDDD(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define MEX_ERROR_MESSAGE_11 "y'''''(x) = Spline1DblendMexWrapper('eval_DDDDD',OBJ,x,s): "
    #define CMD MEX_ERROR_MESSAGE_11

    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 3, CMD ": expected 3 input, nrhs = {}\n", nrhs );

    Spline1Dblend const * ptr{ Utils::mex_convert_mx_to_ptr<Spline1Dblend>( arg_in_1 ) };

    mwSize nx;
    real_type const * x{ Utils::mex_vector_pointer( arg_in_2, nx, CMD ": error in reading `x`" ) };
    real_type const   s{ Utils::mex_get_scalar_value( arg_in_3, CMD ": error in reading `s`" ) };
    real_type       * y{ Utils::mex_create_matrix_value( arg_out_0, nx, 1 ) };

    for ( mwSize i{0}; i < nx; ++i ) *y++ = ptr->eval_DDDDD( *x++, s );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_x_begin(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define MEX_ERROR_MESSAGE_12 "x = Spline1DblendMexWrapper('x_begin',OBJ,s)"
    #define CMD MEX_ERROR_MESSAGE_12

    UTILS_MEX_ASSERT( nlhs == 0, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 3, CMD ": expected 3 input, nrhs = {}\n", nrhs );

    Spline1Dblend const * ptr{ Utils::mex_convert_mx_to_ptr<Spline1Dblend>( arg_in_1 ) };
    real_type const s{ Utils::mex_get_scalar_value( arg_in_2, CMD ": error in reading `s`" ) };
    Utils::mex_set_scalar_value( arg_out_0, ptr->x_begin(s) );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_y_begin(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define MEX_ERROR_MESSAGE_13 "y = Spline1DblendMexWrapper('y_begin',OBJ,s)"
    #define CMD MEX_ERROR_MESSAGE_13

    UTILS_MEX_ASSERT( nlhs == 0, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 3, CMD ": expected 3 input, nrhs = {}\n", nrhs );

    Spline1Dblend const * ptr{ Utils::mex_convert_mx_to_ptr<Spline1Dblend>( arg_in_1 ) };
    real_type const s{ Utils::mex_get_scalar_value( arg_in_2, CMD ": error in reading `s`" ) };
    Utils::mex_set_scalar_value( arg_out_0, ptr->y_begin(s) );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_x_end(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define MEX_ERROR_MESSAGE_14 "x = Spline1DblendMexWrapper('x_end',OBJ,s)"
    #define CMD MEX_ERROR_MESSAGE_14

    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 3, CMD ": expected 3 input, nrhs = {}\n", nrhs );

    Spline1Dblend const * ptr{ Utils::mex_convert_mx_to_ptr<Spline1Dblend>( arg_in_1 ) };
    real_type const s{ Utils::mex_get_scalar_value( arg_in_2, CMD ": error in reading `s`" ) };
    Utils::mex_set_scalar_value( arg_out_0, ptr->x_end(s) );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_y_end(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define MEX_ERROR_MESSAGE_15 "y = Spline1DblendMexWrapper('y_end',OBJ,s)"
    #define CMD MEX_ERROR_MESSAGE_15

    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 3, CMD ": expected 3 input, nrhs = {}\n", nrhs );

    Spline1Dblend const * ptr{ Utils::mex_convert_mx_to_ptr<Spline1Dblend>( arg_in_1 ) };
    real_type const s{ Utils::mex_get_scalar_value( arg_in_2, CMD ": error in reading `s`" ) };
    Utils::mex_set_scalar_value( arg_out_0, ptr->y_end(s) );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  typedef void (*DO_CMD)( int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[] );

  static unordered_map<string,DO_CMD> cmd_to_fun = {
    {"new",do_new},
    {"delete",do_delete},
    {"build",do_build},
    {"degree",do_degree},
    {"order",do_order},
    {"eval",do_eval},
    {"eval_D",do_eval_D},
    {"eval_DD",do_eval_DD},
    {"eval_DDD",do_eval_DDD},
    {"eval_DDDD",do_eval_DDDD},
    {"eval_DDDDD",do_eval_DDDDD},
    {"x_begin",do_x_begin},
    {"y_begin",do_y_begin},
    {"x_end",do_x_end},
    {"y_end",do_y_end}
  };

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  #define MEX_ERROR_MESSAGE \
"%======================================================================%\n" \
"Spline1DblendMexWrapper:  Compute spline curve\n" \
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
"\n" \
"%======================================================================%\n" \
"%                                                                      %\n" \
"%  Autor: Enrico Bertolazzi                                            %\n" \
"%         Department of Industrial Engineering                         %\n" \
"%         University of Trento                                         %\n" \
"%         enrico.bertolazzi@unitn.it                                   %\n" \
"%                                                                      %\n" \
"%======================================================================%\n"

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

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
      mexErrMsgTxt( fmt::format( "Spline1DblendMexWrapper Error: {}", e.what() ).c_str() );
    } catch (...) {
      mexErrMsgTxt("Spline1DblendMexWrapper failed\n");
    }

  }

}
