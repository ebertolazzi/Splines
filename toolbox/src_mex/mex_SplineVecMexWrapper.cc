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

#define ASSERT(COND,MSG)                           \
  if ( !(COND) ) {                                 \
    std::ostringstream ost;                        \
    ost << "SplineVecMexWrapper: " << MSG << '\n'; \
    mexErrMsgTxt(ost.str().c_str());               \
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
          int nrhs, mxArray const *[] ) {

    #define MEX_ERROR_MESSAGE_1 "SplineVecMexWrapper( 'new' )"
    #define CMD MEX_ERROR_MESSAGE_1

    UTILS_MEX_ASSERT( nrhs == 1, CMD ": expected 1 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );

    arg_out_0 = Utils::mex_convert_ptr_to_mx<SplineVec>( new SplineVec() );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_delete( int nlhs, mxArray       *[],
             int nrhs, mxArray const *prhs[] ) {

    #define MEX_ERROR_MESSAGE_2 "SplineVecMexWrapper( 'delete', OBJ ): "
    #define CMD MEX_ERROR_MESSAGE_2

    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 0, CMD ": expected 0 output, nlhs = {}\n", nlhs );

    // Destroy the C++ object
    Utils::mex_destroy_object<SplineVec>(arg_in_1);

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_setup( int nlhs, mxArray       *[],
            int nrhs, mxArray const *prhs[] ) {

    #define MEX_ERROR_MESSAGE_3 "SplineVecMexWrapper( 'setup', obj, Y )"
    #define CMD MEX_ERROR_MESSAGE_3

    UTILS_MEX_ASSERT( nrhs == 3, CMD ": expected 3 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 0, CMD ": expected 0 output, nlhs = {}\n", nlhs );

    SplineVec * ptr = Utils::mex_convert_mx_to_ptr<SplineVec>( arg_in_1 );

    mwSize dim, npts;
    real_type const * Y = Utils::mex_matrix_pointer( arg_in_2, dim, npts, CMD ": error in reading 'Y'" );
    ptr->setup( dim, npts, Y, dim );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_knots( int nlhs, mxArray       *[],
            int nrhs, mxArray const *prhs[] ) {

    #define MEX_ERROR_MESSAGE_4 "SplineVecMexWrapper( 'knots', obj, X )"
    #define CMD MEX_ERROR_MESSAGE_4

    UTILS_MEX_ASSERT( nrhs == 3, CMD ": expected 3 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 0, CMD ": expected 0 output, nlhs = {}\n", nlhs );

    SplineVec * ptr = Utils::mex_convert_mx_to_ptr<SplineVec>( arg_in_1 );

    mwSize npts;
    real_type const * X = Utils::mex_vector_pointer( arg_in_2, npts, CMD ": error in reading 'X'" );
    UTILS_MEX_ASSERT(
      npts == mwSize(ptr->num_points()),
      CMD "size(X) = {} must be = {}",
      npts, ptr->dimension()
    );
    ptr->set_knots( X );
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_chordal( int nlhs, mxArray       *[],
              int nrhs, mxArray const *prhs[] ) {

    #define MEX_ERROR_MESSAGE_5 "SplineVecMexWrapper( 'chordal', obj )"
    #define CMD MEX_ERROR_MESSAGE_5

    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 0, CMD ": expected 0 output, nlhs = {}\n", nlhs );

    SplineVec * ptr = Utils::mex_convert_mx_to_ptr<SplineVec>( arg_in_1 );

    ptr->set_knots_chord_length();
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_centripetal( int nlhs, mxArray       *[],
                  int nrhs, mxArray const *prhs[] ) {

    #define MEX_ERROR_MESSAGE_6 "SplineVecMexWrapper( 'centripetal', obj )"
    #define CMD MEX_ERROR_MESSAGE_6

    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 0, CMD ": expected 0 output, nlhs = {}\n", nlhs );

    SplineVec * ptr = Utils::mex_convert_mx_to_ptr<SplineVec>( arg_in_1 );

    ptr->set_knots_centripetal();
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_CatmullRom( int nlhs, mxArray       *[],
                 int nrhs, mxArray const *prhs[] ) {

    #define MEX_ERROR_MESSAGE_7 "SplineVecMexWrapper( 'CatmullRom', obj )"
    #define CMD MEX_ERROR_MESSAGE_7

    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 0, CMD ": expected 0 output, nlhs = {}\n", nlhs );

    SplineVec * ptr = Utils::mex_convert_mx_to_ptr<SplineVec>( arg_in_1 );

    ptr->catmull_rom();
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_get_knots( int nlhs, mxArray       *plhs[],
                int nrhs, mxArray const *prhs[] ) {

    #define MEX_ERROR_MESSAGE_8 "knots = SplineVecMexWrapper( 'get_knots', obj )"
    #define CMD MEX_ERROR_MESSAGE_8

    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );

    SplineVec * ptr = Utils::mex_convert_mx_to_ptr<SplineVec>( arg_in_1 );

    //! return the number of support points of the splines
    integer N = ptr->num_points();
    real_type * X = Utils::mex_create_matrix_value( arg_out_0, 1, N );
    for ( integer i{0}; i < N; ++i ) X[i] = ptr->x_node( i );
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_eval( int nlhs, mxArray       *plhs[],
           int nrhs, mxArray const *prhs[] ) {

    #define MEX_ERROR_MESSAGE_9 "y(x) = SplineVecMexWrapper( 'eval', obj, x )"
    #define CMD MEX_ERROR_MESSAGE_9

    UTILS_MEX_ASSERT( nrhs == 3, CMD ": expected 3 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );

    SplineVec * ptr = Utils::mex_convert_mx_to_ptr<SplineVec>( arg_in_1 );

    mwSize nx;
    real_type const * x = Utils::mex_vector_pointer( arg_in_2, nx, CMD ": error in reading `x`" );

    mwSize dim = ptr->dimension();
    real_type * Y = Utils::mex_create_matrix_value( arg_out_0, dim, nx );

    for ( mwSize i{0}; i < nx; ++i, Y += dim )
      ptr->eval( x[i], Y, 1 );
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_eval_D( int nlhs, mxArray       *plhs[],
             int nrhs, mxArray const *prhs[] ) {

    #define MEX_ERROR_MESSAGE_10 "y'(x) = SplineVecMexWrapper( 'eval_D', obj, x )"
    #define CMD MEX_ERROR_MESSAGE_10

    UTILS_MEX_ASSERT( nrhs == 3, CMD ": expected 3 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );

    SplineVec * ptr = Utils::mex_convert_mx_to_ptr<SplineVec>( arg_in_1 );

    mwSize nx;
    real_type const * x = Utils::mex_vector_pointer( arg_in_2, nx, CMD ": error in reading `x`" );

    mwSize dim = ptr->dimension();
    real_type * Y = Utils::mex_create_matrix_value( arg_out_0, dim, nx );

    for ( mwSize i{0}; i < nx; ++i, Y += dim )
      ptr->eval_D( x[i], Y, 1 );
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_eval_DD( int nlhs, mxArray       *plhs[],
              int nrhs, mxArray const *prhs[] ) {

    #define MEX_ERROR_MESSAGE_11 "y''(x) = SplineVecMexWrapper( 'eval_DD', obj, x )"
    #define CMD MEX_ERROR_MESSAGE_11

    UTILS_MEX_ASSERT( nrhs == 3, CMD ": expected 3 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );

    SplineVec * ptr = Utils::mex_convert_mx_to_ptr<SplineVec>( arg_in_1 );

    mwSize nx;
    real_type const * x = Utils::mex_vector_pointer( arg_in_2, nx, CMD ": error in reading `x`" );

    mwSize dim = ptr->dimension();
    real_type * Y = Utils::mex_create_matrix_value( arg_out_0, dim, nx );

    for ( mwSize i{0}; i < nx; ++i, Y += dim )
      ptr->eval_DD( x[i], Y, 1 );
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_eval_DDD( int nlhs, mxArray       *plhs[],
               int nrhs, mxArray const *prhs[] ) {

    #define MEX_ERROR_MESSAGE_12 "y'''(x) = SplineVecMexWrapper( 'eval_DDD', obj, x )"
    #define CMD MEX_ERROR_MESSAGE_12

    UTILS_MEX_ASSERT( nrhs == 3, CMD ": expected 3 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );

    SplineVec * ptr = Utils::mex_convert_mx_to_ptr<SplineVec>( arg_in_1 );

    mwSize nx;
    real_type const * x = Utils::mex_vector_pointer( arg_in_2, nx, CMD ": error in reading `x`" );

    mwSize dim = ptr->dimension();
    real_type * Y = Utils::mex_create_matrix_value( arg_out_0, dim, nx );

    for ( mwSize i{0}; i < nx; ++i, Y += dim )
      ptr->eval_DDD( x[i], Y, 1 );
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_eval_curvature( int nlhs, mxArray       *plhs[],
                     int nrhs, mxArray const *prhs[] ) {

    #define MEX_ERROR_MESSAGE_13 "kappa(x) = SplineVecMexWrapper( 'eval_curvature', obj, x )"
    #define CMD MEX_ERROR_MESSAGE_13

    UTILS_MEX_ASSERT( nrhs == 3, CMD ": expected 3 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );

    SplineVec * ptr = Utils::mex_convert_mx_to_ptr<SplineVec>( arg_in_1 );

    mwSize nx;
    real_type const * x = Utils::mex_vector_pointer( arg_in_2, nx, CMD ": error in reading `x`" );

    real_type * curvature = Utils::mex_create_matrix_value( arg_out_0, 1, nx );

    for ( mwSize i{0}; i < nx; ++i )
      curvature[i] = ptr->curvature( x[i] );
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_eval_curvature_D( int nlhs, mxArray       *plhs[],
                       int nrhs, mxArray const *prhs[] ) {

    #define MEX_ERROR_MESSAGE_14 "kappa'(x) = SplineVecMexWrapper( 'eval_curvature_D', obj, x )"
    #define CMD MEX_ERROR_MESSAGE_14

    UTILS_MEX_ASSERT( nrhs == 3, CMD ": expected 3 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );

    SplineVec * ptr = Utils::mex_convert_mx_to_ptr<SplineVec>( arg_in_1 );

    mwSize nx;
    real_type const * x = Utils::mex_vector_pointer( arg_in_2, nx, CMD ": error in reading `x`" );

    real_type * curvature_D = Utils::mex_create_matrix_value( arg_out_0, 1, nx );

    for ( mwSize i{0}; i < nx; ++i )
      curvature_D[i] = ptr->curvature_D( x[i] );
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_tmin( int nlhs, mxArray       *plhs[],
           int nrhs, mxArray const *prhs[] ) {

    #define MEX_ERROR_MESSAGE_15 "tmin = SplineVecMexWrapper('tmin',OBJ)"
    #define CMD MEX_ERROR_MESSAGE_15

    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 input, nrhs = {}\n", nrhs );

    SplineVec * ptr = Utils::mex_convert_mx_to_ptr<SplineVec>( arg_in_1 );

    Utils::mex_set_scalar_value( arg_out_0, ptr->x_min() );
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_tmax( int nlhs, mxArray       *plhs[],
           int nrhs, mxArray const *prhs[] ) {

    #define MEX_ERROR_MESSAGE_16 "tmax = SplineVecMexWrapper('tmax',OBJ)"
    #define CMD MEX_ERROR_MESSAGE_16

    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 input, nrhs = {}\n", nrhs );

    SplineVec * ptr = Utils::mex_convert_mx_to_ptr<SplineVec>( arg_in_1 );

    Utils::mex_set_scalar_value( arg_out_0, ptr->x_max() );
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  typedef void (*DO_CMD)( int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[] );

  static unordered_map<string,DO_CMD> cmd_to_fun = {
    {"new",do_new},
    {"delete",do_delete},
    {"setup",do_setup},
    {"knots",do_knots},
    {"chordal",do_chordal},
    {"centripetal",do_centripetal},
    {"CatmullRom",do_CatmullRom},
    {"get_knots",do_get_knots},
    {"eval",do_eval},
    {"eval_D",do_eval_D},
    {"eval_DD",do_eval_DD},
    {"eval_DDD",do_eval_DDD},
    {"eval_curvature",do_eval_curvature},
    {"eval_curvature_D",do_eval_curvature_D},
    {"tmin",do_tmin},
    {"tmax",do_tmax}
  };

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  #define MEX_ERROR_MESSAGE \
"%======================================================================%\n" \
"SplineVecMexWrapper:  Compute Hermite base\n" \
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
      mexErrMsgTxt( fmt::format( "SplineVecMexWrapper Error: {}", e.what() ).c_str() );
    } catch (...) {
      mexErrMsgTxt("SplineVecMexWrapper failed\n");
    }

  }

}
