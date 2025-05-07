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
    ost << "Spline1DMexWrapper: " << MSG << '\n'; \
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
  do_new(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define MEX_ERROR_MESSAGE_1 "Spline1DMexWrapper( 'new', kind )"
    #define CMD MEX_ERROR_MESSAGE_1

    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );

    UTILS_MEX_ASSERT(
      mxIsChar(arg_in_1),
      CMD ": second argument must be a string, found ``{}''\n",
      mxGetClassName(arg_in_1)
    );

    // the first argument must be a string
    string tname = mxArrayToString(arg_in_1);

    if      ( tname == "linear"  ) arg_out_0 = Utils::mex_convert_ptr_to_mx<Spline>( new Splines::LinearSpline()  );
    else if ( tname == "cubic"   ) arg_out_0 = Utils::mex_convert_ptr_to_mx<Spline>( new Splines::CubicSpline()   );
    else if ( tname == "akima"   ) arg_out_0 = Utils::mex_convert_ptr_to_mx<Spline>( new Splines::AkimaSpline()   );
    else if ( tname == "bessel"  ) arg_out_0 = Utils::mex_convert_ptr_to_mx<Spline>( new Splines::BesselSpline()  );
    else if ( tname == "pchip"   ) arg_out_0 = Utils::mex_convert_ptr_to_mx<Spline>( new Splines::PchipSpline()   );
    else if ( tname == "hermite" ) arg_out_0 = Utils::mex_convert_ptr_to_mx<Spline>( new Splines::HermiteSpline() );
    else if ( tname == "quintic" ) arg_out_0 = Utils::mex_convert_ptr_to_mx<Spline>( new Splines::QuinticSpline() );
    else {
      UTILS_MEX_ASSERT0(
        false,
        CMD ": second argument must be one of the strings:\n"
        "'linear', 'cubic', 'akima', 'bessel', 'hermite', 'pchip', 'quintic'"
      );
    }

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_delete(
    int nlhs, mxArray       *[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_2 "Spline1DMexWrapper( 'delete', OBJ )"
    #define CMD MEX_ERROR_MESSAGE_2

    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 0, CMD ": expected 0 output, nlhs = {}\n", nlhs );

    // Destroy the C++ object
    Utils::mex_destroy_object<Spline>(arg_in_1);
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_build(
    int nlhs, mxArray       *[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define MEX_ERROR_MESSAGE_3 "Spline1DMexWrapper('build',OBJ,x,y[,yp or subtype])"
    #define CMD MEX_ERROR_MESSAGE_3

    UTILS_MEX_ASSERT( nlhs == 0, CMD ": expected no output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT(
      nrhs == 4 || nrhs == 5,
      CMD "expected 4 or 5 input, nrhs = {}\n", nrhs
    );

    mwSize nx, ny;
    real_type const * x  = nullptr;
    real_type const * y  = nullptr;
    real_type const * yp = nullptr;
    x = Utils::mex_vector_pointer( arg_in_2, nx, CMD ": error in reading 'x'" );
    y = Utils::mex_vector_pointer( arg_in_3, ny, CMD ": error in reading 'y'" );

    UTILS_MEX_ASSERT0( nx == ny, CMD ": lenght of 'x' must be the lenght of 'y'" );

    Spline * ptr = Utils::mex_convert_mx_to_ptr<Spline>( arg_in_1 );

    if ( nrhs == 5 ) {
      if ( mxIsChar(arg_in_4) ) {
        UTILS_MEX_ASSERT0(
          ptr->type() == SplineType1D::CUBIC || ptr->type() == SplineType1D::QUINTIC,
          "subtype can be specifiewd only for Cubic or Quintic spline"
        );
        string subtype = mxArrayToString(arg_in_4);
        switch ( ptr->type() ) {
        case SplineType1D::CONSTANT:
        case SplineType1D::LINEAR:
        case SplineType1D::AKIMA:
        case SplineType1D::BESSEL:
        case SplineType1D::PCHIP:
        case SplineType1D::HERMITE:
          break;
        case SplineType1D::CUBIC:
          if ( subtype == "extrapolate" ) {
            static_cast<CubicSpline*>(ptr)->set_initial_BC( CubicSpline_BC::EXTRAPOLATE );
            static_cast<CubicSpline*>(ptr)->set_final_BC( CubicSpline_BC::EXTRAPOLATE );
          } else if ( subtype == "natural" ) {
            static_cast<CubicSpline*>(ptr)->set_initial_BC( CubicSpline_BC::NATURAL );
            static_cast<CubicSpline*>(ptr)->set_final_BC( CubicSpline_BC::NATURAL );
          } else if ( subtype == "parabolic"  ) {
            static_cast<CubicSpline*>(ptr)->set_initial_BC( CubicSpline_BC::PARABOLIC_RUNOUT );
            static_cast<CubicSpline*>(ptr)->set_final_BC( CubicSpline_BC::PARABOLIC_RUNOUT );
          } else if ( subtype == "not_a_knot" ) {
            static_cast<CubicSpline*>(ptr)->set_initial_BC( CubicSpline_BC::NOT_A_KNOT );
            static_cast<CubicSpline*>(ptr)->set_final_BC( CubicSpline_BC::NOT_A_KNOT );
          } else {
            UTILS_MEX_ASSERT(
              false,
              CMD "subtype: {} must be in:\n"
              "['extrapolate','natural','parabolic','not_a_knot']\n",
              subtype
            );
          }
          break;
        case SplineType1D::QUINTIC:
          if ( subtype == "cubic" ) {
            static_cast<QuinticSpline*>(ptr)->setQuinticType( QuinticSpline_sub_type::CUBIC );
          } else if ( subtype == "pchip"  ) {
            static_cast<QuinticSpline*>(ptr)->setQuinticType( QuinticSpline_sub_type::PCHIP );
          } else if ( subtype == "akima" ) {
            static_cast<QuinticSpline*>(ptr)->setQuinticType( QuinticSpline_sub_type::AKIMA );
          } else if ( subtype == "bessel" ) {
            static_cast<QuinticSpline*>(ptr)->setQuinticType( QuinticSpline_sub_type::BESSEL );
          } else {
            UTILS_MEX_ASSERT(
              false,
              CMD "subtype: {} must be in:\n"
              "['cubic','pchip','akima','bessel']\n",
              subtype
            );
          }
          break;
        default:
          break;
        }
      } else {
        UTILS_MEX_ASSERT0(
          ptr->type() == SplineType1D::HERMITE,
          CMD ": yp can be specified only for Hermite type Spline"
        );
        mwSize nyp;
        yp = Utils::mex_vector_pointer( arg_in_4, nyp, CMD "error in reading 'yp'" );
        UTILS_MEX_ASSERT0(
          ny == nyp,
          CMD ": lenght of 'yp' must be the lenght of 'y'"
        );
      }
    }
    switch ( ptr->type() ) {
    case SplineType1D::HERMITE:
      UTILS_MEX_ASSERT0(
        x != nullptr && y != nullptr && yp != nullptr,
        CMD ": something go wrong in reading x, y, or yp"
      );
      static_cast<HermiteSpline*>(ptr)->build( x, y, yp, nx );
      break;
    default:
      UTILS_MEX_ASSERT0(
        x != nullptr && y != nullptr,
        CMD ": something go wrong in reading x or y"
      );
      ptr->build( x, y, nx );
      break;
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

    #define MEX_ERROR_MESSAGE_4 "degree = Spline1DMexWrapper('degree',OBJ)"
    #define CMD MEX_ERROR_MESSAGE_4

    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 input, nrhs = {}\n", nrhs );
    Spline * ptr = Utils::mex_convert_mx_to_ptr<Spline>( arg_in_1 );
    Utils::mex_set_scalar_int32( arg_out_0, ptr->order()-1 );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_order(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define MEX_ERROR_MESSAGE_5 "order = Spline1DMexWrapper('order',OBJ)"
    #define CMD MEX_ERROR_MESSAGE_5

    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 input, nrhs = {}\n", nrhs );
    Spline * ptr = Utils::mex_convert_mx_to_ptr<Spline>( arg_in_1 );
    Utils::mex_set_scalar_int32( arg_out_0, ptr->order() );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_coeffs(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define MEX_ERROR_MESSAGE_6 "[coeffs,nodes] = Spline1DMexWrapper('coeffs',OBJ)"
    #define CMD MEX_ERROR_MESSAGE_6

    UTILS_MEX_ASSERT( nlhs == 2, CMD ": expected 2 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 input, nrhs = {}\n", nrhs );

    Spline * ptr = Utils::mex_convert_mx_to_ptr<Spline>( arg_in_1 );

    integer ord = ptr->order();
    integer npt = ptr->num_points();

    real_type * cfs   = Utils::mex_create_matrix_value( arg_out_0, npt-1, ord );
    real_type * nodes = Utils::mex_create_matrix_value( arg_out_1, npt,   1   );

    ptr->coeffs( cfs, nodes, false );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_y_min_max(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_7 "struct = Spline1DMexWrapper('y_min_max',OBJ)"
    #define CMD MEX_ERROR_MESSAGE_7

    UTILS_MEX_ASSERT( nlhs == 1 || nlhs == 6, CMD ": expected 1 or 6 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 input, nrhs = {}\n", nrhs );
    Spline * ptr = Utils::mex_convert_mx_to_ptr<Spline>( arg_in_1 );

    if ( nlhs == 1 ) {

      vector<integer>   i_min_pos, i_max_pos;
      vector<real_type> x_min_pos, y_min, x_max_pos, y_max;

      ptr->y_min_max(
        i_min_pos, x_min_pos, y_min,
        i_max_pos, x_max_pos, y_max
      );

      static char const * fieldnames[] = {
        "i_min_pos",
        "x_min_pos",
        "y_min",
        "i_max_pos",
        "x_max_pos",
        "y_max"
      };

      arg_out_0 = mxCreateStructMatrix(1,1,6,fieldnames);

      mxArray * mx_i_min_pos;
      mxArray * mx_x_min_pos;
      mxArray * mx_y_min;
      mxArray * mx_i_max_pos;
      mxArray * mx_x_max_pos;
      mxArray * mx_y_max;

      real_type * ptr_i_min_pos = Utils::mex_create_matrix_value( mx_i_min_pos, i_min_pos.size(), 1 );
      real_type * ptr_x_min_pos = Utils::mex_create_matrix_value( mx_x_min_pos, x_min_pos.size(), 1 );
      real_type * ptr_y_min     = Utils::mex_create_matrix_value( mx_y_min,     y_min.size(),     1 );
      real_type * ptr_i_max_pos = Utils::mex_create_matrix_value( mx_i_max_pos, i_max_pos.size(), 1 );
      real_type * ptr_x_max_pos = Utils::mex_create_matrix_value( mx_x_max_pos, x_max_pos.size(), 1 );
      real_type * ptr_y_max     = Utils::mex_create_matrix_value( mx_y_max,     y_max.size(),     1 );

      std::copy( i_min_pos.begin(), i_min_pos.end(), ptr_i_min_pos );
      std::copy( x_min_pos.begin(), x_min_pos.end(), ptr_x_min_pos );
      std::copy( y_min.begin(),     y_min.end(),     ptr_y_min     );
      std::copy( i_max_pos.begin(), i_max_pos.end(), ptr_i_max_pos );
      std::copy( x_max_pos.begin(), x_max_pos.end(), ptr_x_max_pos );
      std::copy( y_max.begin(),     y_max.end(),     ptr_y_max     );

      mxSetFieldByNumber( arg_out_0, 0, 0, mx_i_min_pos );
      mxSetFieldByNumber( arg_out_0, 0, 1, mx_x_min_pos );
      mxSetFieldByNumber( arg_out_0, 0, 2, mx_y_min     );
      mxSetFieldByNumber( arg_out_0, 0, 3, mx_i_max_pos );
      mxSetFieldByNumber( arg_out_0, 0, 4, mx_x_max_pos );
      mxSetFieldByNumber( arg_out_0, 0, 5, mx_y_max     );

    } else {

      integer   i_min_pos, i_max_pos;
      real_type x_min_pos, y_min, x_max_pos, y_max;

      ptr->y_min_max(
        i_min_pos, x_min_pos, y_min,
        i_max_pos, x_max_pos, y_max
      );

      Utils::mex_set_scalar_int32  ( arg_out_0, i_min_pos );
      Utils::mex_set_scalar_value( arg_out_1, x_min_pos );
      Utils::mex_set_scalar_value( arg_out_2, y_min     );
      Utils::mex_set_scalar_int32  ( arg_out_3, i_max_pos );
      Utils::mex_set_scalar_value( arg_out_4, x_max_pos );
      Utils::mex_set_scalar_value( arg_out_5, y_max     );

    }
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_eval(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define MEX_ERROR_MESSAGE_8 "y = Spline1DMexWrapper('eval',OBJ,x)"
    #define CMD MEX_ERROR_MESSAGE_8

    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 3, CMD ": expected 3 input, nrhs = {}\n", nrhs );

    Spline * ptr = Utils::mex_convert_mx_to_ptr<Spline>( arg_in_1 );

    mwSize nx;
    real_type const * x = Utils::mex_vector_pointer(
      arg_in_2, nx, CMD ": error in reading `x`"
    );
    real_type * y = Utils::mex_create_matrix_value( arg_out_0, nx, 1 );

    for ( mwSize i{0}; i < nx; ++i ) *y++ = ptr->eval( *x++ );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_eval_D(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define MEX_ERROR_MESSAGE_9 "y'(x) = Spline1DMexWrapper('eval_D',OBJ,x)"
    #define CMD MEX_ERROR_MESSAGE_9

    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 3, CMD ": expected 3 input, nrhs = {}\n", nrhs );

    Spline * ptr = Utils::mex_convert_mx_to_ptr<Spline>( arg_in_1 );

    mwSize nx;
    real_type const * x = Utils::mex_vector_pointer(
      arg_in_2, nx, CMD ": error in reading `x`"
    );
    real_type * y = Utils::mex_create_matrix_value( arg_out_0, nx, 1 );

    for ( mwSize i{0}; i < nx; ++i ) *y++ = ptr->eval_D( *x++ );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_eval_DD(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define MEX_ERROR_MESSAGE_10 "y''(x) = Spline1DMexWrapper('eval_DD',OBJ,x): "
    #define CMD MEX_ERROR_MESSAGE_10

    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 3, CMD ": expected 3 input, nrhs = {}\n", nrhs );

    Spline * ptr = Utils::mex_convert_mx_to_ptr<Spline>( arg_in_1 );

    mwSize nx;
    real_type const * x = Utils::mex_vector_pointer(
      arg_in_2, nx, CMD ": error in reading `x`"
    );
    real_type * y = Utils::mex_create_matrix_value( arg_out_0, nx, 1 );

    for ( mwSize i{0}; i < nx; ++i ) *y++ = ptr->eval_DD( *x++ );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_eval_DDD(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define MEX_ERROR_MESSAGE_11 "y'''(x) = Spline1DMexWrapper('eval_DDD',OBJ,x): "
    #define CMD MEX_ERROR_MESSAGE_11

    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 3, CMD ": expected 3 input, nrhs = {}\n", nrhs );

    Spline * ptr = Utils::mex_convert_mx_to_ptr<Spline>( arg_in_1 );

    mwSize nx;
    real_type const * x = Utils::mex_vector_pointer(
      arg_in_2, nx, CMD ": error in reading `x`"
    );
    real_type * y = Utils::mex_create_matrix_value( arg_out_0, nx, 1 );

    for ( mwSize i{0}; i < nx; ++i ) *y++ = ptr->eval_DDD( *x++ );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_eval_DDDD(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define MEX_ERROR_MESSAGE_12 "y''''(x) = Spline1DMexWrapper('eval_DDDD',OBJ,x)"
    #define CMD MEX_ERROR_MESSAGE_12

    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 3, CMD ": expected 3 input, nrhs = {}\n", nrhs );

    Spline * ptr = Utils::mex_convert_mx_to_ptr<Spline>( arg_in_1 );

    mwSize nx;
    real_type const * x = Utils::mex_vector_pointer(
      arg_in_2, nx, CMD ": error in reading `x`"
    );
    real_type * y = Utils::mex_create_matrix_value( arg_out_0, nx, 1 );

    for ( mwSize i{0}; i < nx; ++i ) *y++ = ptr->eval_DDDD( *x++ );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_eval_DDDDD(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define MEX_ERROR_MESSAGE_13 "y'''''(x) = Spline1DMexWrapper('eval_DDDDD',OBJ,x): "
    #define CMD MEX_ERROR_MESSAGE_13

    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 3, CMD ": expected 3 input, nrhs = {}\n", nrhs );

    Spline * ptr = Utils::mex_convert_mx_to_ptr<Spline>( arg_in_1 );

    mwSize nx;
    real_type const * x = Utils::mex_vector_pointer(
      arg_in_2, nx, CMD ": error in reading `x`"
    );
    real_type * y = Utils::mex_create_matrix_value( arg_out_0, nx, 1 );

    for ( mwSize i{0}; i < nx; ++i ) *y++ = ptr->eval_DDDDD( *x++ );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_make_closed(
    int nlhs, mxArray       *[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define MEX_ERROR_MESSAGE_14 "Spline1DMexWrapper('make_closed',OBJ)"
    #define CMD MEX_ERROR_MESSAGE_14

    UTILS_MEX_ASSERT( nlhs == 0, CMD ": expected 0 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 input, nrhs = {}\n", nrhs );

    Spline * ptr = Utils::mex_convert_mx_to_ptr<Spline>( arg_in_1 );
    ptr->make_closed();

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_make_opened(
    int nlhs, mxArray       *[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define MEX_ERROR_MESSAGE_15 "Spline1DMexWrapper('make_opened',OBJ)"
    #define CMD MEX_ERROR_MESSAGE_15

    UTILS_MEX_ASSERT( nlhs == 0, CMD ": expected 0 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 input, nrhs = {}\n", nrhs );

    Spline * ptr = Utils::mex_convert_mx_to_ptr<Spline>( arg_in_1 );
    ptr->make_opened();

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_is_closed(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define MEX_ERROR_MESSAGE_16 "Spline1DMexWrapper('is_closed',OBJ)"
    #define CMD MEX_ERROR_MESSAGE_16

    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 input, nrhs = {}\n", nrhs );

    Spline * ptr = Utils::mex_convert_mx_to_ptr<Spline>( arg_in_1 );
    Utils::mex_set_scalar_bool( arg_out_0, ptr->is_closed() );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_make_bounded(
    int nlhs, mxArray       *[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define MEX_ERROR_MESSAGE_17 "Spline1DMexWrapper('make_bounded',OBJ)"
    #define CMD MEX_ERROR_MESSAGE_17

    UTILS_MEX_ASSERT( nlhs == 0, CMD ": expected 0 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 input, nrhs = {}\n", nrhs );

    Spline * ptr = Utils::mex_convert_mx_to_ptr<Spline>( arg_in_1 );
    ptr->make_bounded();

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_make_unbounded(
    int nlhs, mxArray       *[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define MEX_ERROR_MESSAGE_18 "Spline1DMexWrapper('make_unbounded',OBJ)"
    #define CMD MEX_ERROR_MESSAGE_18

    UTILS_MEX_ASSERT( nlhs == 0, CMD ": expected 0 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 input, nrhs = {}\n", nrhs );

    Spline * ptr = Utils::mex_convert_mx_to_ptr<Spline>( arg_in_1 );
    ptr->make_unbounded();

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_is_bounded(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define MEX_ERROR_MESSAGE_19 "res = Spline1DMexWrapper('is_bounded',OBJ)"
    #define CMD MEX_ERROR_MESSAGE_19

    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 input, nrhs = {}\n", nrhs );

    Spline * ptr = Utils::mex_convert_mx_to_ptr<Spline>( arg_in_1 );
    Utils::mex_set_scalar_bool( arg_out_0, ptr->is_bounded() );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_make_extended_constant(
    int nlhs, mxArray       *[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define MEX_ERROR_MESSAGE_20 "Spline1DMexWrapper('make_extended_constant',OBJ)"
    #define CMD MEX_ERROR_MESSAGE_20

    UTILS_MEX_ASSERT( nlhs == 0, CMD ": expected 0 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 input, nrhs = {}\n", nrhs );

    Spline * ptr = Utils::mex_convert_mx_to_ptr<Spline>( arg_in_1 );
    ptr->make_extended_constant();

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_make_extended_not_constant(
    int nlhs, mxArray       *[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define MEX_ERROR_MESSAGE_21 "Spline1DMexWrapper('make_extended_not_constant',OBJ)"
    #define CMD MEX_ERROR_MESSAGE_21

    UTILS_MEX_ASSERT( nlhs == 0, CMD ": expected 0 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 input, nrhs = {}\n", nrhs );

    Spline * ptr = Utils::mex_convert_mx_to_ptr<Spline>( arg_in_1 );
    ptr->make_extended_not_constant();

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_is_extended_constant(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define MEX_ERROR_MESSAGE_22 "res = Spline1DMexWrapper('is_extended_constant',OBJ)"
    #define CMD MEX_ERROR_MESSAGE_22

    UTILS_MEX_ASSERT( nlhs == 0, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 input, nrhs = {}\n", nrhs );

    Spline * ptr = Utils::mex_convert_mx_to_ptr<Spline>( arg_in_1 );
    Utils::mex_set_scalar_bool( arg_out_0, ptr->is_extended_constant() );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_x_begin(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define MEX_ERROR_MESSAGE_23 "x = Spline1DMexWrapper('x_begin',OBJ)"
    #define CMD MEX_ERROR_MESSAGE_23

    UTILS_MEX_ASSERT( nlhs == 0, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 input, nrhs = {}\n", nrhs );

    Spline * ptr = Utils::mex_convert_mx_to_ptr<Spline>( arg_in_1 );
    Utils::mex_set_scalar_value( arg_out_0, ptr->x_begin() );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_y_begin(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define MEX_ERROR_MESSAGE_24 "y = Spline1DMexWrapper('y_begin',OBJ)"
    #define CMD MEX_ERROR_MESSAGE_24

    UTILS_MEX_ASSERT( nlhs == 0, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 input, nrhs = {}\n", nrhs );

    Spline * ptr = Utils::mex_convert_mx_to_ptr<Spline>( arg_in_1 );
    Utils::mex_set_scalar_value( arg_out_0, ptr->y_begin() );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_x_end(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define MEX_ERROR_MESSAGE_25 "x = Spline1DMexWrapper('x_end',OBJ)"
    #define CMD MEX_ERROR_MESSAGE_25

    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 input, nrhs = {}\n", nrhs );

    Spline * ptr = Utils::mex_convert_mx_to_ptr<Spline>( arg_in_1 );
    Utils::mex_set_scalar_value( arg_out_0, ptr->x_end() );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_y_end(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define MEX_ERROR_MESSAGE_26 "y = Spline1DMexWrapper('y_end',OBJ)"
    #define CMD MEX_ERROR_MESSAGE_26

    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 input, nrhs = {}\n", nrhs );

    Spline * ptr = Utils::mex_convert_mx_to_ptr<Spline>( arg_in_1 );
    Utils::mex_set_scalar_value( arg_out_0, ptr->y_end() );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_x_min(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define MEX_ERROR_MESSAGE_27 "x = Spline1DMexWrapper('x_min',OBJ)"
    #define CMD MEX_ERROR_MESSAGE_27

    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 input, nrhs = {}\n", nrhs );

    Spline * ptr = Utils::mex_convert_mx_to_ptr<Spline>( arg_in_1 );
    Utils::mex_set_scalar_value( arg_out_0, ptr->x_min() );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_y_min(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define MEX_ERROR_MESSAGE_28 "y = Spline1DMexWrapper('y_min',OBJ)"
    #define CMD MEX_ERROR_MESSAGE_28

    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 input, nrhs = {}\n", nrhs );

    Spline * ptr = Utils::mex_convert_mx_to_ptr<Spline>( arg_in_1 );
    Utils::mex_set_scalar_value( arg_out_0, ptr->y_min() );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_x_max(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define MEX_ERROR_MESSAGE_29 "x = Spline1DMexWrapper('x_max',OBJ)"
    #define CMD MEX_ERROR_MESSAGE_29

    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 input, nrhs = {}\n", nrhs );

    Spline * ptr = Utils::mex_convert_mx_to_ptr<Spline>( arg_in_1 );
    Utils::mex_set_scalar_value( arg_out_0, ptr->x_max() );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_y_max(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define MEX_ERROR_MESSAGE_30 "y = Spline1DMexWrapper('y_max',OBJ)"
    #define CMD MEX_ERROR_MESSAGE_30

    UTILS_MEX_ASSERT( nlhs == 1, CMD ": expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD ": expected 2 input, nrhs = {}\n", nrhs );

    Spline * ptr = Utils::mex_convert_mx_to_ptr<Spline>( arg_in_1 );
    Utils::mex_set_scalar_value( arg_out_0, ptr->y_max() );

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
    {"coeffs",do_coeffs},
    {"y_min_max",do_y_min_max},
    {"eval",do_eval},
    {"eval_D",do_eval_D},
    {"eval_DD",do_eval_DD},
    {"eval_DDD",do_eval_DDD},
    {"eval_DDDD",do_eval_DDDD},
    {"eval_DDDDD",do_eval_DDDDD},
    {"make_closed",do_make_closed},
    {"make_opened",do_make_opened},
    {"is_closed",do_is_closed},
    {"make_bounded",do_make_bounded},
    {"make_unbounded",do_make_unbounded},
    {"is_bounded",do_is_bounded},
    {"make_extended_constant",do_make_extended_constant},
    {"make_extended_not_constant",do_make_extended_not_constant},
    {"is_extended_constant",do_is_extended_constant},
    {"x_begin",do_x_begin},
    {"y_begin",do_y_begin},
    {"x_end",do_x_end},
    {"y_end",do_y_end},
    {"x_min",do_x_min},
    {"y_min",do_y_min},
    {"x_max",do_x_max},
    {"y_max",do_y_max}
  };

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  #define MEX_ERROR_MESSAGE \
"%======================================================================%\n" \
"Spline1DMexWrapper:  Compute spline curve\n" \
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
MEX_ERROR_MESSAGE_30 "\n" \
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
      mexErrMsgTxt( fmt::format( "SplineSetMexWrapper Error: {}", e.what() ).c_str() );
    } catch (...) {
      mexErrMsgTxt("SplineSetMexWrapper failed\n");
    }

  }

}
