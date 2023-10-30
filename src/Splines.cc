/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2016                                                      |
 |                                                                          |
 |         , __                 , __                                        |
 |        /|/  \               /|/  \                                       |
 |         | __/ _   ,_         | __/ _   ,_                                |
 |         |   \|/  /  |  |   | |   \|/  /  |  |   |                        |
 |         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                       |
 |                           /|                   /|                        |
 |                           \|                   \|                        |
 |                                                                          |
 |      Enrico Bertolazzi                                                   |
 |      Dipartimento di Ingegneria Industriale                              |
 |      Universita` degli Studi di Trento                                   |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#include "SplinesUtils.hh"

#include <cmath>
#include <limits> // std::numeric_limits

#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
#pragma clang diagnostic ignored "-Wpoison-system-directories"
#pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
#endif

#ifdef SPLINES_OS_OSX
  #define UNW_LOCAL_ONLY
  #include <cxxabi.h>
  #include <libunwind.h>
#endif

namespace Splines {

  using std::abs;
  using std::max;
  using std::min;

  // cbrt is not available on WINDOWS? or C++ < C++11?
  #ifdef _MSC_VER
    using std::sqrt;
    using std::pow;
    static inline real_type cbrt( real_type x ) { return pow( x, 1.0/3.0 ); }
  #else
    using std::sqrt;
    using std::pow;
    using std::cbrt;
  #endif

  /*\
   |  ____                                _        _          _   _
   | |  _ \ __ _ _ __ __ _ _ __ ___   ___| |_ _ __(_)______ _| |_(_) ___  _ __
   | | |_) / _` | '__/ _` | '_ ` _ \ / _ \ __| '__| |_  / _` | __| |/ _ \| '_ \
   | |  __/ (_| | | | (_| | | | | | |  __/ |_| |  | |/ / (_| | |_| | (_) | | | |
   | |_|   \__,_|_|  \__,_|_| |_| |_|\___|\__|_|  |_/___\__,_|\__|_|\___/|_| |_|
   |
  \*/

  void
  uniform(
    integer           /* dim */,
    integer           npts,
    real_type const * /* pnts    */,
    integer           /* ld_pnts */,
    real_type       * t
  ) {
    t[0]      = 0;
    t[npts-1] = 1;
    for ( integer k = 1; k < npts-1; ++k )
      t[k] = static_cast<real_type>(k)/static_cast<real_type>(npts);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  chordal(
    integer           dim,
    integer           npts,
    real_type const * pnts,
    integer           ld_pnts,
    real_type       * t
  ) {
    t[0] = 0;
    real_type const * p0 = pnts;
    for ( integer k = 1; k < npts; ++k ) {
      real_type const * p1 = p0 + ld_pnts;
      real_type dst = 0;
      for ( integer j = 0; j < dim; ++j ) {
        real_type c = p1[j] - p0[j];
        dst += c*c;
      }
      t[k] = t[k-1] + sqrt(dst);
    }
    for ( integer k = 1; k < npts-1; ++k ) t[k] /= t[npts-1];
    t[npts-1] = 1;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  centripetal(
    integer           dim,
    integer           npts,
    real_type const * pnts,
    integer           ld_pnts,
    real_type         alpha,
    real_type       * t
  ) {
    t[0] = 0;
    real_type const * p0 = pnts;
    for ( integer k = 1; k < npts; ++k ) {
      real_type const * p1 = p0 + ld_pnts;
      real_type dst = 0;
      for ( integer j = 0; j < dim; ++j ) {
        real_type c = p1[j] - p0[j];
        dst += c*c;
      }
      t[k] = t[k-1] + pow(dst,alpha/2);
    }
    for ( integer k = 1; k < npts-1; ++k ) t[k] /= t[npts-1];
    t[npts-1] = 1;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  #if 0

  void
  universal(
    integer           dim,
    integer           npts,
    real_type const * pnts,
    integer           ld_pnts,
    real_type       * t
  ); // to be done

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  FoleyNielsen(
    integer           dim,
    integer           npts,
    real_type const * pnts,
    integer           ld_pnts,
    real_type       * t
  ); // to be done

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  FangHung(
    integer           dim,
    integer           npts,
    real_type const * pnts,
    integer           ld_pnts,
    real_type       * t
  ); // to be done

  #endif

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  #ifndef DOXYGEN_SHOULD_SKIP_THIS

  SplineType1D
  string_to_splineType1D( std::string const & nin ) {
    std::string n = nin;
    std::transform(n.begin(), n.end(), n.begin(), ::tolower);
    if      ( n == "constant" )   return SplineType1D::CONSTANT;
    else if ( n == "linear" )     return SplineType1D::LINEAR;
    else if ( n == "cubic" )      return SplineType1D::CUBIC;
    else if ( n == "akima" )      return SplineType1D::AKIMA;
    else if ( n == "bessel" )     return SplineType1D::BESSEL;
    else if ( n == "pchip" )      return SplineType1D::PCHIP;
    else if ( n == "quintic" )    return SplineType1D::QUINTIC;
    else if ( n == "hermite" )    return SplineType1D::HERMITE;
    else if ( n == "spline_set" ) return SplineType1D::SPLINE_SET;
    else if ( n == "spline_vec" ) return SplineType1D::SPLINE_VEC;
    throw std::runtime_error(fmt::format( "string_to_splineType1D({}) unknown type\n", n ));
  }

  SplineType2D
  string_to_splineType2D( std::string const & nin ) {
    std::string n = nin;
    std::transform(n.begin(), n.end(), n.begin(), ::tolower);
    if      ( n == "bilinear"  ) return SplineType2D::BILINEAR;
    else if ( n == "bicubic"   ) return SplineType2D::BICUBIC;
    else if ( n == "biquintic" ) return SplineType2D::BIQUINTIC;
    else if ( n == "akima"     ) return SplineType2D::AKIMA2D;
    throw std::runtime_error(fmt::format( "string_to_splineType2D({}) unknown type\n", n ));
  }

  char const *
  to_string( SplineType1D t ) {
    switch( t ) {
      case SplineType1D::CONSTANT:   return "SPLINE_CONSTANT";
      case SplineType1D::LINEAR:     return "SPLINE_LINEAR";
      case SplineType1D::CUBIC:      return "SPLINE_CUBIC";
      case SplineType1D::AKIMA:      return "SPLINE_AKIMA";
      case SplineType1D::BESSEL:     return "SPLINE_BESSEL";
      case SplineType1D::PCHIP:      return "SPLINE_PCHIP";
      case SplineType1D::QUINTIC:    return "SPLINE_QUINTIC";
      case SplineType1D::HERMITE:    return "SPLINE_HERMITE";
      case SplineType1D::SPLINE_SET: return "SPLINE_SPLINE_SET";
      case SplineType1D::SPLINE_VEC: return "SPLINE_SPLINE_VEC";
    }
    return "NO_TYPE";
  };

  char const *
  to_string( SplineType2D t ) {
    switch( t ) {
      case SplineType2D::BILINEAR:  return "SPLINE2D_BILINEAR";
      case SplineType2D::BICUBIC:   return "SPLINE2D_BICUBIC";
      case SplineType2D::BIQUINTIC: return "SPLINE2D_BIQUINTIC";
      case SplineType2D::AKIMA2D:   return "SPLINE2D_AKIMA2D";
    }
    return "NO_TYPE";
  };

  #endif

  /*\
   |   ____        _ _
   |  / ___| _ __ | (_)_ __   ___
   |  \___ \| '_ \| | | '_ \ / _ \
   |   ___) | |_) | | | | | |  __/
   |  |____/| .__/|_|_|_| |_|\___|
   |        |_|
  \*/

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer
  Spline::search( real_type & x ) const {
    UTILS_ASSERT( m_npts > 0, "in Spline[{}]::search(...), npts == 0!", m_name );
    #ifdef SPLINES_USE_THREADS
    bool ok{true};
    integer & last_interval = *m_last_interval.search( std::this_thread::get_id(), ok );
    if ( !ok ) last_interval = 0;
    #else
    integer & last_interval = m_last_interval;
    #endif
    Utils::searchInterval(
      m_npts,
      m_X,
      x,
      last_interval,
      m_curve_is_closed,
      m_curve_can_extend
    );
    return last_interval;
  }

  void
  Spline::init_last_interval() {
    #ifdef SPLINES_USE_THREADS
    bool ok;
    integer & last_interval = *m_last_interval.search( std::this_thread::get_id(), ok );
    #else
    integer & last_interval = m_last_interval;
    #endif
    last_interval = 0;
  }

  integer
  SplineSurf::search_x( real_type & x ) const {
    #ifdef SPLINES_USE_THREADS
    bool ok{true};
    integer & last_interval = *m_last_interval_x.search( std::this_thread::get_id(), ok );
    if ( !ok ) last_interval = 0;
    #else
    integer & last_interval = m_last_interval_x;
    #endif
    Utils::searchInterval(
      m_nx,
      m_X,
      x,
      last_interval,
      m_x_closed,
      m_x_can_extend
    );
    return last_interval;
  }

  void
  SplineSurf::init_last_interval_x() {
    #ifdef SPLINES_USE_THREADS
    bool ok;
    integer & last_interval = *m_last_interval_x.search( std::this_thread::get_id(), ok );
    #else
    integer & last_interval = m_last_interval_x;
    #endif
    last_interval = 0;
  }

  integer
  SplineSurf::search_y( real_type & y ) const {
    #ifdef SPLINES_USE_THREADS
    bool ok{true};
    integer & last_interval = *m_last_interval_y.search( std::this_thread::get_id(), ok );
    if ( !ok ) last_interval = 0;
    #else
    integer & last_interval = m_last_interval_y;
    #endif
    Utils::searchInterval(
      m_ny,
      m_Y,
      y,
      last_interval,
      m_y_closed,
      m_y_can_extend
    );
    return last_interval;
  }

  void
  SplineSurf::init_last_interval_y() {
    #ifdef SPLINES_USE_THREADS
    bool ok;
    integer & last_interval = *m_last_interval_y.search( std::this_thread::get_id(), ok );
    #else
    integer & last_interval = m_last_interval_y;
    #endif
    last_interval = 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  //! Check if cubic spline with this data is monotone, return -1 no, 0 yes, 1 strictly monotone
  integer
  checkCubicSplineMonotonicity(
    real_type const * X,
    real_type const * Y,
    real_type const * Yp,
    integer           npts
  ) {
    // check monotonicity of data: (assuming X monotone)
    integer flag = 1;
    for ( size_t i = 1; i < size_t(npts); ++i ) {
      if ( Y[i-1] > Y[i] ) return -2; // non monotone data
      if ( Utils::is_zero(Y[i-1]-Y[i]) && X[i-1] < X[i] ) flag = 0; // non strict monotone
    }
    // pag 146 Methods of Shape-Preserving Spline Approximation, K
    for ( size_t i = 1; i < size_t(npts); ++i ) {
      if ( X[i] <= X[i-1] ) continue; // skip duplicate points
      real_type dd = (Y[i]-Y[i-1])/(X[i]-X[i-1]);
      real_type m0 = Yp[i-1]/dd;
      real_type m1 = Yp[i]/dd;
      if ( m0 < 0 || m1 < 0 ) return -1; // non monotone
      if ( m0 <= 3 && m1 <= 3 ) {
        if ( flag > 0 && i > 1 &&
             (Utils::is_zero(m0) || Utils::is_zero(m0-3) ) ) flag = 0;
        if ( flag > 0 && i < size_t(npts-1) &&
             (Utils::is_zero(m1) || Utils::is_zero(m1-3) ) ) flag = 0;
      } else {
        real_type tmp1 = 2*m0+m1-3;
        real_type tmp2 = 2*(m0+m1-2);
        real_type tmp3 = m0*tmp2-(tmp1*tmp1);
        if ( tmp2 >= 0 ) {
          if ( tmp3 < 0 ) return -1; // non monotone spline
        } else {
          if ( tmp3 > 0 ) return -1;
        }
        if ( Utils::is_zero(tmp3) ) flag = 0;
      }
    }
    return flag; // passed all check
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Spline::build(
    real_type const * x, integer incx,
    real_type const * y, integer incy,
    integer n
  ) {
    reserve( n );
    for ( integer i = 0; i < n; ++i ) m_X[i] = x[i*incx];
    for ( integer i = 0; i < n; ++i ) m_Y[i] = y[i*incy];
    m_npts = n;
    build();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  string
  Spline::info() const {
    string res = fmt::format(
      "Spline `{}` of type: {} of order: {}",
      m_name, type_name(), order()
    );
    if ( m_npts > 0 )
      res += fmt::format(
        "\nx_min = {} x_max = {} y_min = {} y_max = {}",
        x_min(), x_max(), y_min(), y_max()
      );
    return res;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Spline::push_back( real_type x, real_type y ) {
    if ( m_npts > 0 ) {
      UTILS_ASSERT(
        x >= m_X[size_t(m_npts-1)], // ammetto punti doppi
        "Spline[{}]::push_back, non monotone insert at insert N.{}"
        "\nX[{}] = {}\nX[{}] = {}\n",
        m_name, m_npts, m_npts-1, m_X[size_t(m_npts-1)], m_npts, x
      );
    }
    if ( m_npts_reserved == 0 ) {
      reserve( 2 );
    } else if ( m_npts >= m_npts_reserved ) {
      // riallocazione & copia
      integer saved_npts = m_npts; // salvo npts perche reserve lo azzera
      Malloc_real mem("Spline::push_back");
      mem.allocate( size_t(2*m_npts) );
      real_type * Xsaved = mem( size_t(m_npts) );
      real_type * Ysaved = mem( size_t(m_npts) );

      std::copy_n( m_X, m_npts, Xsaved );
      std::copy_n( m_Y, m_npts, Ysaved );
      reserve( (m_npts+1) * 2 );
      m_npts = saved_npts;
      std::copy_n( Xsaved, m_npts, m_X );
      std::copy_n( Ysaved, m_npts, m_Y );
    }
    m_X[m_npts] = x;
    m_Y[m_npts] = y;
    ++m_npts;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Spline::set_origin( real_type x0 ) {
    real_type Tx = x0 - m_X[0];
    real_type *ix = m_X;
    while ( ix < m_X+m_npts ) *ix++ += Tx;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Spline::set_range( real_type xmin, real_type xmax ) {
    UTILS_ASSERT(
      xmax > xmin,
      "Spline[{}]::set_range({},{}) bad range ", m_name, xmin, xmax
    );
    real_type S  = (xmax - xmin) / ( m_X[m_npts-1] - m_X[0] );
    real_type Tx = xmin - S * m_X[0];
    for( real_type *ix = m_X; ix < m_X+m_npts; ++ix ) *ix = *ix * S + Tx;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Spline::dump(
    ostream_type & s,
    integer        nintervals,
    char const *   header
  ) const {
    s << header << '\n';
    real_type dx = (x_max()-x_min())/nintervals;
    for ( integer i = 0; i <= nintervals; ++i ) {
      real_type x = x_min() + i*dx;
      fmt::print( s, "{}\t{}\n", x, (*this)(x) );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  curvature( real_type s, Spline const & X, Spline const & Y ) {
    real_type x_1 = X.D(s);
    real_type x_2 = X.DD(s);
    real_type y_1 = Y.D(s);
    real_type y_2 = Y.DD(s);
    real_type t4  = x_1 * x_1;
    real_type t5  = y_1 * y_1;
    real_type t6  = t4 + t5;
    real_type t7  = sqrt(t6);
    return (x_1 * y_2 - y_1 * x_2) / ( t6 * t7 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  curvature_D( real_type s, Spline const & X, Spline const & Y ) {
    real_type x_1 = X.D(s);
    real_type x_2 = X.DD(s);
    real_type x_3 = X.DDD(s);
    real_type y_1 = Y.D(s);
    real_type y_2 = Y.DD(s);
    real_type y_3 = Y.DDD(s);
    real_type t1  = x_1 * x_1;
    real_type t9  = x_2 * x_2;
    real_type t13 = y_1 * y_1;
    real_type t17 = y_2 * y_2;
    real_type t26 = t1 + t13;
    real_type t27 = t26 * t26;
    real_type t28 = sqrt(t26);
    real_type aa  = y_3 * x_1 - y_1 * x_3;
    real_type bb  = 3 * y_2 * x_2;
    return ( t1 * ( aa - bb ) + t13 * ( aa + bb )
             + 3 * x_1 * y_1 * ( t9 - t17 ) ) / ( t28 * t27 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  curvature_DD( real_type s, Spline const & X, Spline const & Y ) {
    real_type x_1 = X.D(s);
    real_type x_2 = X.DD(s);
    real_type x_3 = X.DDD(s);
    real_type x_4 = X.DDDD(s);
    real_type y_1 = Y.D(s);
    real_type y_2 = Y.DD(s);
    real_type y_3 = Y.DDD(s);
    real_type y_4 = Y.DDDD(s);

    real_type t1  = y_1 * y_1;
    real_type t2  = t1 * t1;
    real_type t12 = x_2 * x_2;
    real_type t13 = t12 * x_2;
    real_type t15 = x_1 * x_3;
    real_type t16 = 9 * t15;
    real_type t17 = y_2 * y_2;
    real_type t21 = x_1 * x_1;
    real_type t22 = x_4 * t21;
    real_type t26 = 9 * x_1 * y_2 * y_3;
    real_type t30 = t21 * x_1;
    real_type t40 = t17 * y_2;
    real_type t64 = t21 + t1;
    real_type t65 = t64 * t64;
    real_type t67 = sqrt(t64);
    return (
      t2 * (x_1 * y_4 + 4 * x_2 * y_3 + 5 * x_3 * y_2 - y_1 * x_4)
      + t1 * y_1 * (3 * t13 + x_2 * (t16 - 12 * t17) - 2 * t22 - t26)
      + t1 * ( x_1 * ( 12 * t40 - 33 * y_2 * t12 )
               + t21 * ( y_2 * x_3 - y_3 * x_2)
               + 2 * y_4 * t30 )
      - y_1 * t21 * ( 12 * t13 - x_2 * (t16+33*t17) + t22 + t26 )
      + t30 * ( y_2 * (12 * t12 - 4 * t15) + y_4 * t21
                -5 * x_1 * x_2 * y_3 - 3 * t40)
    ) / (t67*t65*t64);
  }

  /*
  //    ____  ____   ____                               _
  //   / ___|/ ___| / ___| _   _ _ __  _ __   ___  _ __| |_
  //  | |  _| |     \___ \| | | | '_ \| '_ \ / _ \| '__| __|
  //  | |_| | |___   ___) | |_| | |_) | |_) | (_) | |  | |_
  //   \____|\____| |____/ \__,_| .__/| .__/ \___/|_|   \__|
  //                            |_|   |_|
  */

  using GC_namespace::GC_type;
  using GC_namespace::vec_real_type;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Spline::setup( GenericContainer const & gc ) {
    /*
    // gc["xdata"]
    // gc["ydata"]
    //
    */
    string where = fmt::format("Spline[{}]::setup( gc ):", m_name );
    GenericContainer const & gc_x = gc("xdata",where.c_str());
    GenericContainer const & gc_y = gc("ydata",where.c_str());

    vec_real_type x, y;
    {
      std::string ff = fmt::format( "{}, field `xdata'", where );
      gc_x.copyto_vec_real( x, ff.c_str() );
    }
    {
      std::string ff = fmt::format( "{}, field `ydata'", where );
      gc_y.copyto_vec_real ( y, ff.c_str() );
    }
    build( x, y );
  }

}
