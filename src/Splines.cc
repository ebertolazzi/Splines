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
 |      Università degli Studi di Trento                                    |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/
#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wpoison-system-directories"
#pragma clang diagnostic ignored "-Wundefined-func-template"
#endif

#include "SplinesUtils.hh"
#include "Utils_fmt.hh"

#include <cmath>
#include <limits> // std::numeric_limits

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

  using std::copy_n;

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
    integer               /* dim */,
    integer   const npts,
    real_type const [] /* pnts    */,
    integer            /* ld_pnts */,
    real_type       t[]
  ) {
    t[0]      = 0;
    t[npts-1] = 1;
    for ( integer k{1}; k < npts-1; ++k )
      t[k] = static_cast<real_type>(k)/static_cast<real_type>(npts);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  chordal(
    integer   const dim,
    integer   const npts,
    real_type const pnts[],
    integer   const ld_pnts,
    real_type       t[]
  ) {
    t[0] = 0;
    real_type const * p0{pnts};
    for ( integer k{1}; k < npts; ++k ) {
      real_type const * p1{p0 + ld_pnts};
      real_type dst = 0;
      for ( integer j{0}; j < dim; ++j ) {
        real_type const c{p1[j] - p0[j]};
        dst += c*c;
      }
      t[k] = t[k-1] + sqrt(dst);
    }
    for ( integer k{1}; k < npts-1; ++k ) t[k] /= t[npts-1];
    t[npts-1] = 1;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  centripetal(
    integer   const dim,
    integer   const npts,
    real_type const pnts[],
    integer   const ld_pnts,
    real_type const alpha,
    real_type       t[]
  ) {
    t[0] = 0;
    real_type const * p0{pnts};
    for ( integer k{1}; k < npts; ++k ) {
      real_type const * p1{p0 + ld_pnts};
      real_type dst{0};
      for ( integer j{0}; j < dim; ++j ) {
        real_type const c{p1[j] - p0[j]};
        dst += c*c;
      }
      t[k] = t[k-1] + pow(dst,alpha/2);
    }
    for ( integer k{1}; k < npts-1; ++k ) t[k] /= t[npts-1];
    t[npts-1] = 1;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  #if 0

  void
  universal(
    integer         dim,
    integer         npts,
    real_type const pnts[],
    integer         ld_pnts,
    real_type       t[]
  ); // to be done

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  FoleyNielsen(
    integer         dim,
    integer         npts,
    real_type const pnts[],
    integer         ld_pnts,
    real_type       t[]
  ); // to be done

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  FangHung(
    integer         dim,
    integer         npts,
    real_type const pnts[],
    integer         ld_pnts,
    real_type       t[]
  ); // to be done

  #endif

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  #ifndef DOXYGEN_SHOULD_SKIP_THIS

  SplineType1D
  string_to_splineType1D( string_view nin ) {
    string n{nin};
    std::transform(n.begin(), n.end(), n.begin(), ::tolower);
    if ( n == "constant" )   return SplineType1D::CONSTANT;
    if ( n == "linear" )     return SplineType1D::LINEAR;
    if ( n == "cubic" )      return SplineType1D::CUBIC;
    if ( n == "akima" )      return SplineType1D::AKIMA;
    if ( n == "bessel" )     return SplineType1D::BESSEL;
    if ( n == "pchip" )      return SplineType1D::PCHIP;
    if ( n == "quintic" )    return SplineType1D::QUINTIC;
    if ( n == "hermite" )    return SplineType1D::HERMITE;
    if ( n == "spline_set" ) return SplineType1D::SPLINE_SET;
    if ( n == "spline_vec" ) return SplineType1D::SPLINE_VEC;
    throw std::runtime_error(fmt::format( "string_to_splineType1D({}) unknown type\n", n ));
  }

  SplineType2D
  string_to_splineType2D( string_view nin ) {
    string n{nin};
    std::transform(n.begin(), n.end(), n.begin(), ::tolower);
    if ( n == "bilinear"  ) return SplineType2D::BILINEAR;
    if ( n == "bicubic"   ) return SplineType2D::BICUBIC;
    if ( n == "biquintic" ) return SplineType2D::BIQUINTIC;
    if ( n == "akima"     ) return SplineType2D::AKIMA2D;
    throw std::runtime_error(fmt::format( "string_to_splineType2D({}) unknown type\n", n ));
  }

  char const *
  to_string( SplineType1D const t ) {
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
  }

  char const *
  to_string( SplineType2D const t ) {
    switch( t ) {
      case SplineType2D::BILINEAR:  return "SPLINE2D_BILINEAR";
      case SplineType2D::BICUBIC:   return "SPLINE2D_BICUBIC";
      case SplineType2D::BIQUINTIC: return "SPLINE2D_BIQUINTIC";
      case SplineType2D::AKIMA2D:   return "SPLINE2D_AKIMA2D";
    }
    return "NO_TYPE";
  }

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

  void
  Spline::search( std::pair<integer,real_type> & res ) const {
    UTILS_ASSERT( m_npts > 0, "in Spline[{}]::search(...), npts == 0!", m_name );
    #ifdef SPLINES_USE_THREADS
    std::unique_lock<std::mutex> lock(m_last_interval_mutex);
    auto id = std::this_thread::get_id();
    auto it = m_last_interval.find(id);
    if ( it == m_last_interval.end() ) {
      it = m_last_interval.insert( {id,std::make_shared<integer>()} ).first;
      *it->second = 0;
    }
    integer & last_interval{ *it->second };
    lock.unlock();
    #else
    integer & last_interval = m_last_interval;
    #endif
    Utils::search_interval(
      m_npts,
      m_X,
      res.second,
      last_interval,
      m_curve_is_closed,
      m_curve_can_extend
    );
    res.first = last_interval;
  }

  #ifndef DOXYGEN_SHOULD_SKIP_THIS

  void
  Spline::init_last_interval() const {
    #ifdef SPLINES_USE_THREADS
    std::unique_lock<std::mutex> lock(m_last_interval_mutex);
    auto id = std::this_thread::get_id();
    auto it = m_last_interval.find(id);
    if ( it == m_last_interval.end() ) it = m_last_interval.insert( {id,std::make_shared<integer>()} ).first;
    integer & last_interval{ *it->second };
    #else
    integer & last_interval{ m_last_interval };
    #endif
    last_interval = 0;
  }

  integer
  SplineSurf::search_x( real_type & x ) const {
    #ifdef SPLINES_USE_THREADS
    std::unique_lock<std::mutex> lock(m_last_interval_x_mutex);
    auto id = std::this_thread::get_id();
    auto it = m_last_interval_x.find(id);
    if ( it == m_last_interval_x.end() ) {
      it = m_last_interval_x.insert( {id,std::make_shared<integer>()} ).first;
      *it->second = 0;
    }
    integer & last_interval{ *it->second };
    lock.unlock();
    #else
    integer & last_interval{ m_last_interval_x };
    #endif
    Utils::search_interval(
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
  SplineSurf::init_last_interval_x() const {
    #ifdef SPLINES_USE_THREADS
    std::unique_lock<std::mutex> lock(m_last_interval_x_mutex);
    auto id = std::this_thread::get_id();
    auto it = m_last_interval_x.find(id);
    if ( it == m_last_interval_x.end() ) it = m_last_interval_x.insert( {id,std::make_shared<integer>()} ).first;
    integer & last_interval{ *it->second };
    #else
    integer & last_interval{ m_last_interval_x };
    #endif
    last_interval = 0;
  }

  integer
  SplineSurf::search_y( real_type & y ) const {
    #ifdef SPLINES_USE_THREADS
    std::unique_lock<std::mutex> lock(m_last_interval_y_mutex);
    auto id = std::this_thread::get_id();
    auto it = m_last_interval_y.find(id);
    if ( it == m_last_interval_y.end() ) {
      it = m_last_interval_y.insert( {id,std::make_shared<integer>()} ).first;
      *it->second = 0;
    }
    integer & last_interval{ *it->second };
    lock.unlock();
    #else
    integer & last_interval{ m_last_interval_y };
    #endif
    Utils::search_interval(
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
  SplineSurf::init_last_interval_y() const {
    #ifdef SPLINES_USE_THREADS
    std::unique_lock<std::mutex> lock(m_last_interval_y_mutex);
    auto id = std::this_thread::get_id();
    auto it = m_last_interval_y.find(id);
    if ( it == m_last_interval_y.end() ) it = m_last_interval_y.insert( {id,std::make_shared<integer>()} ).first;
    integer & last_interval{ *it->second };
    #else
    integer & last_interval{ m_last_interval_y };
    #endif
    last_interval = 0;
  }

  #endif

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  //! Check if cubic spline with this data is monotone, return -1 no, 0 yes, 1 strictly monotone
  integer
  check_cubic_spline_monotonicity(
    real_type const X[],
    real_type const Y[],
    real_type const Yp[],
    integer   const npts
  ) {
    // check monotonicity of data: (assuming X monotone)
    integer flag{1};
    for ( integer i{1}; i < npts; ++i ) {
      if ( Y[i-1] > Y[i] ) return -2; // non monotone data
      if ( Utils::is_zero(Y[i-1]-Y[i]) && X[i-1] < X[i] ) flag = 0; // non strict monotone
    }
    // pag 146 Methods of Shape-Preserving Spline Approximation, K
    for ( integer i{1}; i < npts; ++i ) {
      if ( X[i] <= X[i-1] ) continue; // skip duplicate points
      real_type const dd = (Y[i]-Y[i-1])/(X[i]-X[i-1]);
      real_type const m0 = Yp[i-1]/dd;
      real_type const m1 = Yp[i]/dd;
      if ( m0 < 0 || m1 < 0 ) return -1; // non monotone
      if ( m0 <= 3 && m1 <= 3 ) {
        if ( flag > 0 && i > 1 &&
             (Utils::is_zero(m0) || Utils::is_zero(m0-3) ) ) flag = 0;
        if ( flag > 0 && i < npts-1 &&
             (Utils::is_zero(m1) || Utils::is_zero(m1-3) ) ) flag = 0;
      } else {
        real_type const tmp1 = 2*m0+m1-3;
        real_type const tmp2 = 2*(m0+m1-2);
        real_type const tmp3 = m0*tmp2-(tmp1*tmp1);
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
    real_type const x[], integer const incx,
    real_type const y[], integer const incy,
    integer const n
  ) {
    reserve( n );
    for ( integer i{0}; i < n; ++i ) m_X[i] = x[i*incx];
    for ( integer i{0}; i < n; ++i ) m_Y[i] = y[i*incy];
    m_npts = n;
    build();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Spline::setup( string const & file_name ) {
    GenericContainer gc;
    UTILS_ASSERT( gc.from_file( file_name ), "Spline::setup( '{}' ) failed to read\n", file_name );
    setup( gc );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  string
  Spline::info() const {
    string res { fmt::format( "Spline `{}` of type: {} of order: {}", m_name, type_name(), order() ) };
    if ( m_npts > 0 )
      res += fmt::format(
        "\nx_min={:.5} x_max={:.5} y_min={:.5} y_max={:.5}",
        x_min(), x_max(), y_min(), y_max()
      );
    return res;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Spline::push_back( real_type const x, real_type const y ) {
    if ( m_npts > 0 ) {
      UTILS_ASSERT(
        x >= m_X[m_npts-1], // ammetto punti doppi
        "Spline[{}]::push_back, non monotone insert at insert N.{}"
        "\nX[{}] = {:.5}\nX[{}] = {:>5}\n",
        m_name, m_npts, m_npts-1, m_X[m_npts-1], m_npts, x
      );
    }
    if ( m_npts_reserved == 0 ) {
      reserve( 2 );
    } else if ( m_npts >= m_npts_reserved ) {
      // riallocazione & copia
      integer const saved_npts{m_npts}; // salvo npts perche reserve lo azzera
      Malloc_real mem("Spline::push_back");
      mem.allocate( 2*m_npts );
      real_type * Xsaved{ mem( m_npts ) };
      real_type * Ysaved{ mem( m_npts ) };

      copy_n( m_X, m_npts, Xsaved );
      copy_n( m_Y, m_npts, Ysaved );
      reserve( (m_npts+1) * 2 );
      m_npts = saved_npts;
      copy_n( Xsaved, m_npts, m_X );
      copy_n( Ysaved, m_npts, m_Y );
    }
    m_X[m_npts] = x;
    m_Y[m_npts] = y;
    ++m_npts;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Spline::set_origin( real_type const x0 ) const {
    real_type const Tx{x0 - m_X[0]};
    real_type * ix{m_X};
    while ( ix < m_X+m_npts ) *ix++ += Tx;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Spline::set_range( real_type xmin, real_type xmax ) {
    UTILS_ASSERT(
      xmax > xmin,
      "Spline[{}]::set_range({},{}) bad range ", m_name, xmin, xmax
    );
    real_type const S  = (xmax - xmin) / ( m_X[m_npts-1] - m_X[0] );
    real_type const Tx = xmin - S * m_X[0];
    for( real_type *ix = m_X; ix < m_X+m_npts; ++ix ) *ix = *ix * S + Tx;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Spline::dump(
    ostream_type &    s,
    integer     const nintervals,
    string_view const header
  ) const {
    s << header << '\n';
    real_type const dx{ (x_max()-x_min())/nintervals };
    for ( integer i{0}; i <= nintervals; ++i ) {
      real_type x{ x_min() + i*dx };
      fmt::print( s, "{}\t{}\n", x, this->eval(x) );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Spline::y_min_max(
    integer   & i_min_pos,
    real_type & x_min_pos,
    real_type & y_min,
    integer   & i_max_pos,
    real_type & x_max_pos,
    real_type & y_max
  ) const {
    i_min_pos = i_max_pos = 0;
    x_min_pos = y_min = x_max_pos = y_max = 0;
    UTILS_ERROR(
      "In spline: {} y_min_max not implemented\n",
      info()
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Spline::y_min_max(
    vector<integer>   & i_min_pos,
    vector<real_type> & x_min_pos,
    vector<real_type> & y_min,
    vector<integer>   & i_max_pos,
    vector<real_type> & x_max_pos,
    vector<real_type> & y_max
  ) const {
    i_min_pos.clear();
    i_max_pos.clear();
    x_min_pos.clear();
    x_max_pos.clear();
    y_min.clear();
    y_max.clear();
    UTILS_ERROR(
      "In spline: {} y_min_max not implemented\n",
      info()
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  curvature( real_type const s, Spline const & X, Spline const & Y ) {
    real_type const x_1 { X.D(s)    };
    real_type const x_2 { X.DD(s)   };
    real_type const y_1 { Y.D(s)    };
    real_type const y_2 { Y.DD(s)   };
    real_type const t4  { x_1 * x_1 };
    real_type const t5  { y_1 * y_1 };
    real_type const t6  { t4 + t5   };
    real_type const t7  { sqrt(t6)  };
    return (x_1 * y_2 - y_1 * x_2) / ( t6 * t7 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  curvature_D( real_type const s, Spline const & X, Spline const & Y ) {
    real_type const x_1 { X.D(s)    };
    real_type const x_2 { X.DD(s)   };
    real_type const x_3 { X.DDD(s)  };
    real_type const y_1 { Y.D(s)    };
    real_type const y_2 { Y.DD(s)   };
    real_type const y_3 { Y.DDD(s)  };
    real_type const t1  { x_1 * x_1 };
    real_type const t9  { x_2 * x_2 };
    real_type const t13 { y_1 * y_1 };
    real_type const t17 { y_2 * y_2 };
    real_type const t26 { t1 + t13  };
    real_type const t27 { t26 * t26 };
    real_type const t28 { sqrt(t26) };
    real_type const aa  { y_3 * x_1 - y_1 * x_3 };
    real_type const bb  { 3 * y_2 * x_2 };
    return ( t1 * ( aa - bb ) + t13 * ( aa + bb )
             + 3 * x_1 * y_1 * ( t9 - t17 ) ) / ( t28 * t27 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  curvature_DD( real_type const s, Spline const & X, Spline const & Y ) {
    real_type const x_1 { X.D(s)    };
    real_type const x_2 { X.DD(s)   };
    real_type const x_3 { X.DDD(s)  };
    real_type const x_4 { X.DDDD(s) };
    real_type const y_1 { Y.D(s)    };
    real_type const y_2 { Y.DD(s)   };
    real_type const y_3 { Y.DDD(s)  };
    real_type const y_4 { Y.DDDD(s) };

    real_type const t1  { y_1 * y_1 };
    real_type const t2  { t1 * t1   };
    real_type const t12 { x_2 * x_2 };
    real_type const t13 { t12 * x_2 };
    real_type const t15 { x_1 * x_3 };
    real_type const t16 { 9 * t15   };
    real_type const t17 { y_2 * y_2 };
    real_type const t21 { x_1 * x_1 };
    real_type const t22 { x_4 * t21 };
    real_type const t26 { 9 * x_1 * y_2 * y_3 };
    real_type const t30 { t21 * x_1 };
    real_type const t40 { t17 * y_2 };
    real_type const t64 { t21 + t1  };
    real_type const t65 { t64 * t64 };
    real_type const t67 { sqrt(t64) };
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
    string const where{ fmt::format("Spline[{}]::setup( gc ):", m_name ) };
    GenericContainer const & gc_x{ gc("xdata",where) };
    GenericContainer const & gc_y{ gc("ydata",where) };

    vec_real_type x, y;
    {
      string const ff{ fmt::format( "{}, field `xdata'", where ) };
      gc_x.copyto_vec_real( x, ff );
    }
    {
      string const ff{ fmt::format( "{}, field `ydata'", where ) };
      gc_y.copyto_vec_real ( y, ff );
    }
    build( x, y );
  }

}
