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

#pragma once

#ifndef SPLINES_HH
#define SPLINES_HH

#ifdef __GNUC__
#pragma GCC diagnostic push
#endif

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
#pragma clang diagnostic ignored "-Wpoison-system-directories"
#endif

#include "SplinesConfig.hh"
#include <fstream>

//!
//! Namespace of Splines library
//!
namespace Splines {

  using std::vector;
  using std::string;
  using std::string_view;
  using std::exception;
  using std::runtime_error;
  using std::basic_ostream;
  using std::basic_istream;
  using std::lower_bound;
  using std::pair;
  using std::cout;
  using std::cin;
  using std::cerr;

  typedef double real_type; //!< Floating point type for splines
  typedef int    integer;   //!< Signed integer type for splines

  using Malloc_real  = Utils::Malloc<real_type>;
  using ostream_type = basic_ostream<char>;
  using istream_type = basic_istream<char>;

  void backtrace( ostream_type & );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  //! Associate a number for each type of splines implemented
  using SplineType1D = enum class SplineType1D : integer {
    CONSTANT   = 0,
    LINEAR     = 1,
    CUBIC      = 2,
    AKIMA      = 3,
    BESSEL     = 4,
    PCHIP      = 5,
    QUINTIC    = 6,
    HERMITE    = 7,
    SPLINE_SET = 8,
    SPLINE_VEC = 9
  };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  //! Associate a number for each type of splines implemented
  using SplineType2D = enum class SplineType2D : integer {
    BILINEAR  = 0,
    BICUBIC   = 1,
    BIQUINTIC = 2,
    AKIMA2D   = 3
  };

  extern SplineType1D string_to_splineType1D( string_view n );
  extern SplineType2D string_to_splineType2D( string_view n );
  extern char const * to_string( SplineType2D t );
  extern char const * to_string( SplineType1D t );

  using GC_namespace::GenericContainer;
  using GC_namespace::vec_real_type;
  using GC_namespace::vec_string_type;
  using GC_namespace::vector_type;
  using GC_namespace::map_type;

  /*
  //   _   _                     _ _
  //  | | | | ___ _ __ _ __ ___ (_) |_ ___
  //  | |_| |/ _ \ '__| '_ ` _ \| | __/ _ \
  //  |  _  |  __/ |  | | | | | | | ||  __/
  //  |_| |_|\___|_|  |_| |_| |_|_|\__\___|
  */

  #ifdef AUTODIFF_SUPPORT
  template <typename T>
  inline
  void
  Hermite3( T const & x, real_type const H, T base[4] ) {
    T const X{x/H};
    base[1] = X*X*(3-2*X);
    base[0] = 1-base[1];
    base[2] = x*(X*(X-2)+1);
    base[3] = x*X*(X-1);
  }
  #endif

  void Hermite3     ( real_type const x, real_type const H, real_type base[4] );
  void Hermite3_D   ( real_type const x, real_type const H, real_type base_D[4] );
  void Hermite3_DD  ( real_type const x, real_type const H, real_type base_DD[4] );
  void Hermite3_DDD ( real_type const x, real_type const H, real_type base_DDD[4] );

  #ifdef AUTODIFF_SUPPORT
  template <typename T>
  inline
  void
  Hermite5( T const & x, real_type const H, T base[6] ) {
    auto const t1  { H*H   };
    auto const t4  { x*x   };
    auto const t7  { H-x   };
    auto const t8  { t7*t7 };
    auto const t9  { t8*t7 };
    auto const t11 { t1*t1 };
    auto const t2  { 1/t11 };
    auto const t3  { 1/H   };
    auto const t13 { t3*t2 };
    auto const t14 { t4*x  };
    auto const t17 { t4*t4 };
    base[0] = t13*t9*(3.0*x*H+t1+6.0*t4);
    base[1] = t13*(-15.0*H*t17+6.0*t17*x+10.0*t1*t14);
    base[2] = t2*t9*x*(H+3*x);
    base[3] = t2*(3*x-4*H)*t7*t14;
    real_type const t36 = t3/t1/2;
    base[4] = t36*t9*t4;
    base[5] = t36*t8*t14;
  }
  #endif

  void Hermite5       ( real_type const x, real_type const H, real_type base[6] );
  void Hermite5_D     ( real_type const x, real_type const H, real_type base_D[6] );
  void Hermite5_DD    ( real_type const x, real_type const H, real_type base_DD[6] );
  void Hermite5_DDD   ( real_type const x, real_type const H, real_type base_DDD[6] );
  void Hermite5_DDDD  ( real_type const x, real_type const H, real_type base_DDDD[6] );
  void Hermite5_DDDDD ( real_type const x, real_type const H, real_type base_DDDDD[6] );

  //!
  //! Convert polynomial defined using Hermite base
  //!
  //! \f[ p(x) = p_0 h_0(t) + p_1 h_1(t) +  p'_0 h_2(t) + p'_1 h_3(t) \f]
  //!
  //! to standard form
  //!
  //! \f[ p(x) = A t^3 + B t^2 + C t + D \f]
  //!
  //!
  static
  inline
  void
  Hermite3_to_poly(
    real_type const H,
    real_type const P0,
    real_type const P1,
    real_type const DP0,
    real_type const DP1,
    real_type     & A,
    real_type     & B,
    real_type     & C,
    real_type     & D
  ) {
    real_type H2{ H*H };
    real_type P10{ P1-P0 };
    A = (DP0+DP1-2*P10/H)/H2;
    B = (3*P10/H-(2*DP0+DP1))/H;
    C = DP0;
    D = P0;
  }

  //!
  //! Convert polynomial defined using Hermite base
  //!
  //! \f[
  //! p(x) = p_0 h_0(t) + p_1 h_1(t) +
  //! p'_0 h_2(t) + p'_1 h_3(t) +
  //! p''_0 h_4(t) + p''_1 h_5(t)
  //! \f]
  //!
  //! to standard form
  //!
  //! \f[ p(x) = A t^5 + B t^4 + C t^3 + D t^2 + E t + F \f]
  //!
  //!
  static
  inline
  void
  Hermite5_to_poly(
    real_type const h,
    real_type const P0,
    real_type const P1,
    real_type const DP0,
    real_type const DP1,
    real_type const DDP0,
    real_type const DDP1,
    real_type     & A,
    real_type     & B,
    real_type     & C,
    real_type     & D,
    real_type     & E,
    real_type     & F
  ) {
    real_type h2{ h*h };
    real_type h3{ h*h2 };
    real_type P10{ P1-P0 };
    A = ( (DDP1-DDP0)/2+(6*P10/h-3*(DP0+DP1))/h)/h3;
    B = ( (1.5*DDP0-DDP1)+ ((8*DP0+7*DP1)-15*P10/h)/h )/h2;
    C = ( 0.5*DDP1-1.5*DDP0 + (10*P10/h -(6*DP0+4*DP1))/h )/h;
    D = DDP0/2;
    E = DP0;
    F = P0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  /*
  //   ____  _ _ _
  //  | __ )(_) (_)_ __   ___  __ _ _ __
  //  |  _ \| | | | '_ \ / _ \/ _` | '__|
  //  | |_) | | | | | | |  __/ (_| | |
  //  |____/|_|_|_|_| |_|\___|\__,_|_|
  */

  real_type
  bilinear3(
    real_type const p[4],
    real_type const M[4][4],
    real_type const q[4]
  );

  real_type
  bilinear5(
    real_type const p[6],
    real_type const M[6][6],
    real_type const q[6]
  );

  integer
  check_cubic_spline_monotonicity(
    real_type const X[],
    real_type const Y[],
    real_type const Yp[],
    integer         npts
  );

  /*\
   |  ____                                _        _          _   _
   | |  _ \ __ _ _ __ __ _ _ __ ___   ___| |_ _ __(_)______ _| |_(_) ___  _ __
   | | |_) / _` | '__/ _` | '_ ` _ \ / _ \ __| '__| |_  / _` | __| |/ _ \| '_ \
   | |  __/ (_| | | | (_| | | | | | |  __/ |_| |  | |/ / (_| | |_| | (_) | | | |
   | |_|   \__,_|_|  \__,_|_| |_| |_|\___|\__|_|  |_/___\__,_|\__|_|\___/|_| |_|
   |
  \*/

  //!
  //! Compute nodes for the spline using uniform distribution
  //!
  //! \param[in]  dim     dimension of the points
  //! \param[in]  npts    number of points
  //! \param[in]  pnts    matrix whose columns are the points
  //! \param[in]  ld_pnts leading dimension of the matrix (fortran storage)
  //! \param[out] t       vector of the computed nodes
  //!
  void
  uniform(
    integer   const dim,
    integer   const npts,
    real_type const pnts[],
    integer   const ld_pnts,
    real_type       t[]
  );

  //!
  //! Compute nodes for the spline using chordal distribution
  //!
  //! \param[in]  dim     dimension of the points
  //! \param[in]  npts    number of points
  //! \param[in]  pnts    matrix whose columns are the points
  //! \param[in]  ld_pnts leading dimension of the matrix (fortran storage)
  //! \param[out] t       vector of the computed nodes
  //!
  void
  chordal(
    integer   const dim,
    integer   const npts,
    real_type const pnts[],
    integer   const ld_pnts,
    real_type       t[]
  );

  //!
  //! Compute nodes for the spline using centripetal distribution
  //!
  //! \param[in]  dim     dimension of the points
  //! \param[in]  npts    number of points
  //! \param[in]  pnts    matrix whose columns are the points
  //! \param[in]  ld_pnts leading dimension of the matrix (fortran storage)
  //! \param[in]  alpha   power factor
  //! \param[out] t       vector of the computed nodes
  //!
  void
  centripetal(
    integer   const dim,
    integer   const npts,
    real_type const pnts[],
    integer   const ld_pnts,
    real_type const alpha,
    real_type       t[]
  );

  //!
  //! Compute nodes for the spline using universal distribution
  //!
  //! \param[in]  dim     dimension of the points
  //! \param[in]  npts    number of points
  //! \param[in]  pnts    matrix whose columns are the points
  //! \param[in]  ld_pnts leading dimension of the matrix (fortran storage)
  //! \param[out] t       vector of the computed nodes
  //!
  void
  universal(
    integer   const dim,
    integer   const npts,
    real_type const pnts[],
    integer   const ld_pnts,
    real_type       t[]
  );

  //!
  //! Compute nodes for the spline using `FoleyNielsen` distribution
  //!
  //! \param[in]  dim     dimension of the points
  //! \param[in]  npts    number of points
  //! \param[in]  pnts    matrix whose columns are the points
  //! \param[in]  ld_pnts leading dimension of the matrix (fortran storage)
  //! \param[out] t       vector of the computed nodes
  //!
  void
  FoleyNielsen(
    integer   const dim,
    integer   const npts,
    real_type const pnts[],
    integer   const ld_pnts,
    real_type       t[]
  );

  //!
  //! Compute nodes for the spline using `FangHung` distribution
  //!
  //! \param[in]  dim     dimension of the points
  //! \param[in]  npts    number of points
  //! \param[in]  pnts    matrix whose columns are the points
  //! \param[in]  ld_pnts leading dimension of the matrix (fortran storage)
  //! \param[out] t       vector of the computed nodes
  //!
  void
  FangHung(
    integer   const dim,
    integer   const npts,
    real_type const pnts[],
    integer   const ld_pnts,
    real_type       t[]
  );

  /*\
   |   ____                      _     ___       _                       _
   |  / ___|  ___  __ _ _ __ ___| |__ |_ _|_ __ | |_ ___ _ ____   ____ _| |
   |  \___ \ / _ \/ _` | '__/ __| '_ \ | || '_ \| __/ _ \ '__\ \ / / _` | |
   |   ___) |  __/ (_| | | | (__| | | || || | | | ||  __/ |   \ V / (_| | |
   |  |____/ \___|\__,_|_|  \___|_| |_|___|_| |_|\__\___|_|    \_/ \__,_|_|
  \*/

  //!
  //! Manage Search intervals
  //!
  class SearchInterval {

    static integer const m_table_size{ 400 };

    string const * p_name{nullptr};
    integer      * p_npts{nullptr};
    bool         * p_curve_is_closed{nullptr};
    bool         * p_curve_can_extend{nullptr};

    mutable real_type ** p_X{nullptr};
    mutable real_type    m_x_min{0};
    mutable real_type    m_x_max{0};
    mutable real_type    m_x_range{0};
    mutable real_type    m_dx{0};
    mutable integer      m_LO[m_table_size+2]; // to avoid overflow and replicate last point
    mutable integer      m_HI[m_table_size+2]; // to avoid overflow and replicate last point
    mutable bool         m_must_reset{ true };
    mutable std::mutex   m_mutex;
    void reset() const;

  public:

    SearchInterval( SearchInterval const & ) = delete;
    SearchInterval const & operator = ( SearchInterval const & ) = delete;

    SearchInterval() {}

    void
    setup( string const * name, integer * n, real_type ** X, bool * is_closed, bool * can_extend ) {
      p_name             = name;
      p_npts             = n;
      p_X                = X;
      p_curve_is_closed  = is_closed;
      p_curve_can_extend = can_extend;
      m_must_reset       = true;
    }

    //!
    //! Find interval containing `res.second` using binary search.
    //! Return result in `res.first`
    //!
    void find( std::pair<integer,real_type> & res ) const;
    void must_reset() { m_must_reset = true; }
  };

  /*\
   |   ____        _ _
   |  / ___| _ __ | (_)_ __   ___
   |  \___ \| '_ \| | | '_ \ / _ \
   |   ___) | |_) | | | | | |  __/
   |  |____/| .__/|_|_|_| |_|\___|
   |        |_|
  \*/
  //!
  //! Spline Management Class
  //!
  class Spline {
    friend class SplineSet;
  protected:

    string const m_name;
    bool         m_curve_is_closed{false};
    bool         m_curve_can_extend{true};
    bool         m_curve_extended_constant{false};

    integer     m_npts{0};
    integer     m_npts_reserved{0};
    real_type * m_X{nullptr}; // allocated in the derived class!
    real_type * m_Y{nullptr}; // allocated in the derived class!

    SearchInterval m_search;

  protected:

    void
    copy_flags( Spline const & S ) {
      m_curve_is_closed         = S.m_curve_is_closed;
      m_curve_can_extend        = S.m_curve_can_extend;
      m_curve_extended_constant = S.m_curve_extended_constant;
    }

  public:

    Spline( Spline const & ) = delete;
    Spline const & operator = ( Spline const & ) = delete;

    //! \name Constructors
    ///@{

    //!
    //! spline constructor
    //!
    explicit
    Spline( string_view const name = "Spline" )
    : m_name(name)
    {
      m_search.setup( &m_name, &m_npts, &m_X, &m_curve_is_closed, &m_curve_can_extend );
    }

    //!
    //! spline destructor
    //!
    virtual ~Spline() = default;

    ///@}

    //!
    //! \name Open/Close
    //!
    ///@{

    //!
    //! \return string with the name of the spline
    //!
    string_view name() const { return m_name; }

    //! \return `true` if spline is a closed spline
    bool is_closed() const { return m_curve_is_closed;  }
    //!
    //! Set spline as a closed spline.
    //! When evaluated if parameter is outside the domain
    //! is wrapped cyclically before evalation.
    //!
    void make_closed() { m_curve_is_closed = true;  }
    //!
    //! Set spline as an opened spline.
    //! When evaluated if parameter is outside the domain
    //! an error is produced.
    //!
    void make_opened() { m_curve_is_closed = false; }

    //!
    //! \return `true` if spline cannot extend outside interval of definition
    //!
    bool is_bounded() const { return !m_curve_can_extend; }
    //!
    //! Set spline as unbounded.
    //! When evaluated if parameter is outside the domain
    //! an extrapolated value is used.
    //!
    void make_unbounded() { m_curve_can_extend = true;  }
    //!
    //! Set spline as bounded.
    //! When evaluated if parameter is outside the domain
    //! an error is issued.
    //!
    void make_bounded() { m_curve_can_extend = false; }

    //!
    //! \return `true` if the spline extend with a constant value
    //!
    bool is_extended_constant() const { return m_curve_extended_constant;  }
    //!
    //! Set spline to extend constant.
    //! When evaluated if parameter is outside the domain
    //! the value returned is the value of the closed border.
    //!
    void make_extended_constant()     { m_curve_extended_constant = true;  }
    //!
    //! Set spline to extend NOT constant.
    //! When evaluated if parameter is outside the domain
    //! teh value returned is extrapolated using the last spline polynomial.
    //!
    void make_extended_not_constant() { m_curve_extended_constant = false; }

    ///@}

    //! \name Spline Data Info
    ///@{

    //!
    //! return true if spline is empty (no points)
    //!
    bool empty() const { return m_npts == 0; }

    //!
    //! the number of support points of the spline.
    //!
    integer num_points() const { return m_npts; }

    //!
    //! Return the pointer of values of x-nodes.
    //!
    real_type const * x_nodes() const { return m_X; }

    //!
    //! Return the pointer of values of y-nodes.
    //!
    real_type const * y_nodes() const { return m_Y; }

    //!
    //! the i-th node of the spline (x component).
    //!
    real_type x_node( integer const i ) const { return m_X[i]; }

    //!
    //! the i-th node of the spline (y component).
    //!
    real_type y_node( integer const i ) const { return m_Y[i]; }

    //!
    //! first node of the spline (x component).
    //!
    real_type x_begin() const { return m_X[0]; }

    //!
    //! first node of the spline (y component).
    //!
    real_type y_begin() const { return m_Y[0]; }

    //!
    //! last node of the spline (x component).
    //!
    real_type x_end() const { return m_X[m_npts-1]; }

    //!
    //! last node of the spline (y component).
    //!
    real_type y_end() const { return m_Y[m_npts-1]; }

    //!
    //! x-minumum spline value
    //!
    real_type x_min() const { return m_X[0]; }

    //!
    //! x-maximum spline value
    //!
    real_type x_max() const { return m_X[m_npts-1]; }

    //!
    //! y-minumum spline value
    //!
    real_type
    y_min() const {
      integer N{m_npts};
      if ( type() == SplineType1D::CONSTANT ) --N;
      return *std::min_element(m_Y,m_Y+N);
    }

    //!
    //! return y-maximum spline value
    //!
    real_type
    y_max() const {
      integer N{m_npts};
      if ( type() == SplineType1D::CONSTANT ) --N;
      return *std::max_element(m_Y,m_Y+N);
    }

    //!
    //! Search the max and min values of `y` along the spline
    //! with the corresponding `x` position
    //!
    //! \param[out] i_min_pos interval where is the minimum
    //! \param[out] x_min_pos where is the minimum
    //! \param[out] y_min     the minimum value
    //! \param[out] i_max_pos interval where is the maximum
    //! \param[out] x_max_pos where is the maximum
    //! \param[out] y_max     the maximum value
    //!
    virtual
    void
    y_min_max(
      integer   & i_min_pos,
      real_type & x_min_pos,
      real_type & y_min,
      integer   & i_max_pos,
      real_type & x_max_pos,
      real_type & y_max
    ) const;

    //!
    //! Search the max and min values of `y` along the spline
    //! with the corresponding `x` position
    //!
    //! \param[out] i_min_pos interval where is the minimum
    //! \param[out] x_min_pos where is the minimum
    //! \param[out] y_min     the minimum value
    //! \param[out] i_max_pos interval where is the maximum
    //! \param[out] x_max_pos where is the maximum
    //! \param[out] y_max     the maximum value
    //!
    virtual
    void
    y_min_max(
      vector<integer>   & i_min_pos,
      vector<real_type> & x_min_pos,
      vector<real_type> & y_min,
      vector<integer>   & i_max_pos,
      vector<real_type> & x_max_pos,
      vector<real_type> & y_max
    ) const;

    ///@}

    //! \name Build
    ///@{

    //!
    //! Build a spline using data in `GenericContainer`
    //!
    void
    build( GenericContainer const & gc )
    { setup( gc ); }

    //!
    //! Build a spline using data in a file `file_name`
    //!
    void
    build( string const & file_name )
    { setup( file_name ); }

    //!
    //! Build a spline.
    //!
    //! \param x    vector of x-coordinates
    //! \param incx access elements as `x[0]`, `x[incx]`, `x[2*incx]`,...
    //! \param y    vector of y-coordinates
    //! \param incy access elements as `y[0]`, `y[incx]`, `y[2*incx]`,...
    //! \param n    total number of points
    //!
    virtual
    void
    build(
      real_type const x[], integer incx,
      real_type const y[], integer incy,
      integer   const n
    );

    //!
    //! Build a spline.
    //!
    //! \param x vector of x-coordinates
    //! \param y vector of y-coordinates
    //! \param n total number of points
    //!
    inline
    void
    build(
      real_type const x[],
      real_type const y[],
      integer   const n
    )
    { this->build( x, 1, y, 1, n ); }

    //!
    //! Build a spline.
    //!
    //! \param x vector of x-coordinates
    //! \param y vector of y-coordinates
    //!
    void
    build( vector<real_type> const & x, vector<real_type> const & y ) {
      integer N{ integer(x.size()) };
      if ( N > integer(y.size()) ) N = integer(y.size());
      this->build( x.data(), 1, y.data(), 1, N );
    }

    //!
    //! Build a spline using internal stored data
    //!
    virtual void build() = 0;

    //!
    //! Setup a spline using a `GenericContainer`
    //!
    //! - gc("xdata") vector with the `x` coordinate of the data
    //! - gc("ydata") vector with the `y` coordinate of the data
    //!
    virtual void setup( GenericContainer const & gc );

    //!
    //! Setup a spline using a `GenericContainer` readed from file
    //!
    //! - file_name file name of the file with data
    //!
    //! - `.json`
    //! - `.yaml` or `.yml`
    //! - `.toml`
    //!
    void setup( string const & file_name );

    ///@}

    //! \name Incremental Build
    ///@{

    //!
    //! Allocate memory for `npts` points
    //!
    virtual void reserve( integer const npts ) = 0;

    //!
    //! Add a support point (x,y) to the spline.
    //!
    void push_back( real_type const x, real_type const y );

    //!
    //! Drop last inserted point of the spline.
    //!
    void drop_back() { if ( m_npts > 0 ) --m_npts; m_search.must_reset(); }

    //!
    //! Delete the support points, empty the spline.
    //!
    virtual void clear() = 0;

    ///@}

    //! \name Manipulate
    ///@{

    //!
    //! change X-origin of the spline
    //!
    void set_origin( real_type const x0 ) const;

    //!
    //! change X-range of the spline
    //!
    void set_range( real_type const xmin, real_type const xmax );

    ///@}

    //! \name Dump Data
    ///@{

    //!
    //! dump a sample of the spline
    //!
    void
    dump(
      ostream_type & s,
      integer const  nintervals,
      string_view    header = "x\ty"
    ) const;

    //!
    //! dump a sample of the spline
    //!
    void
    dump(
      string_view   fname,
      integer const nintervals,
      string_view   header = "x\ty"
    ) const {
      std::ofstream file(fname.data());
      this->dump( file, nintervals, header );
      file.close();
    }

    //!
    //! Print spline coefficients
    //!
    virtual void write_to_stream( ostream_type & s ) const = 0;

    ///@}

    //! \name Evaluation
    ///@{

    #ifdef AUTODIFF_SUPPORT
    //!
    //! Evaluate spline value
    //!
    virtual real_type         eval( real_type         const   x ) const = 0;
    virtual autodiff::dual1st eval( autodiff::dual1st const & x ) const = 0;
    virtual autodiff::dual2nd eval( autodiff::dual2nd const & x ) const = 0;
    #endif

    //!
    //! First derivative
    //!
    virtual real_type D( real_type const x ) const = 0;

    //!
    //! Second derivative
    //!
    virtual real_type DD( real_type const x ) const = 0;

    //!
    //! Third derivative
    //!
    virtual real_type DDD( real_type const x ) const = 0;

    //!
    //! 4th derivative
    //!
    virtual real_type DDDD( real_type const ) const { return real_type(0); }

    //!
    //! 5th derivative
    //!
    virtual real_type DDDDD( real_type const ) const { return real_type(0); }

    virtual void D  ( real_type const x, real_type dd[2] ) const = 0;
    virtual void DD ( real_type const x, real_type dd[3] ) const = 0;

    ///@}

    //!
    //! \name Evaluation Aliases
    //!
    ///@{
    //!
    //! Alias for `real_type eval( real_type x )`
    //!
    real_type operator () ( real_type const x ) const { return this->eval(x); }
    //!
    //! Alias for `real_type D( real_type x )`
    //!
    real_type eval_D( real_type const x ) const { return this->D(x); }
    //!
    //! Alias for `real_type DD( real_type x )`
    //!
    real_type eval_DD( real_type const x ) const { return this->DD(x); }
    //!
    //! Alias for `real_type DDD( real_type x )`
    //!
    real_type eval_DDD( real_type const x ) const { return this->DDD(x); }
    //!
    //! Alias for `real_type DDDD( real_type x )`
    //!
    real_type eval_DDDD( real_type const x ) const { return this->DDDD(x); }
    //!
    //! Alias for `real_type DDDD( real_type x )`
    //!
    real_type eval_DDDDD( real_type const x ) const { return this->DDDDD(x); }
    ///@}

    //!
    //! \name Evaluation when segment is known
    //!
    ///@{

    //!
    //! Evaluate spline value
    //!
    virtual real_type id_eval( integer const ni, real_type const x ) const = 0;

    //!
    //! First derivative
    //!
    virtual real_type id_D( integer const ni, real_type const x ) const = 0;

    //!
    //! Second derivative
    //!
    virtual real_type id_DD( integer const ni, real_type const x ) const = 0;

    //!
    //! Third derivative
    //!
    virtual real_type id_DDD( integer const ni, real_type const x ) const = 0;

    //!
    //! 4th derivative
    //!
    virtual real_type id_DDDD( integer const, real_type const) const { return real_type(0); }

    //!
    //! 5th derivative
    //!
    virtual real_type id_DDDDD( integer const, real_type const) const { return real_type(0); }

    ///@}

    //! \name Get Info
    ///@{

    //!
    //! get the piecewise polinomials of the spline
    //!
    virtual
    integer // order
    coeffs(
      real_type cfs[],
      real_type nodes[],
      bool      transpose = false
    ) const = 0;

    //!
    //! Spline order of the piecewise polynomial
    //!
    virtual integer order() const = 0;

    //!
    //! spline type returned as a string
    //!
    char const *
    type_name() const
    { return to_string(type()); }

    //!
    //! spline type returned as integer
    //!
    virtual SplineType1D type() const = 0;

    //!
    //! String information of the kind and order of the spline
    //!
    string info() const;

    //!
    //! Print information of the kind and order of the spline
    //!
    void
    info( ostream_type & stream ) const
    { stream << this->info() << '\n'; }

    ///@}

    #ifdef SPLINES_BACK_COMPATIBILITY
    void pushBack( real_type x, real_type y ) { push_back(x,y); }
    void dropBack() { drop_back(); }
    void setOrigin( real_type x0 ) { set_origin(x0); }
    void setRange( real_type xmin, real_type xmax ) { set_range( xmin, xmax ); }
    void writeToStream( ostream_type & s ) const { write_to_stream(s); }
    #endif

  };

  //!
  //! compute curvature of a planar curve
  //!
  real_type curvature( real_type const s, Spline const & X, Spline const & Y );

  //!
  //! compute curvature derivative of a planar curve
  //!
  real_type curvature_D( real_type const s, Spline const & X, Spline const & Y );

  //!
  //! compute curvature second derivative of a planar curve
  //!
  real_type curvature_DD( real_type const s, Spline const & X, Spline const & Y );

  /*\
   |    ____      _     _        ____        _ _              ____
   |   / ___|   _| |__ (_) ___  / ___| _ __ | (_)_ __   ___  | __ )  __ _ ___  ___
   |  | |  | | | | '_ \| |/ __| \___ \| '_ \| | | '_ \ / _ \ |  _ \ / _` / __|/ _ \
   |  | |__| |_| | |_) | | (__   ___) | |_) | | | | | |  __/ | |_) | (_| \__ \  __/
   |   \____\__,_|_.__/|_|\___| |____/| .__/|_|_|_| |_|\___| |____/ \__,_|___/\___|
   |                                  |_|
  \*/
  //!
  //! cubic spline base class
  //!
  class CubicSplineBase : public Spline {
  protected:

    Malloc_real m_mem_cubic;
    real_type * m_Yp{nullptr};
    bool        m_external_alloc{false};

  public:

    using Spline::build;

    //!
    //! \name Contructors/Destructors
    ///@{
    //!
    //! Spline constructor.
    //!
    explicit
    CubicSplineBase( string_view name = "CubicSplineBase" );

    ~CubicSplineBase() override {}
    ///@}

    //!
    //! Build a copy of spline `S`
    //!
    void copy_spline( CubicSplineBase const & S );

    //!
    //! Return the pointer of values of yp-nodes.
    //!
    real_type const * yp_nodes() const { return m_Yp; }

    //!
    //! Return the i-th node of the spline (y' component).
    //!
    real_type yp_node( integer i ) const { return m_Yp[i]; }

    //!
    //! Change X-range of the spline.
    //!
    void set_range( real_type xmin, real_type xmax );

    //!
    //! Use externally allocated memory for `npts` points.
    //!
    void
    reserve_external(
      integer       n,
      real_type * & p_x,
      real_type * & p_y,
      real_type * & p_dy
    );

    void
    y_min_max(
      integer   & i_min_pos,
      real_type & x_min_pos,
      real_type & y_min,
      integer   & i_max_pos,
      real_type & x_max_pos,
      real_type & y_max
    ) const override;

    void
    y_min_max(
      vector<integer>   & i_min_pos,
      vector<real_type> & x_min_pos,
      vector<real_type> & y_min,
      vector<integer>   & i_max_pos,
      vector<real_type> & x_max_pos,
      vector<real_type> & y_max
    ) const override;

    // --------------------------- VIRTUALS -----------------------------------

    //!
    //! \name Evaluation
    //!
    ///@{
    real_type eval  ( real_type const x ) const override;
    real_type D     ( real_type const x ) const override;
    real_type DD    ( real_type const x ) const override;
    real_type DDD   ( real_type const x ) const override;
    real_type DDDD  ( real_type const   ) const override { return 0; }
    real_type DDDDD ( real_type const   ) const override { return 0; }

    void D  ( real_type const x, real_type dd[2] ) const override;
    void DD ( real_type const x, real_type dd[3] ) const override;

    ///@}

    //!
    //! \name Evaluation when segment is known
    ///@{
    real_type id_eval  ( integer const ni, real_type const x ) const override;
    real_type id_D     ( integer const ni, real_type const x ) const override;
    real_type id_DD    ( integer const ni, real_type const x ) const override;
    real_type id_DDD   ( integer const ni, real_type const x ) const override;
    real_type id_DDDD  ( integer const   , real_type const   ) const override { return 0; }
    real_type id_DDDDD ( integer const   , real_type const   ) const override { return 0; }
    ///@}

    #ifdef AUTODIFF_SUPPORT
    autodiff::dual1st eval( autodiff::dual1st const & x ) const override;
    autodiff::dual2nd eval( autodiff::dual2nd const & x ) const override;

    template <typename T>
    autodiff::HigherOrderDual<autodiff::detail::DualOrder<T>::value,real_type>
    eval( T const & x ) const { return eval( autodiff::detail::to_dual(x) ); }

    template <typename T>
    autodiff::HigherOrderDual<autodiff::detail::DualOrder<T>::value,real_type>
    operator () ( T const & x ) const { return eval( autodiff::detail::to_dual(x) ); }
    #endif

    void write_to_stream( ostream_type & s ) const override;

    // --------------------------- VIRTUALS -----------------------------------

    //!
    //! \name Build
    //!
    ///@{

    //!
    //! Build a spline.
    //!
    //! \param[in] x     vector of x-coordinates
    //! \param[in] incx  access elements as x[0], x[incx], x[2*incx],...
    //! \param[in] y     vector of y-coordinates
    //! \param[in] incy  access elements as y[0], y[incy], y[2*incy],...
    //! \param[in] yp    vector of y-defivative
    //! \param[in] incyp access elements as yp[0], yp[incyp], yp[2*incyy],...
    //! \param[in] n     total number of points
    //!
    void
    build(
      real_type const x[],  integer const incx,
      real_type const y[],  integer const incy,
      real_type const yp[], integer const incyp,
      integer n
    );

    //!
    //! Build a spline.
    //!
    //! \param[in] x  vector of x-coordinates
    //! \param[in] y  vector of y-coordinates
    //! \param[in] yp vector of y'-coordinates
    //! \param[in] n  total number of points
    //!
    inline
    void
    build(
      real_type const x[],
      real_type const y[],
      real_type const yp[],
      integer   const n
    ) {
      this->build( x, 1, y, 1, yp, 1, n );
    }

    //!
    //! Build a spline.
    //!
    //! \param[in] x  vector of x-coordinates
    //! \param[in] y  vector of y-coordinates
    //! \param[in] yp vector of y'-coordinates
    //!
    void
    build(
      vector<real_type> const & x,
      vector<real_type> const & y,
      vector<real_type> const & yp
    );

    void reserve( integer npts ) override;

    ///@}

    void clear() override;

    integer // order
    coeffs(
      real_type cfs[],
      real_type nodes[],
      bool      transpose = false
    ) const override;

    integer order() const override;

    #ifdef SPLINES_BACK_COMPATIBILITY
    void copySpline( CubicSplineBase const & S ) { this->copy_spline(S); }
    integer numPoints() const { return m_npts; }
    real_type xNode( integer i ) const { return m_X[i]; }
    real_type yNode( integer i ) const { return m_Y[i]; }
    real_type ypNode( integer i ) const { return this->yp_node(i); }
    real_type xBegin() const { return m_X[0]; }
    real_type yBegin() const { return m_Y[0]; }
    real_type xEnd() const { return m_X[m_npts-1]; }
    real_type yEnd() const { return m_Y[m_npts-1]; }
    real_type xMin() const { return m_X[0]; }
    real_type xMax() const { return m_X[m_npts-1]; }
    real_type yMin() const { return y_min(); }
    real_type yMax() const { return y_max(); }
    #endif

  };

  /*\
   |   ____        _ _            ____              __
   |  / ___| _ __ | (_)_ __   ___/ ___| _   _ _ __ / _|
   |  \___ \| '_ \| | | '_ \ / _ \___ \| | | | '__| |_
   |   ___) | |_) | | | | | |  __/___) | |_| | |  |  _|
   |  |____/| .__/|_|_|_| |_|\___|____/ \__,_|_|  |_|
   |        |_|
  \*/

  //!
  //! Spline Management Class
  //!
  class SplineSurf {

    Malloc_real m_mem;

  protected:

    string const m_name;
    bool         m_x_closed{false};
    bool         m_y_closed{false};
    bool         m_x_can_extend{true};
    bool         m_y_can_extend{true};

    integer      m_nx{0};
    integer      m_ny{0};

    real_type *  m_X{nullptr};
    real_type *  m_Y{nullptr};
    real_type *  m_Z{nullptr};

    real_type    m_Z_min{0};
    real_type    m_Z_max{0};

    SearchInterval m_search_x;
    SearchInterval m_search_y;

    static
    integer
    ipos_C( integer const i, integer const j, integer const ldZ )
    { return i*ldZ + j; }

    static
    integer
    ipos_F( integer const i, integer const j, integer const ldZ )
    { return i + ldZ*j; }

    integer
    ipos_C( integer const i, integer const j ) const
    { return this->ipos_C(i,j,m_ny); }

    integer
    ipos_F( integer const i, integer const j ) const
    { return this->ipos_F(i,j,m_nx); }

    real_type &
    z_node_ref( integer const i, integer const j )
    { return m_Z[this->ipos_C(i,j)]; }

    void
    load_Z(
      real_type const z[],
      integer   const ldZ,
      bool            fortran_storage,
      bool            transposed
    );

    virtual void make_spline() = 0;

    void make_derivative_x  ( real_type const z[], real_type dz[] );
    void make_derivative_y  ( real_type const z[], real_type dz[] );
    void make_derivative_xy ( real_type const dx[], real_type const dy[], real_type dxy[] );

  public:

    SplineSurf( SplineSurf const & ) = delete; // block copy constructor
    SplineSurf const & operator = ( SplineSurf const & ) = delete; // block copy method

    //!
    //! Spline constructor
    //!
    explicit
    SplineSurf( string_view name = "Spline" )
    : m_mem(name.data())
    , m_name(name)
    {
      m_search_x.setup( &m_name, &m_nx, &m_X, &m_x_closed, &m_x_can_extend );
      m_search_y.setup( &m_name, &m_ny, &m_Y, &m_y_closed, &m_y_can_extend );
    }

    //!
    //! Spline destructor
    //!
    virtual
    ~SplineSurf();

    //!
    //! \name Open/Close
    //!
    ///@{

    //!
    //! Return `true` if the surface is assumed closed in the `x` direction.
    //!
    bool is_x_closed() const { return m_x_closed; }

    //!
    //! Setup the surface as closed in the `x` direction.
    //!
    void make_x_closed() { m_x_closed = true; }

    //!
    //! Setup the surface as open in the `x` direction.
    //!
    void make_x_opened() { m_x_closed = false; }

    //!
    //! Return `true` if the surface is assumed closed in the `y` direction.
    //!
    bool is_y_closed() const { return m_y_closed; }

    //!
    //! Setup the surface as closed in the `y` direction.
    //!
    void make_y_closed() { m_y_closed = true; }

    //!
    //! Setup the surface as open in the `y` direction.
    //!
    void make_y_opened() { m_y_closed = false; }

    //!
    //! Return `true` if the parameter `x` assumed bounded.
    //! If false the spline is estrapolated for `x` values
    //! outside the range.
    //!
    bool is_x_bounded() const { return m_x_can_extend; }

    //!
    //! Make the spline surface unbounded in the `x` direction.
    //!
    void make_x_unbounded() { m_x_can_extend = true; }

    //!
    //! Make the spline surface bounded in the `x` direction.
    //!
    void make_x_bounded() { m_x_can_extend = false; }

    //!
    //! Return `true` if the parameter `y` assumed bounded.
    //! If false the spline is extrapolated for `y` values
    //! outside the range.
    //!
    bool is_y_bounded() const { return m_y_can_extend; }

    //!
    //! Make the spline surface unbounded in the `y` direction
    //!
    void make_y_unbounded() { m_y_can_extend = true; }

    //!
    //! Make the spline surface bounded in the `x` direction.
    //!
    void make_y_bounded() { m_y_can_extend = false; }

    ///@}

    //!
    //! Cancel the support points, empty the spline.
    //!
    void clear();

    //!
    //! \name Info
    //!
    ///@{

    //!
    //! \return string with the name of the spline
    //!
    string_view name() const { return m_name; }

    //!
    //! Return the number of support points of the spline along x direction.
    //!
    integer num_point_x() const { return m_nx; }

    //!
    //! Return the number of support points of the spline along y direction.
    //!
    integer num_point_y() const { return m_ny; }

    //!
    //! Return the i-th node of the spline (x component).
    //!
    real_type x_node( integer const i ) const { return m_X[i]; }

    //!
    //! Return the i-th node of the spline (y component).
    //!
    real_type y_node( integer const i ) const { return m_Y[i]; }

    //!
    //! Return the i-th node of the spline (y component).
    //!
    real_type z_node( integer const i, integer const j ) const { return m_Z[this->ipos_C(i,j)]; }

    //!
    //! Return x-minumum spline value.
    //!
    real_type x_min() const { return m_X[0]; }

    //!
    //! Return x-maximum spline value.
    //!
    real_type x_max() const { return m_X[m_nx-1]; }

    //!
    //! Return y-minumum spline value.
    //!
    real_type y_min() const { return m_Y[0]; }

    //!
    //! Return y-maximum spline value.
    //!
    real_type y_max() const { return m_Y[m_ny-1]; }

    //!
    //! Return z-minumum spline value.
    //!
    real_type z_min() const { return m_Z_min; }

    //!
    //! Return z-maximum spline value.
    //!
    real_type z_max() const { return m_Z_max; }

    ///@}

    //!
    //! \name Build Spline
    //!
    ///@{

    //!
    //! Build surface spline
    //!
    //! \param x               vector of `x`-coordinates
    //! \param incx            access elements as `x[0]`, `x[incx]`, `x[2*incx]`,...
    //! \param y               vector of `y`-coordinates
    //! \param incy            access elements as `y[0]`, `y[incy]`, `y[2*incy]`,...
    //! \param z               matrix of `z`-values. Elements are stored
    //!                        by row Z(i,j) = z[i*ny+j] as C-matrix
    //! \param ldZ             leading dimension of `z`
    //! \param nx              number of points in `x` direction
    //! \param ny              number of points in `y` direction
    //! \param fortran_storage if true elements are stored by column
    //!                        i.e. Z(i,j) = z[i+j*nx] as Fortran-matrix
    //! \param transposed      if true matrix Z is stored transposed
    //!
    void
    build(
      real_type const x[], integer const incx,
      real_type const y[], integer const incy,
      real_type const z[], integer const ldZ,
      integer   const nx,  integer const ny,
      bool fortran_storage = false,
      bool transposed      = false
    );

    //!
    //! Build surface spline
    //!
    //! \param x               vector of x-coordinates, nx = x.size()
    //! \param y               vector of y-coordinates, ny = y.size()
    //! \param z               matrix of z-values. Elements are stored
    //!                        by row Z(i,j) = z[i*ny+j] as C-matrix
    //! \param fortran_storage if true elements are stored by column
    //!                        i.e. Z(i,j) = z[i+j*nx] as Fortran-matrix
    //! \param transposed      if true matrix Z is stored transposed
    //!
    void
    build(
      vector<real_type> const & x,
      vector<real_type> const & y,
      vector<real_type> const & z,
      bool fortran_storage = false,
      bool transposed      = false
    ) {
      this->build(
        x.data(), 1,
        y.data(), 1,
        z.data(), integer(fortran_storage ? y.size() : x.size()),
        integer(x.size()), integer(y.size()),
        fortran_storage, transposed
      );
    }

    void
    build(
      real_type const z[],
      integer         ldZ,
      integer         nx,
      integer         ny,
      bool fortran_storage = false,
      bool transposed      = false
    );

    //!
    //! Build surface spline
    //!
    //! \param z               matrix of z-values. Elements are stored
    //!                        by row Z(i,j) = z[i*ny+j] as C-matrix.
    //!                        ldZ leading dimension of the matrix is ny for C-storage
    //!                        and nx for Fortran storage.
    //! \param nx              x-dimension
    //! \param ny              y-dimension
    //! \param fortran_storage if true elements are stored by column
    //!                        i.e. Z(i,j) = z[i+j*nx] as Fortran-matrix
    //! \param transposed      if true matrix Z is stored transposed
    //!
    void
    build(
      vector<real_type> const & z,
      integer           const   nx,
      integer           const   ny,
      bool fortran_storage = false,
      bool transposed      = false
    ) {
      this->build( z.data(), nx, ny, fortran_storage ? nx : ny, fortran_storage, transposed );
    }

    //!
    //! Build spline using data in `gc`
    //!
    void
    setup( GenericContainer const & gc );

    //!
    //! Setup a spline using a `GenericContainer` readed from file
    //!
    //! - file_name file name of the file with data
    //!
    //! - `.json`
    //! - `.yaml` or `.yml`
    //! - `.toml`
    //!
    void setup( string const & file_name );

    //!
    //! Build a spline using data in `GenericContainer`
    //!
    void
    build( GenericContainer const & gc )
    { setup( gc ); }

    //!
    //! Build a spline using data in a file `file_name`
    //!
    void
    build( string const & file_name )
    { setup( file_name ); }

    ///@}

    //!
    //! \name Evaluate
    //!
    ///@{

    //!
    //! Evaluate spline value at point \f$ (x,y) \f$.
    //!
    virtual real_type eval( real_type const x, real_type const y ) const = 0;

    #ifdef AUTODIFF_SUPPORT
    autodiff::dual1st eval( autodiff::dual1st const & x, autodiff::dual1st const & y ) const;
    autodiff::dual2nd eval( autodiff::dual2nd const & x, autodiff::dual2nd const & y ) const;

    template <typename T1, typename T2>
    autodiff::HigherOrderDual<autodiff::detail::DualOrder<T1,T2>::value,real_type>
    eval( T1 const & x, T2 const & y ) const {
      autodiff::HigherOrderDual<autodiff::detail::DualOrder<T1,T2>::value,real_type> X{x}, Y{y};
      return eval( X, Y );
    }
    #endif

    //!
    //! Value and first derivatives at point \f$ (x,y) \f$:
    //!
    //! - d[0] value of the spline \f$ S(x,y) \f$
    //! - d[1] derivative respect to \f$ x \f$ of the spline: \f$ S_x(x,y) \f$
    //! - d[2] derivative respect to \f$ y \f$ of the spline: \f$ S_y(x,y) \f$
    //!
    virtual void D( real_type const x, real_type const y, real_type d[3] ) const = 0;

    //!
    //! First derivatives respect to \f$ x \f$ at point \f$ (x,y) \f$
    //! of the spline: \f$ S_x(x,y) \f$.
    //!
    virtual real_type Dx( real_type const x, real_type const y ) const = 0;

    //!
    //! First derivatives respect to \f$ y \f$ at point \f$ (x,y) \f$
    //! of the spline: \f$ S_y(x,y) \f$.
    //!
    virtual real_type Dy( real_type const x, real_type const y ) const = 0;

    //!
    //! Value, first and second derivatives at point \f$ (x,y) \f$:
    //!
    //! - dd[0] value of the spline \f$ S(x,y) \f$
    //! - dd[1] derivative respect to \f$ x \f$ of the spline: \f$ S_x(x,y) \f$
    //! - dd[2] derivative respect to \f$ y \f$ of the spline: \f$ S_y(x,y) \f$
    //! - dd[3] second derivative respect to \f$ x \f$ of the spline: \f$ S_{xx}(x,y) \f$
    //! - dd[4] mixed second derivative: \f$ S_{xy}(x,y) \f$
    //! - dd[5] second derivative respect to \f$ y \f$ of the spline: \f$ S_{yy}(x,y) \f$
    //!
    virtual void DD( real_type const x, real_type const y, real_type dd[6] ) const = 0;

    //!
    //! Second derivatives respect to \f$ x \f$ at point \f$ (x,y) \f$
    //! of the spline: \f$ S_{xx}(x,y) \f$.
    //!
    virtual real_type Dxx( real_type const x, real_type const y ) const = 0;

    //!
    //! Mixed second derivatives: \f$ S_{xy}(x,y) \f$.
    //!
    virtual real_type Dxy( real_type const x, real_type const y ) const = 0;

    //!
    //! Second derivatives respect to \f$ y \f$ at point \f$ (x,y) \f$
    //! of the spline: \f$ S_{yy}(x,y) \f$.
    //!
    virtual real_type Dyy( real_type const x, real_type const y ) const = 0;

    //!
    //! Evaluate spline value at point \f$ (x,y) \f$.
    //!
    real_type operator () ( real_type const x, real_type const y ) const { return this->eval(x,y); }

    //!
    //! Alias for `Dx(x,y)`
    //!
    real_type eval_D_1( real_type const x, real_type const y ) const { return this->Dx(x,y); }

    //!
    //! Alias for `Dy(x,y)`
    //!
    real_type eval_D_2( real_type const x, real_type const y ) const { return this->Dy(x,y); }

    //!
    //! Alias for `Dxx(x,y)`
    //!
    real_type eval_D_1_1( real_type const x, real_type const y ) const { return this->Dxx(x,y); }

    //!
    //! Alias for `Dxy(x,y)`
    //!
    real_type eval_D_1_2( real_type const x, real_type const y ) const { return this->Dxy(x,y); }

    //!
    //! Alias for `Dyy(x,y)`
    //!
    real_type eval_D_2_2( real_type const x, real_type const y ) const { return this->Dyy(x,y); }

    ///@}

    //!
    //! Print spline coefficients.
    //!
    virtual void write_to_stream( ostream_type & s ) const = 0;

    //!
    //! Return spline type as a string pointer.
    //!
    virtual char const * type_name() const = 0;

    //!
    //! String information of the kind and order of the spline
    //!
    virtual string info() const;

    //!
    //! Print information of the kind and order of the spline
    //!
    void
    info( ostream_type & stream ) const
    { stream << this->info() << '\n'; }

    //!
    //! Print stored data x, y, and matrix z.
    //!
    void dump_data( ostream_type & s ) const;

    #ifdef SPLINES_BACK_COMPATIBILITY
    integer numPointX() const { return m_nx; }
    integer numPointY() const { return m_ny; }
    real_type xNode( integer i ) const { return m_X[i]; }
    real_type yNode( integer i ) const { return m_Y[i]; }
    real_type zNode( integer i, integer j ) const { return z_node(i,j); }
    real_type xMin() const { return this->x_min(); }
    real_type xMax() const { return this->x_max(); }
    real_type yMin() const { return this->y_min(); }
    real_type yMax() const { return this->y_max(); }
    real_type zMin() const { return m_Z_min; }
    real_type zMax() const { return m_Z_max; }
    void writeToStream( ostream_type & s ) const { write_to_stream(s); }
    #endif

  };

}

#include "Splines/SplineAkima.hxx"
#include "Splines/SplineBessel.hxx"
#include "Splines/SplineConstant.hxx"
#include "Splines/SplineLinear.hxx"
#include "Splines/SplineCubic.hxx"
#include "Splines/SplineHermite.hxx"
#include "Splines/SplinePchip.hxx"
#include "Splines/SplineQuinticBase.hxx"
#include "Splines/SplineQuintic.hxx"
#include "Splines/SplineBilinear.hxx"
#include "Splines/SplineBiCubic.hxx"
#include "Splines/SplineAkima2D.hxx"
#include "Splines/SplineBiQuintic.hxx"

#include "Splines/SplineVec.hxx"
#include "Splines/SplineSet.hxx"
#include "Splines/Splines1D.hxx"
#include "Splines/Splines2D.hxx"

#include "Splines/Splines1Dblend.hxx"
#include "Splines/Splines2Dblend.hxx"

namespace SplinesLoad {

  using Splines::Spline;
  using Splines::CubicSplineBase;
  using Splines::CubicSpline;
  using Splines::AkimaSpline;
  using Splines::BesselSpline;
  using Splines::PchipSpline;
  using Splines::LinearSpline;
  using Splines::ConstantSpline;
  using Splines::QuinticSpline;
  using Splines::Spline1D;

  using Splines::BilinearSpline;
  using Splines::BiCubicSpline;
  using Splines::BiQuinticSpline;
  using Splines::Akima2Dspline;
  using Splines::Spline2D;

  using Splines::SplineVec;
  using Splines::SplineSet;

  using Splines::SplineType1D;
  using Splines::SplineType2D;
}

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
#ifdef __clang__
#pragma clang diagnostic pop
#endif

#endif
