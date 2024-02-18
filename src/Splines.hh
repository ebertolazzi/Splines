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

namespace Splines {

  using std::vector;
  using std::string;
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

  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  extern SplineType1D string_to_splineType1D( string const & n );
  extern SplineType2D string_to_splineType2D( string const & n );
  extern char const * to_string( SplineType2D t );
  extern char const * to_string( SplineType1D t );
  #endif

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

  #ifndef DOXYGEN_SHOULD_SKIP_THIS

  void Hermite3( real_type x, real_type H, real_type base[4] );
  void Hermite3_D( real_type x, real_type H, real_type base_D[4] );
  void Hermite3_DD( real_type x, real_type H, real_type base_DD[4] );
  void Hermite3_DDD( real_type x, real_type H, real_type base_DDD[4] );

  void Hermite5( real_type x, real_type H, real_type base[6] );
  void Hermite5_D( real_type x, real_type H, real_type base_D[6] );
  void Hermite5_DD( real_type x, real_type H, real_type base_DD[6] );
  void Hermite5_DDD( real_type x, real_type H, real_type base_DDD[6] );
  void Hermite5_DDDD( real_type x, real_type H, real_type base_DDDD[6] );
  void Hermite5_DDDDD( real_type x, real_type H, real_type base_DDDDD[6] );

  #endif

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
    real_type   H,
    real_type   P0,
    real_type   P1,
    real_type   DP0,
    real_type   DP1,
    real_type & A,
    real_type & B,
    real_type & C,
    real_type & D
  ) {
    real_type H2  = H*H;
    real_type P10 = P1-P0;
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
    real_type   h,
    real_type   P0,
    real_type   P1,
    real_type   DP0,
    real_type   DP1,
    real_type   DDP0,
    real_type   DDP1,
    real_type & A,
    real_type & B,
    real_type & C,
    real_type & D,
    real_type & E,
    real_type & F
  ) {
    real_type h2  = h*h;
    real_type h3  = h*h2;
    real_type P10 = P1-P0;
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

  #ifndef DOXYGEN_SHOULD_SKIP_THIS

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
  checkCubicSplineMonotonicity(
    real_type const * X,
    real_type const * Y,
    real_type const * Yp,
    integer           npts
  );

  #endif

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
    integer           dim,
    integer           npts,
    real_type const * pnts,
    integer           ld_pnts,
    real_type       * t
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
    integer           dim,
    integer           npts,
    real_type const * pnts,
    integer           ld_pnts,
    real_type       * t
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
    integer           dim,
    integer           npts,
    real_type const * pnts,
    integer           ld_pnts,
    real_type         alpha,
    real_type       * t
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
    integer           dim,
    integer           npts,
    real_type const * pnts,
    integer           ld_pnts,
    real_type       * t
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
    integer           dim,
    integer           npts,
    real_type const * pnts,
    integer           ld_pnts,
    real_type       * t
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
    integer           dim,
    integer           npts,
    real_type const * pnts,
    integer           ld_pnts,
    real_type       * t
  );

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

    string m_name;
    bool   m_curve_is_closed{false};
    bool   m_curve_can_extend{true};
    bool   m_curve_extended_constant{false};

    integer     m_npts{0};
    integer     m_npts_reserved{0};
    real_type * m_X{nullptr}; // allocated in the derived class!
    real_type * m_Y{nullptr}; // allocated in the derived class!

    #ifdef SPLINES_USE_THREADS
    mutable Utils::BinarySearch<integer> m_last_interval;
    #else
    mutable integer m_last_interval;
    #endif

    void init_last_interval();

    Spline( Spline const & ) = delete;
    Spline const & operator = ( Spline const & ) = delete;

  public:

    //! \name Constructors
    ///@{

    //!
    //! spline constructor
    //!
    Spline( string const & name = "Spline" )
    : m_name(name)
    {
      this->init_last_interval();
    }

    //!
    //! spline destructor
    //!
    virtual
    ~Spline()
    {}

    ///@}

    //!
    //! Find interval containing `x` using binary search.
    //!
    integer search( real_type & x ) const;

    //!
    //! \name Open/Close
    //!
    ///@{

    //!
    //! Return the name of the spline used in the constructor
    //!
    string const & name() const { return m_name; }

    bool is_closed() const { return m_curve_is_closed;  }
    void make_closed()     { m_curve_is_closed = true;  }
    void make_opened()     { m_curve_is_closed = false; }

    bool is_bounded() const { return !m_curve_can_extend; }
    void make_unbounded()   { m_curve_can_extend = true;  }
    void make_bounded()     { m_curve_can_extend = false; }

    bool is_extended_constant() const { return m_curve_extended_constant;  }
    void make_extended_constant()     { m_curve_extended_constant = true;  }
    void make_extended_not_constant() { m_curve_extended_constant = false; }

    ///@}

    //! \name Spline Data Info
    ///@{

    //!
    //! the number of support points of the spline.
    //!
    integer num_points() const { return m_npts; }

    //!
    //! the i-th node of the spline (x component).
    //!
    real_type x_node( integer i ) const { return m_X[size_t(i)]; }

    //!
    //! the i-th node of the spline (y component).
    //!
    real_type y_node( integer i ) const { return m_Y[size_t(i)]; }

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
    real_type x_end() const { return m_X[size_t(m_npts-1)]; }

    //!
    //! last node of the spline (y component).
    //!
    real_type y_end() const { return m_Y[size_t(m_npts-1)]; }

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
      integer N = m_npts;
      if ( type() == SplineType1D::CONSTANT ) --N;
      return *std::min_element(m_Y,m_Y+N);
    }

    //!
    //! return y-maximum spline value
    //!
    real_type
    y_max() const {
      integer N = m_npts;
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
    ) const {
      i_min_pos = i_max_pos = 0;
      x_min_pos = y_min = x_max_pos = y_max = 0;
      UTILS_ERROR(
        "In spline: {} y_min_max not implemented\n",
        info()
      );
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

    ///@}

    //! \name Build
    ///@{

    //!
    //! Build a spline using a generic container.
    //!
    void
    build( GenericContainer const & gc )
    { setup(gc); }

    //!
    //! Build a spline.
    //!
    //! \param x    vector of x-coordinates
    //! \param incx access elements as x[0], x[incx], x[2*incx],...
    //! \param y    vector of y-coordinates
    //! \param incy access elements as y[0], y[incy], y[2*incy],...
    //! \param n    total number of points
    //!
    virtual
    void
    build(
      real_type const * x, integer incx,
      real_type const * y, integer incy,
      integer n
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
      real_type const * x,
      real_type const * y,
      integer           n
    )
    { this->build( x, 1, y, 1, n ); }

    //!
    //! Build a spline.
    //!
    //! \param x vector of x-coordinates
    //! \param y vector of y-coordinates
    //!
    inline
    void
    build( vector<real_type> const & x, vector<real_type> const & y ) {
      integer N = integer(x.size());
      if ( N > integer(y.size()) ) N = integer(y.size());
      this->build( &x.front(), 1, &y.front(), 1, N );
    }

    //!
    //! Build a spline using internal stored data
    //!
    virtual
    void build() = 0;

    //!
    //! Setup a spline using a `GenericContainer`
    //!
    //! - gc("xdata") vector with the `x` coordinate of the data
    //! - gc("ydata") vector with the `y` coordinate of the data
    //!
    virtual
    void setup( GenericContainer const & gc );

    ///@}

    //! \name Incremental Build
    ///@{

    //!
    //! Allocate memory for `npts` points
    //!
    virtual void reserve( integer npts ) = 0;

    //!
    //! Add a support point (x,y) to the spline.
    //!
    void push_back( real_type x, real_type y );

    //!
    //! Drop last inserted point of the spline.
    //!
    void drop_back() { if ( m_npts > 0 ) --m_npts; }

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
    void set_origin( real_type x0 );

    //!
    //! change X-range of the spline
    //!
    void set_range( real_type xmin, real_type xmax );

    ///@}

    //! \name Dump Data
    ///@{

    //!
    //! dump a sample of the spline
    //!
    void
    dump(
      ostream_type & s,
      integer        nintervals,
      char const *   header = "x\ty"
    ) const;

    //!
    //! dump a sample of the spline
    //!
    void
    dump(
      char const * fname,
      integer      nintervals,
      char const * header = "x\ty"
    ) const {
      std::ofstream file(fname);
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

    //!
    //! Evaluate spline value
    //!
    virtual real_type operator () ( real_type x ) const = 0;

    //!
    //! First derivative
    //!
    virtual real_type D( real_type x ) const = 0;

    //!
    //! Second derivative
    //!
    virtual real_type DD( real_type x ) const = 0;

    //!
    //! Third derivative
    //!
    virtual real_type DDD( real_type x ) const = 0;

    //!
    //! 4th derivative
    //!
    virtual real_type DDDD( real_type ) const { return real_type(0); }

    //!
    //! 5th derivative
    //!
    virtual real_type DDDDD( real_type ) const { return real_type(0); }

    ///@}

    //!
    //! \name Evaluation Aliases
    //!
    ///@{
    //!
    //! Alias for `operator () ( real_type x )`
    //!
    real_type eval( real_type x ) const { return (*this)(x); }
    //!
    //! Alias for `real_type D( real_type x )`
    //!
    real_type eval_D( real_type x ) const { return this->D(x); }
    //!
    //! Alias for `real_type DD( real_type x )`
    //!
    real_type eval_DD( real_type x ) const { return this->DD(x); }
    //!
    //! Alias for `real_type DDD( real_type x )`
    //!
    real_type eval_DDD( real_type x ) const { return this->DDD(x); }
    //!
    //! Alias for `real_type DDDD( real_type x )`
    //!
    real_type eval_DDDD( real_type x ) const { return this->DDDD(x); }
    //!
    //! Alias for `real_type DDDD( real_type x )`
    //!
    real_type eval_DDDDD( real_type x ) const { return this->DDDDD(x); }
    ///@}

    //!
    //! \name Evaluation when segment is known
    //!
    ///@{

    //!
    //! Evaluate spline value
    //!
    virtual real_type id_eval( integer ni, real_type x ) const = 0;

    //!
    //! First derivative
    //!
    virtual real_type id_D( integer ni, real_type x ) const = 0;

    //!
    //! Second derivative
    //!
    virtual real_type id_DD( integer ni, real_type x ) const = 0;

    //!
    //! Third derivative
    //!
    virtual real_type id_DDD( integer ni, real_type x ) const = 0;

    //!
    //! 4th derivative
    //!
    virtual real_type id_DDDD( integer, real_type ) const { return real_type(0); }

    //!
    //! 5th derivative
    //!
    virtual real_type id_DDDDD( integer, real_type ) const { return real_type(0); }

    ///@}

    //! \name Get Info
    ///@{

    //!
    //! get the piecewise polinomials of the spline
    //!
    virtual
    integer // order
    coeffs(
      real_type * const cfs,
      real_type * const nodes,
      bool              transpose = false
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
    void info( ostream_type & stream ) const { stream << this->info() << '\n'; }

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
  real_type curvature( real_type s, Spline const & X, Spline const & Y );

  //!
  //! compute curvature derivative of a planar curve
  //!
  real_type curvature_D( real_type s, Spline const & X, Spline const & Y );

  //!
  //! compute curvature second derivative of a planar curve
  //!
  real_type curvature_DD( real_type s, Spline const & X, Spline const & Y );

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
    Malloc_real m_baseValue;
    real_type * m_Yp{nullptr};
    bool        m_external_alloc{false};

  public:

    #ifndef DOXYGEN_SHOULD_SKIP_THIS
    using Spline::build;
    #endif

    //!
    //! \name Contructors/Destructors
    ///@{
    //!
    //! Spline constructor.
    //!
    CubicSplineBase( string const & name = "CubicSplineBase" )
    : Spline(name)
    , m_baseValue(name+"_memory")
    {}

    ~CubicSplineBase() override {}
    ///@}

    void copy_spline( CubicSplineBase const & S );

    //!
    //! Return the i-th node of the spline (y' component).
    //!
    real_type yp_node( integer i ) const { return m_Yp[size_t(i)]; }

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
    real_type operator () ( real_type x ) const override;
    real_type D( real_type x ) const override;
    real_type DD( real_type x ) const override;
    real_type DDD( real_type x ) const override;
    real_type DDDD( real_type ) const override { return 0; }
    real_type DDDDD( real_type ) const override { return 0; }
    ///@}

    //!
    //! \name Evaluation when segment is known
    ///@{
    real_type id_eval( integer ni, real_type x ) const override;
    real_type id_D( integer ni, real_type x ) const override;
    real_type id_DD( integer ni, real_type x ) const override;
    real_type id_DDD( integer ni, real_type x ) const override;
    real_type id_DDDD( integer, real_type ) const override { return 0; }
    real_type id_DDDDD( integer, real_type ) const override { return 0; }
    ///@}

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
      real_type const * x,  integer incx,
      real_type const * y,  integer incy,
      real_type const * yp, integer incyp,
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
      real_type const * x,
      real_type const * y,
      real_type const * yp,
      integer n
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
      real_type * const cfs,
      real_type * const nodes,
      bool              transpose = false
    ) const override;

    integer order() const override;

    #ifdef SPLINES_BACK_COMPATIBILITY
    void copySpline( CubicSplineBase const & S ) { this->copy_spline(S); }
    integer numPoints() const { return m_npts; }
    real_type xNode( integer i ) const { return m_X[size_t(i)]; }
    real_type yNode( integer i ) const { return m_Y[size_t(i)]; }
    real_type ypNode( integer i ) const { return this->yp_node(i); }
    real_type xBegin() const { return m_X[0]; }
    real_type yBegin() const { return m_Y[0]; }
    real_type xEnd() const { return m_X[size_t(m_npts-1)]; }
    real_type yEnd() const { return m_Y[size_t(m_npts-1)]; }
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

    SplineSurf( SplineSurf const & ) = delete; // block copy constructor
    SplineSurf const & operator = ( SplineSurf const & ) = delete; // block copy method

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

    #ifdef SPLINES_USE_THREADS
    mutable Utils::BinarySearch<integer> m_last_interval_x;
    mutable Utils::BinarySearch<integer> m_last_interval_y;
    #else
    mutable integer m_last_interval_x;
    mutable integer m_last_interval_y;
    #endif

    integer search_x( real_type & x ) const;
    integer search_y( real_type & y ) const;

    void init_last_interval_x();
    void init_last_interval_y();

    integer
    ipos_C( integer i, integer j, integer ldZ ) const
    { return i*ldZ + j; }

    integer
    ipos_F( integer i, integer j, integer ldZ ) const
    { return i + ldZ*j; }

    integer
    ipos_C( integer i, integer j ) const
    { return this->ipos_C(i,j,m_ny); }

    integer
    ipos_F( integer i, integer j ) const
    { return this->ipos_F(i,j,m_nx); }

    virtual void make_spline() = 0;

    void
    load_Z(
      real_type const * z,
      integer           ldZ,
      bool              fortran_storage,
      bool              transposed
    );

  public:

    //!
    //! Spline constructor
    //!
    SplineSurf( string const & name = "Spline" )
    : m_mem("SplineSurf")
    , m_name(name)
    {
      this->init_last_interval_x();
      this->init_last_interval_y();
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

    string const & name() const { return m_name; }

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
    real_type x_node( integer i ) const { return m_X[size_t(i)]; }

    //!
    //! Return the i-th node of the spline (y component).
    //!
    real_type y_node( integer i ) const { return m_Y[size_t(i)]; }

    //!
    //! Return the i-th node of the spline (y component).
    //!
    real_type
    z_node( integer i, integer j ) const
    { return m_Z[size_t(this->ipos_C(i,j))]; }

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

    void
    build(
      real_type const * x, integer incx,
      real_type const * y, integer incy,
      real_type const * z, integer ldZ,
      integer nx, integer ny,
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
      bool xyswp = (fortran_storage && transposed) ||
                   (!fortran_storage && !transposed);
      this->build(
        &x.front(), 1,
        &y.front(), 1,
        &z.front(), integer(xyswp ? y.size() : x.size()),
        integer(x.size()), integer(y.size()),
        fortran_storage, transposed
      );
    }

    void
    build(
      real_type const * z,
      integer           ldZ,
      integer           nx,
      integer           ny,
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
      integer                   nx,
      integer                   ny,
      bool fortran_storage = false,
      bool transposed      = false
    ) {
      this->build( &z.front(), nx, ny, fortran_storage ? nx : ny, fortran_storage, transposed );
    }

    void
    setup( GenericContainer const & gc );

    void
    build( GenericContainer const & gc )
    { setup(gc); }

    ///@}

    //!
    //! \name Evaluate
    //!
    ///@{

    //!
    //! Evaluate spline value at point \f$ (x,y) \f$.
    //!
    virtual
    real_type
    operator () ( real_type x, real_type y ) const = 0;

    //!
    //! Value and first derivatives at point \f$ (x,y) \f$:
    //!
    //! - d[0] value of the spline \f$ S(x,y) \f$
    //! - d[1] derivative respect to \f$ x \f$ of the spline: \f$ S_x(x,y) \f$
    //! - d[2] derivative respect to \f$ y \f$ of the spline: \f$ S_y(x,y) \f$
    //!
    virtual
    void
    D( real_type x, real_type y, real_type d[3] ) const = 0;

    //!
    //! First derivatives respect to \f$ x \f$ at point \f$ (x,y) \f$
    //! of the spline: \f$ S_x(x,y) \f$.
    //!
    virtual
    real_type
    Dx( real_type x, real_type y ) const = 0;

    //!
    //! First derivatives respect to \f$ y \f$ at point \f$ (x,y) \f$
    //! of the spline: \f$ S_y(x,y) \f$.
    //!
    virtual
    real_type
    Dy( real_type x, real_type y ) const = 0;

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
    virtual
    void
    DD( real_type x, real_type y, real_type dd[6] ) const = 0;

    //!
    //! Second derivatives respect to \f$ x \f$ at point \f$ (x,y) \f$
    //! of the spline: \f$ S_{xx}(x,y) \f$.
    //!
    virtual
    real_type
    Dxx( real_type x, real_type y ) const = 0;

    //!
    //! Mixed second derivatives: \f$ S_{xy}(x,y) \f$.
    //!
    virtual
    real_type
    Dxy( real_type x, real_type y ) const = 0;

    //!
    //! Second derivatives respect to \f$ y \f$ at point \f$ (x,y) \f$
    //! of the spline: \f$ S_{yy}(x,y) \f$.
    //!
    virtual
    real_type
    Dyy( real_type x, real_type y ) const = 0;

    //!
    //! Evaluate spline value at point \f$ (x,y) \f$.
    //!
    real_type
    eval( real_type x, real_type y ) const
    { return (*this)(x,y); }

    //!
    //! Alias for `Dx(x,y)`
    //!
    real_type
    eval_D_1( real_type x, real_type y ) const
    { return this->Dx(x,y); }

    //!
    //! Alias for `Dy(x,y)`
    //!
    real_type
    eval_D_2( real_type x, real_type y ) const
    { return this->Dy(x,y); }

    //!
    //! Alias for `Dxx(x,y)`
    //!
    real_type
    eval_D_1_1( real_type x, real_type y ) const
    { return this->Dxx(x,y); }

    //!
    //! Alias for `Dxy(x,y)`
    //!
    real_type
    eval_D_1_2( real_type x, real_type y ) const
    { return this->Dxy(x,y); }

    //!
    //! Alias for `Dyy(x,y)`
    //!
    real_type
    eval_D_2_2( real_type x, real_type y ) const
    { return this->Dyy(x,y); }

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
    //! Return a string with information about the spline.
    //!
    virtual string info() const;

    //!
    //! write a string with information about the spline.
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
    real_type xNode( integer i ) const { return m_X[size_t(i)]; }
    real_type yNode( integer i ) const { return m_Y[size_t(i)]; }
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

#ifndef DOXYGEN_SHOULD_SKIP_THIS

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

#endif

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
#ifdef __clang__
#pragma clang diagnostic pop
#endif

#endif
