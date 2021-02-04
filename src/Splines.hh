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

//! Various kind of splines
namespace Splines {

  using std::vector;
  using std::string;
  using std::exception;
  using std::runtime_error;
  using std::basic_ostream;
  using std::ostringstream;
  using std::lower_bound;
  using std::pair;
  using std::cerr;

  typedef double real_type; //!< Floating point type for splines
  typedef int    integer;   //!< Signed integer type for splines
  typedef basic_ostream<char> ostream_type;

  void backtrace( ostream_type & );

  //! Associate a number for each type of splines implemented
  typedef enum {
    CONSTANT_TYPE   = 0,
    LINEAR_TYPE     = 1,
    CUBIC_TYPE      = 2,
    AKIMA_TYPE      = 3,
    BESSEL_TYPE     = 4,
    PCHIP_TYPE      = 5,
    QUINTIC_TYPE    = 6,
    HERMITE_TYPE    = 7,
    SPLINE_SET_TYPE = 8,
    SPLINE_VEC_TYPE = 9
  } SplineType1D;

  //! Associate a number for each type of splines implemented
  typedef enum {
    BILINEAR_TYPE  = 0,
    BICUBIC_TYPE   = 1,
    BIQUINTIC_TYPE = 2,
    AKIMA2D_TYPE   = 3
  } SplineType2D;

  extern char const *spline_type_1D[];

  extern SplineType1D string_to_splineType( string const & n );

  using GenericContainerNamespace::GenericContainer;
  using GenericContainerNamespace::vec_real_type;
  using GenericContainerNamespace::vec_string_type;
  using GenericContainerNamespace::vector_type;
  using GenericContainerNamespace::map_type;

  /*
  //   _   _                     _ _
  //  | | | | ___ _ __ _ __ ___ (_) |_ ___
  //  | |_| |/ _ \ '__| '_ ` _ \| | __/ _ \
  //  |  _  |  __/ |  | | | | | | | ||  __/
  //  |_| |_|\___|_|  |_| |_| |_|_|\__\___|
  */
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

  //! Check if cubic spline with this data is monotone, -1 no, 0 yes, 1 strictly monotone
  integer
  checkCubicSplineMonotonicity(
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

  void
  uniform(
    integer         dim,
    integer         npts,
    real_type const pnts[],
    integer         ld_pnts,
    real_type       t[]
  );

  void
  chordal(
    integer         dim,
    integer         npts,
    real_type const pnts[],
    integer         ld_pnts,
    real_type       t[]
  );

  void
  centripetal(
    integer         dim,
    integer         npts,
    real_type const pnts[],
    integer         ld_pnts,
    real_type const alpha,
    real_type       t[]
  );

  void
  universal(
    integer         dim,
    integer         npts,
    real_type const pnts[],
    integer         ld_pnts,
    real_type       t[]
  );

  void
  FoleyNielsen(
    integer         dim,
    integer         npts,
    real_type const pnts[],
    integer         ld_pnts,
    real_type       t[]
  );

  void
  FangHung(
    integer         dim,
    integer         npts,
    real_type const pnts[],
    integer         ld_pnts,
    real_type       t[]
  );

  /*\
   |   ____        _ _
   |  / ___| _ __ | (_)_ __   ___
   |  \___ \| '_ \| | | '_ \ / _ \
   |   ___) | |_) | | | | | |  __/
   |  |____/| .__/|_|_|_| |_|\___|
   |        |_|
  \*/
  //! Spline Management Class
  class Spline {
  protected:

    string m_name;
    bool   m_curve_is_closed;
    bool   m_curve_can_extend;
    bool   m_curve_extended_constant;

    integer   m_npts;
    integer   m_npts_reserved;
    real_type *m_X; // allocated in the derived class!
    real_type *m_Y; // allocated in the derived class!

    mutable Utils::BinarySearch<integer> m_bs;

    void initLastInterval();

    Spline( Spline const & ) = delete;
    Spline const & operator = ( Spline const & ) = delete;

  public:

    //! spline constructor
    Spline( string const & name = "Spline" )
    : m_name(name)
    , m_curve_is_closed(false)
    , m_curve_can_extend(true)
    , m_curve_extended_constant(false)
    , m_npts(0)
    , m_npts_reserved(0)
    , m_X(nullptr)
    , m_Y(nullptr)
    {
      this->initLastInterval();
    }
    //! spline destructor
    virtual
    ~Spline()
    {}

    integer search( real_type & x ) const;

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

    //! return the number of support points of the spline.
    integer
    numPoints(void) const { return m_npts; }

    //! return the i-th node of the spline (x component).
    real_type
    xNode( integer i ) const { return m_X[size_t(i)]; }

    //! return the i-th node of the spline (y component).
    real_type
    yNode( integer i ) const { return m_Y[size_t(i)]; }

    //! return first node of the spline (x component).
    real_type
    xBegin() const { return m_X[0]; }

    //! return first node of the spline (y component).
    real_type
    yBegin() const { return m_Y[0]; }

    //! return last node of the spline (x component).
    real_type
    xEnd() const { return m_X[size_t(m_npts-1)]; }

    //! return last node of the spline (y component).
    real_type
    yEnd() const { return m_Y[size_t(m_npts-1)]; }

    //! Allocate memory for `npts` points
    virtual
    void
    reserve( integer npts ) UTILS_PURE_VIRTUAL;

    //! Add a support point (x,y) to the spline.
    void pushBack( real_type x, real_type y );

    //! Drop a support point to the spline.
    void dropBack() { if ( m_npts > 0 ) --m_npts; }

    //! Build a spline.
    // must be defined in derived classes
    virtual
    void
    build(void) UTILS_PURE_VIRTUAL;

    virtual
    void
    setup( GenericContainer const & gc );

    void
    build( GenericContainer const & gc )
    { setup(gc); }

    //! Build a spline.
    /*!
     | \param x    vector of x-coordinates
     | \param incx access elements as x[0], x[incx], x[2*incx],...
     | \param y    vector of y-coordinates
     | \param incy access elements as y[0], y[incy], y[2*incy],...
     | \param n    total number of points
    \*/
    // must be defined in derived classes
    virtual
    void
    build(
      real_type const x[], integer incx,
      real_type const y[], integer incy,
      integer n
    );

    //! Build a spline.
    /*!
     | \param x vector of x-coordinates
     | \param y vector of y-coordinates
     | \param n total number of points
    \*/
    inline
    void
    build( real_type const x[], real_type const y[], integer n )
    { this->build( x, 1, y, 1, n ); }

    //! Build a spline.
    /*!
     | \param x vector of x-coordinates
     | \param y vector of y-coordinates
    \*/
    inline
    void
    build( vector<real_type> const & x, vector<real_type> const & y ) {
      integer N = integer(x.size());
      if ( N > integer(y.size()) ) N = integer(y.size());
      this->build( &x.front(), 1, &y.front(), 1, N );
    }

    //! Cancel the support points, empty the spline.
    virtual
    void
    clear(void) UTILS_PURE_VIRTUAL;

    //! return x-minumum spline value
    real_type xMin() const { return m_X[0]; }

    //! return x-maximum spline value
    real_type xMax() const { return m_X[m_npts-1]; }

    //! return y-minumum spline value
    real_type
    yMin() const {
      integer N = m_npts;
      if ( type() == CONSTANT_TYPE ) --N;
      return *std::min_element(m_Y,m_Y+N);
    }

    //! return y-maximum spline value
    real_type
    yMax() const {
      integer N = m_npts;
      if ( type() == CONSTANT_TYPE ) --N;
      return *std::max_element(m_Y,m_Y+N);
    }

    ///////////////////////////////////////////////////////////////////////////
    //! change X-origin of the spline
    void setOrigin( real_type x0 );

    //! change X-range of the spline
    void setRange( real_type xmin, real_type xmax );

    ///////////////////////////////////////////////////////////////////////////
    //! dump a sample of the spline
    void
    dump(
      ostream_type & s,
      integer        nintervals,
      char const     header[] = "x\ty"
    ) const;

    void
    dump(
      char const fname[],
      integer    nintervals,
      char const header[] = "x\ty"
    ) const {
      std::ofstream file(fname);
      this->dump( file, nintervals, header );
      file.close();
    }

    ///////////////////////////////////////////////////////////////////////////
    //! Evaluate spline value
    virtual
    real_type
    operator () ( real_type x ) const UTILS_PURE_VIRTUAL;

    //! First derivative
    virtual
    real_type
    D( real_type x ) const UTILS_PURE_VIRTUAL;

    //! Second derivative
    virtual
    real_type
    DD( real_type x ) const UTILS_PURE_VIRTUAL;

    //! Third derivative
    virtual
    real_type
    DDD( real_type x ) const UTILS_PURE_VIRTUAL;

    //! 4th derivative
    virtual
    real_type
    DDDD( real_type ) const
    { return real_type(0); }

    //! 4th derivative
    virtual
    real_type
    DDDDD( real_type ) const
    { return real_type(0); }

    //! Some aliases
    real_type eval( real_type x ) const { return (*this)(x); }
    real_type eval_D( real_type x ) const { return this->D(x); }
    real_type eval_DD( real_type x ) const { return this->DD(x); }
    real_type eval_DDD( real_type x ) const { return this->DDD(x); }
    real_type eval_DDDD( real_type x ) const { return this->DDDD(x); }
    real_type eval_DDDDD( real_type x ) const { return this->DDDDD(x); }

    ///////////////////////////////////////////////////////////////////////////
    //! Evaluate spline value when interval is known
    virtual
    real_type
    id_eval( integer ni, real_type x ) const UTILS_PURE_VIRTUAL;

    //! First derivative
    virtual
    real_type
    id_D( integer ni, real_type x ) const UTILS_PURE_VIRTUAL;

    //! Second derivative
    virtual
    real_type
    id_DD( integer ni, real_type x ) const UTILS_PURE_VIRTUAL;

    //! Third derivative
    virtual
    real_type
    id_DDD( integer ni, real_type x ) const UTILS_PURE_VIRTUAL;

    //! 4th derivative
    virtual
    real_type
    id_DDDD( integer, real_type ) const
    { return real_type(0); }

    //! 4th derivative
    virtual
    real_type
    id_DDDDD( integer, real_type ) const
    { return real_type(0); }

    //! get the piecewise polinomials of the spline
    virtual
    integer // order
    coeffs(
      real_type cfs[],
      real_type nodes[],
      bool      transpose = false
    ) const UTILS_PURE_VIRTUAL;

    virtual
    integer // order
    order() const UTILS_PURE_VIRTUAL;

    //! Print spline coefficients
    virtual
    void
    writeToStream( ostream_type & s ) const UTILS_PURE_VIRTUAL;

    //! Return spline typename
    char const *
    type_name() const
    { return Splines::spline_type_1D[type()]; }

    //! Return spline type (as number)
    virtual
    unsigned
    type() const UTILS_PURE_VIRTUAL;

    void
    info( ostream_type & s ) const;

    friend class SplineSet;

  };

  //! compute curvature of a planar curve
  real_type
  curvature( real_type s, Spline const & X, Spline const & Y );

  //! compute curvature derivative of a planar curve
  real_type
  curvature_D( real_type s, Spline const & X, Spline const & Y );

  //! compute curvature second derivative of a planar curve
  real_type
  curvature_DD( real_type s, Spline const & X, Spline const & Y );

  /*\
   |    ____      _     _        ____        _ _              ____
   |   / ___|   _| |__ (_) ___  / ___| _ __ | (_)_ __   ___  | __ )  __ _ ___  ___
   |  | |  | | | | '_ \| |/ __| \___ \| '_ \| | | '_ \ / _ \ |  _ \ / _` / __|/ _ \
   |  | |__| |_| | |_) | | (__   ___) | |_) | | | | | |  __/ | |_) | (_| \__ \  __/
   |   \____\__,_|_.__/|_|\___| |____/| .__/|_|_|_| |_|\___| |____/ \__,_|___/\___|
   |                                  |_|
  \*/
  //! cubic spline base class
  class CubicSplineBase : public Spline {
  protected:
    Utils::Malloc<real_type> m_baseValue;

    real_type * m_Yp;
    bool        m_external_alloc;

  public:

    using Spline::build;

    //! spline constructor
    CubicSplineBase( string const & name = "CubicSplineBase")
    : Spline(name)
    , m_baseValue(name+"_memory")
    , m_Yp(nullptr)
    , m_external_alloc(false)
    {}

    virtual
    ~CubicSplineBase() UTILS_OVERRIDE
    {}

    void
    copySpline( CubicSplineBase const & S );

    //! return the i-th node of the spline (y' component).
    real_type
    ypNode( integer i ) const
    { return m_Yp[size_t(i)]; }

    //! change X-range of the spline
    void
    setRange( real_type xmin, real_type xmax );

    //! Use externally allocated memory for `npts` points
    void
    reserve_external(
      integer       n,
      real_type * & p_x,
      real_type * & p_y,
      real_type * & p_dy
    );

    // --------------------------- VIRTUALS -----------------------------------
    //! Evaluate spline value
    virtual
    real_type
    operator () ( real_type x ) const UTILS_OVERRIDE;

    //! First derivative
    virtual
    real_type
    D( real_type x ) const UTILS_OVERRIDE;

    //! Second derivative
    virtual
    real_type
    DD( real_type x ) const UTILS_OVERRIDE;

    //! Third derivative
    virtual
    real_type
    DDD( real_type x ) const UTILS_OVERRIDE;

    //! Evaluate spline value knowing interval
    virtual
    real_type
    id_eval( integer ni, real_type x ) const UTILS_OVERRIDE;

    //! First derivative
    virtual
    real_type
    id_D( integer ni, real_type x ) const UTILS_OVERRIDE;

    //! Second derivative
    virtual
    real_type
    id_DD( integer ni, real_type x ) const UTILS_OVERRIDE;

    //! Third derivative
    virtual
    real_type
    id_DDD( integer ni, real_type x ) const UTILS_OVERRIDE;

    //! Print spline coefficients
    virtual
    void
    writeToStream( ostream_type & s ) const UTILS_OVERRIDE;

    // --------------------------- VIRTUALS -----------------------------------

    //! Allocate memory for `npts` points
    virtual
    void
    reserve( integer npts ) UTILS_OVERRIDE;

    //! Build a spline.
    /*!
     | \param x     vector of x-coordinates
     | \param incx  access elements as x[0], x[incx], x[2*incx],...
     | \param y     vector of y-coordinates
     | \param incy  access elements as y[0], y[incy], y[2*incy],...
     | \param yp    vector of y-defivative
     | \param incyp access elements as yp[0], yp[incyp], yp[2*incyy],...
     | \param n     total number of points
    \*/
    // must be defined in derived classes
    void
    build(
      real_type const x[],  integer incx,
      real_type const y[],  integer incy,
      real_type const yp[], integer incyp,
      integer n
    );

    //! Build a spline.
    /*!
     | \param x  vector of x-coordinates
     | \param y  vector of y-coordinates
     | \param yp vector of y'-coordinates
     | \param n  total number of points
    \*/
    inline
    void
    build(
      real_type const x[],
      real_type const y[],
      real_type const yp[],
      integer n
    ) {
      this->build( x, 1, y, 1, yp, 1, n );
    }

    //! Build a spline.
    /*!
     | \param x  vector of x-coordinates
     | \param y  vector of y-coordinates
     | \param yp vector of y'-coordinates
    \*/
    void
    build(
      vector<real_type> const & x,
      vector<real_type> const & y,
      vector<real_type> const & yp
    );

    //! Cancel the support points, empty the spline.
    virtual
    void
    clear(void) UTILS_OVERRIDE;

    //! get the piecewise polinomials of the spline
    virtual
    integer // order
    coeffs(
      real_type cfs[],
      real_type nodes[],
      bool      transpose = false
    ) const UTILS_OVERRIDE;

    virtual
    integer // order
    order() const UTILS_OVERRIDE;

  };

  /*\
   |    ____      _     _      ____        _ _
   |   / ___|   _| |__ (_) ___/ ___| _ __ | (_)_ __   ___
   |  | |  | | | | '_ \| |/ __\___ \| '_ \| | | '_ \ / _ \
   |  | |__| |_| | |_) | | (__ ___) | |_) | | | | | |  __/
   |   \____\__,_|_.__/|_|\___|____/| .__/|_|_|_| |_|\___|
   |                                |_|
  \*/

  typedef enum {
    EXTRAPOLATE_BC = 0,
    NATURAL_BC,
    PARABOLIC_RUNOUT_BC,
    NOT_A_KNOT
  } CUBIC_SPLINE_TYPE_BC;

  void
  CubicSpline_build(
    real_type const      X[],
    real_type const      Y[],
    real_type            Yp[],
    integer              npts,
    CUBIC_SPLINE_TYPE_BC bc0,
    CUBIC_SPLINE_TYPE_BC bcn
  );

  void
  CubicSpline_build(
    real_type const      X[],
    real_type const      Y[],
    real_type            Yp[],
    real_type            Ypp[],
    real_type            L[],
    real_type            D[],
    real_type            U[],
    integer              npts,
    CUBIC_SPLINE_TYPE_BC bc0,
    CUBIC_SPLINE_TYPE_BC bcn
  );

  //! Cubic Spline Management Class
  class CubicSpline : public CubicSplineBase {
  private:
    CUBIC_SPLINE_TYPE_BC m_bc0, m_bcn;
  public:

    using CubicSplineBase::build;
    using CubicSplineBase::reserve;

    //! spline constructor
    CubicSpline( string const & name = "CubicSpline" )
    : CubicSplineBase( name )
    , m_bc0( EXTRAPOLATE_BC )
    , m_bcn( EXTRAPOLATE_BC )
    {}

    //! spline destructor
    virtual
    ~CubicSpline() UTILS_OVERRIDE
    {}

    /*!
     | \param bc0  initial boundary condition.
    \*/
    void
    setInitialBC( CUBIC_SPLINE_TYPE_BC bc0 )
    { m_bc0 = bc0; }

    /*!
     | \param bcn final boundary condition.
    \*/
    void
    setFinalBC( CUBIC_SPLINE_TYPE_BC bcn )
    { m_bcn = bcn; }

    //! Return spline type (as number)
    virtual
    unsigned
    type() const UTILS_OVERRIDE
    { return CUBIC_TYPE; }

    // --------------------------- VIRTUALS -----------------------------------

    virtual
    void
    build(void) UTILS_OVERRIDE;

    virtual
    void
    setup( GenericContainer const & gc ) UTILS_OVERRIDE;

  };

  /*\
   |      _    _    _                   ____        _ _
   |     / \  | | _(_)_ __ ___   __ _  / ___| _ __ | (_)_ __   ___
   |    / _ \ | |/ / | '_ ` _ \ / _` | \___ \| '_ \| | | '_ \ / _ \
   |   / ___ \|   <| | | | | | | (_| |  ___) | |_) | | | | | |  __/
   |  /_/   \_\_|\_\_|_| |_| |_|\__,_| |____/| .__/|_|_|_| |_|\___|
   |                                         |_|
  \*/

  void
  Akima_build(
    real_type const X[],
    real_type const Y[],
    real_type       Yp[],
    integer         npts
  );

  //! Akima spline class
  /*!
   |  Reference
   |  =========
   |  Hiroshi Akima, Journal of the ACM, Vol. 17, No. 4, October 1970, pages 589-602.
  \*/
  class AkimaSpline : public CubicSplineBase {
  public:

    using CubicSplineBase::build;
    using CubicSplineBase::reserve;

    //! spline constructor
    AkimaSpline( string const & name = "AkimaSpline" )
    : CubicSplineBase( name )
    {}

    //! spline destructor
    virtual
    ~AkimaSpline() UTILS_OVERRIDE
    {}

    //! Return spline type (as number)
    virtual
    unsigned
    type() const UTILS_OVERRIDE
    { return AKIMA_TYPE; }

    // --------------------------- VIRTUALS -----------------------------------

    //! Build an Akima spline from previously inserted points
    virtual
    void
    build(void) UTILS_OVERRIDE;

    virtual
    void
    setup( GenericContainer const & gc ) UTILS_OVERRIDE;

  };

  /*
   |   ____                     _ ____        _ _
   |  | __ )  ___  ___ ___  ___| / ___| _ __ | (_)_ __   ___
   |  |  _ \ / _ \/ __/ __|/ _ \ \___ \| '_ \| | | '_ \ / _ \
   |  | |_) |  __/\__ \__ \  __/ |___) | |_) | | | | | |  __/
   |  |____/ \___||___/___/\___|_|____/| .__/|_|_|_| |_|\___|
   |                                   |_|
  \*/

  void
  Bessel_build(
    real_type const X[],
    real_type const Y[],
    real_type       Yp[],
    integer         npts
  );

  //! Bessel spline class
  class BesselSpline : public CubicSplineBase {
  public:

    using CubicSplineBase::build;
    using CubicSplineBase::reserve;

    //! spline constructor
    BesselSpline( string const & name = "BesselSpline" )
    : CubicSplineBase( name )
    {}

    //! spline destructor
    virtual
    ~BesselSpline() UTILS_OVERRIDE
    {}

    //! Return spline type (as number)
    virtual
    unsigned
    type() const UTILS_OVERRIDE
    { return BESSEL_TYPE; }

    // --------------------------- VIRTUALS -----------------------------------

    //! Build a Bessel spline from previously inserted points
    virtual
    void
    build (void) UTILS_OVERRIDE;

    virtual
    void
    setup( GenericContainer const & gc ) UTILS_OVERRIDE;
  };

  /*\
   |   ____      _     _      ____        _ _
   |  |  _ \ ___| |__ (_)_ __/ ___| _ __ | (_)_ __   ___
   |  | |_) / __| '_ \| | '_ \___ \| '_ \| | | '_ \ / _ \
   |  |  __/ (__| | | | | |_) |__) | |_) | | | | | |  __/
   |  |_|   \___|_| |_|_| .__/____/| .__/|_|_|_| |_|\___|
   |                    |_|        |_|
  \*/
  void
  Pchip_build(
    real_type const X[],
    real_type const Y[],
    real_type       Yp[],
    integer         npts
  );

  //! Pchip (Piecewise Cubic Hermite Interpolating Polynomial) spline class
  class PchipSpline : public CubicSplineBase {
  public:

    using CubicSplineBase::build;
    using CubicSplineBase::reserve;

    //! spline constructor
    PchipSpline( string const & name = "PchipSpline" )
    : CubicSplineBase( name )
    {}

    //! spline destructor
    virtual
    ~PchipSpline() UTILS_OVERRIDE
    {}

    //! Return spline type (as number)
    virtual
    unsigned
    type() const UTILS_OVERRIDE
    { return PCHIP_TYPE; }

    // --------------------------- VIRTUALS -----------------------------------

    //! Build a Monotone spline from previously inserted points
    virtual
    void
    build(void) UTILS_OVERRIDE;

    virtual
    void
    setup( GenericContainer const & gc ) UTILS_OVERRIDE;

  };

  /*\
   |   _     _                       ____        _ _
   |  | |   (_)_ __   ___  __ _ _ __/ ___| _ __ | (_)_ __   ___
   |  | |   | | '_ \ / _ \/ _` | '__\___ \| '_ \| | | '_ \ / _ \
   |  | |___| | | | |  __/ (_| | |   ___) | |_) | | | | | |  __/
   |  |_____|_|_| |_|\___|\__,_|_|  |____/| .__/|_|_|_| |_|\___|
   |                                      |_|
  \*/
  //! Linear spline class
  class LinearSpline : public Spline {
    Utils::Malloc<real_type> m_baseValue;
    bool                     m_external_alloc;
  public:

    using Spline::build;

    LinearSpline( string const & name = "LinearSpline" )
    : Spline(name)
    , m_baseValue( name+"_memory")
    , m_external_alloc(false)
    {
      m_curve_extended_constant = true; // by default linear spline extend constant
    }

    virtual
    ~LinearSpline() UTILS_OVERRIDE
    {}

    //! Use externally allocated memory for `npts` points
    void
    reserve_external(
      integer      n,
      real_type *& p_x,
      real_type *& p_y
    );

    // --------------------------- VIRTUALS -----------------------------------

    //! Evalute spline value at `x`
    virtual
    real_type
    operator () ( real_type x ) const UTILS_OVERRIDE;

    //! First derivative
    virtual
    real_type
    D( real_type x ) const UTILS_OVERRIDE;

    //! Second derivative
    virtual
    real_type
    DD( real_type ) const UTILS_OVERRIDE
    { return 0; }

    //! Third derivative
    virtual
    real_type
    DDD( real_type ) const UTILS_OVERRIDE
    { return 0; }

    //! Evaluate spline value knowing interval
    virtual
    real_type
    id_eval( integer ni, real_type x ) const UTILS_OVERRIDE;

    //! First derivative
    virtual
    real_type
    id_D( integer, real_type ) const UTILS_OVERRIDE;

    //! Second derivative
    virtual
    real_type
    id_DD( integer, real_type ) const UTILS_OVERRIDE
    { return 0; }

    //! Third derivative
    virtual
    real_type
    id_DDD( integer, real_type ) const UTILS_OVERRIDE
    { return 0; }

    //! Print spline coefficients
    virtual
    void
    writeToStream( ostream_type & s ) const UTILS_OVERRIDE;

    //! Return spline type (as number)
    virtual
    unsigned
    type() const UTILS_OVERRIDE
    { return LINEAR_TYPE; }

    // --------------------------- VIRTUALS -----------------------------------

    //! Allocate memory for `npts` points
    virtual
    void
    reserve( integer npts ) UTILS_OVERRIDE;

    //! added for compatibility with cubic splines
    virtual
    void
    build(void) UTILS_OVERRIDE
    {}

    //! Cancel the support points, empty the spline.
    virtual
    void
    clear(void) UTILS_OVERRIDE;

    //! get the piecewise polinomials of the spline
    virtual
    integer // order
    coeffs(
      real_type cfs[],
      real_type nodes[],
      bool      transpose = false
    ) const UTILS_OVERRIDE;

    virtual
    integer // order
    order() const UTILS_OVERRIDE;

    virtual
    void
    setup( GenericContainer const & gc ) UTILS_OVERRIDE;

  };

  /*\
   |    ____                _              _       ____        _ _
   |   / ___|___  _ __  ___| |_ __ _ _ __ | |_ ___/ ___| _ __ | (_)_ __   ___
   |  | |   / _ \| '_ \/ __| __/ _` | '_ \| __/ __\___ \| '_ \| | | '_ \ / _ \
   |  | |__| (_) | | | \__ \ || (_| | | | | |_\__ \___) | |_) | | | | | |  __/
   |   \____\___/|_| |_|___/\__\__,_|_| |_|\__|___/____/| .__/|_|_|_| |_|\___|
   |                                                    |_|
  \*/
  //! Picewise constants spline class
  class ConstantSpline : public Spline {
    Utils::Malloc<real_type> m_baseValue;
    bool                     m_external_alloc;
  public:

    using Spline::build;

    ConstantSpline( string const & name = "ConstantSpline" )
    : Spline(name)
    , m_baseValue(name+"_memory")
    , m_external_alloc(false)
    {}

    ~ConstantSpline() UTILS_OVERRIDE
    {}

    //! Use externally allocated memory for `npts` points
    void
    reserve_external(
      integer       n,
      real_type * & p_x,
      real_type * & p_y
    );

    // --------------------------- VIRTUALS -----------------------------------
    //! Build a spline.
    virtual
    void
    build(void) UTILS_OVERRIDE
    {} // nothing to do

    virtual
    void
    build(
      real_type const x[], integer incx,
      real_type const y[], integer incy,
      integer n
    ) UTILS_OVERRIDE;

    //! Evaluate spline value at `x`
    virtual
    real_type
    operator () ( real_type x ) const UTILS_OVERRIDE;

    //! First derivative
    virtual
    real_type
    D( real_type ) const UTILS_OVERRIDE
    { return 0; }

    //! Second derivative
    virtual
    real_type
    DD( real_type ) const UTILS_OVERRIDE
    { return 0; }

    //! Third derivative
    virtual
    real_type
    DDD( real_type ) const UTILS_OVERRIDE
    { return 0; }

    //! Evaluate spline value at `x` knowing interval
    virtual
    real_type
    id_eval( integer ni, real_type x ) const UTILS_OVERRIDE;

    //! First derivative
    virtual
    real_type
    id_D( integer, real_type ) const UTILS_OVERRIDE
    { return 0; }

    //! Second derivative
    virtual
    real_type
    id_DD( integer, real_type ) const UTILS_OVERRIDE
    { return 0; }

    //! Third derivative
    virtual
    real_type
    id_DDD( integer, real_type ) const UTILS_OVERRIDE
    { return 0; }

    //! Print spline coefficients
    virtual
    void
    writeToStream( ostream_type & ) const UTILS_OVERRIDE;

    //! Return spline type (as number)
    virtual
    unsigned
    type() const UTILS_OVERRIDE
    { return CONSTANT_TYPE; }

    // --------------------------- VIRTUALS -----------------------------------

    //! Allocate memory for `npts` points
    virtual
    void
    reserve( integer npts ) UTILS_OVERRIDE;

    //! Cancel the support points, empty the spline.
    virtual
    void
    clear(void) UTILS_OVERRIDE;

    //! get the piecewise polinomials of the spline
    virtual
    integer // order
    coeffs(
      real_type cfs[],
      real_type nodes[],
      bool      transpose = false
    ) const UTILS_OVERRIDE;

    virtual
    integer // order
    order() const UTILS_OVERRIDE;

    virtual
    void
    setup( GenericContainer const & gc ) UTILS_OVERRIDE;

  };

  /*\
   |    _   _                     _ _       ____        _ _
   |   | | | | ___ _ __ _ __ ___ (_) |_ ___/ ___| _ __ | (_)_ __   ___
   |   | |_| |/ _ \ '__| '_ ` _ \| | __/ _ \___ \| '_ \| | | '_ \ / _ \
   |   |  _  |  __/ |  | | | | | | | ||  __/___) | |_) | | | | | |  __/
   |   |_| |_|\___|_|  |_| |_| |_|_|\__\___|____/| .__/|_|_|_| |_|\___|
   |                                             |_|
  \*/
  //! Hermite Spline Management Class
  class HermiteSpline : public CubicSplineBase {
  public:

    using CubicSplineBase::build;
    using CubicSplineBase::reserve;

    //! spline constructor
    HermiteSpline( string const & name = "HermiteSpline" )
    : CubicSplineBase( name )
    {}

    //! spline destructor
    virtual
    ~HermiteSpline() UTILS_OVERRIDE
    {}

    //! Return spline type (as number)
    virtual
    unsigned
    type() const UTILS_OVERRIDE
    { return HERMITE_TYPE; }

    // --------------------------- VIRTUALS -----------------------------------

    virtual
    void
    build(void) UTILS_OVERRIDE
    {} // nothing to do

    // block method!
    virtual
    void
    build(
      real_type const [], integer,
      real_type const [], integer,
      integer
    ) UTILS_OVERRIDE;

    virtual
    void
    setup( GenericContainer const & gc ) UTILS_OVERRIDE;

  };

  /*\
   |    ___        _       _   _      ____        _ _            ____
   |   / _ \ _   _(_)_ __ | |_(_) ___/ ___| _ __ | (_)_ __   ___| __ )  __ _ ___  ___
   |  | | | | | | | | '_ \| __| |/ __\___ \| '_ \| | | '_ \ / _ \  _ \ / _` / __|/ _ \
   |  | |_| | |_| | | | | | |_| | (__ ___) | |_) | | | | | |  __/ |_) | (_| \__ \  __/
   |   \__\_\\__,_|_|_| |_|\__|_|\___|____/| .__/|_|_|_| |_|\___|____/ \__,_|___/\___|
   |                                       |_|
   |
  \*/
  //! cubic quintic base class
  class QuinticSplineBase : public Spline {
  protected:
    Utils::Malloc<real_type> m_baseValue;

    real_type * m_Yp;
    real_type * m_Ypp;
    bool        m_external_alloc;

  public:

    using Spline::build;

    //! spline constructor
    QuinticSplineBase( string const & name = "Spline" )
    : Spline(name)
    , m_baseValue(name+"_memeory")
    , m_Yp(nullptr)
    , m_Ypp(nullptr)
    , m_external_alloc(false)
    {}

    virtual
    ~QuinticSplineBase() UTILS_OVERRIDE
    {}

    void
    copySpline( QuinticSplineBase const & S );

    //! return the i-th node of the spline (y' component).
    real_type
    ypNode( integer i ) const
    { return m_Yp[size_t(i)]; }

    //! return the i-th node of the spline (y'' component).
    real_type
    yppNode( integer i ) const
    { return m_Ypp[size_t(i)]; }

    //! change X-range of the spline
    void
    setRange( real_type xmin, real_type xmax );

    //! Use externally allocated memory for `npts` points
    void
    reserve_external(
      integer       n,
      real_type * & p_x,
      real_type * & p_y,
      real_type * & p_Yp,
      real_type * & p_Ypp
    );

    // --------------------------- VIRTUALS -----------------------------------

    //! Evaluate spline value
    virtual
    real_type
    operator () ( real_type x ) const UTILS_OVERRIDE;

    //! First derivative
    virtual
    real_type
    D( real_type x ) const UTILS_OVERRIDE;

    //! Second derivative
    virtual
    real_type
    DD( real_type x ) const UTILS_OVERRIDE;

    //! Third derivative
    virtual
    real_type
    DDD( real_type x ) const UTILS_OVERRIDE;

    //! Fourth derivative
    virtual
    real_type
    DDDD( real_type x ) const UTILS_OVERRIDE;

    //! Fifth derivative
    virtual
    real_type
    DDDDD( real_type x ) const UTILS_OVERRIDE;

    //! Evaluate spline value knowing interval
    virtual
    real_type
    id_eval( integer ni, real_type x ) const UTILS_OVERRIDE;

    //! First derivative
    virtual
    real_type
    id_D( integer ni, real_type x ) const UTILS_OVERRIDE;

    //! Second derivative
    virtual
    real_type
    id_DD( integer ni, real_type x ) const UTILS_OVERRIDE;

    //! Third derivative
    virtual
    real_type
    id_DDD( integer ni, real_type x ) const UTILS_OVERRIDE;

    //! Fourth derivative
    virtual
    real_type
    id_DDDD( integer ni, real_type x ) const UTILS_OVERRIDE;

    //! Fifth derivative
    virtual
    real_type
    id_DDDDD( integer ni, real_type x ) const UTILS_OVERRIDE;

    //! Print spline coefficients
    virtual
    void
    writeToStream( ostream_type & s ) const UTILS_OVERRIDE;

    //! Return spline type (as number)
    virtual
    unsigned
    type() const UTILS_OVERRIDE
    { return QUINTIC_TYPE; }

    //! Allocate memory for `npts` points
    virtual
    void
    reserve( integer npts ) UTILS_OVERRIDE;

    //! Cancel the support points, empty the spline.
    virtual
    void
    clear(void) UTILS_OVERRIDE;

    //! get the piecewise polinomials of the spline
    virtual
    integer // order
    coeffs(
      real_type cfs[],
      real_type nodes[],
      bool      transpose = false
    ) const UTILS_OVERRIDE;

    virtual
    integer // order
    order() const UTILS_OVERRIDE;

  };

  /*\
   |    ___        _       _   _      ____        _ _
   |   / _ \ _   _(_)_ __ | |_(_) ___/ ___| _ __ | (_)_ __   ___
   |  | | | | | | | | '_ \| __| |/ __\___ \| '_ \| | | '_ \ / _ \
   |  | |_| | |_| | | | | | |_| | (__ ___) | |_) | | | | | |  __/
   |   \__\_\\__,_|_|_| |_|\__|_|\___|____/| .__/|_|_|_| |_|\___|
   |                                       |_|
   |
  \*/

  typedef enum {
    CUBIC_QUINTIC = 0,
    PCHIP_QUINTIC,
    AKIMA_QUINTIC,
    BESSEL_QUINTIC
  } QUINTIC_SPLINE_TYPE;

  //! Quintic spline class
  class QuinticSpline : public QuinticSplineBase {
    QUINTIC_SPLINE_TYPE m_q_sub_type;
  public:

    using Spline::build;
    using QuinticSplineBase::reserve;

    //! spline constructor
    QuinticSpline( string const & name = "Spline" )
    : QuinticSplineBase( name )
    , m_q_sub_type(CUBIC_QUINTIC)
    {}

    //! spline destructor
    virtual
    ~QuinticSpline() UTILS_OVERRIDE
    {}

    void
    setQuinticType( QUINTIC_SPLINE_TYPE qt )
    { m_q_sub_type = qt; }

    // --------------------------- VIRTUALS -----------------------------------
    //! Build a Monotone quintic spline from previously inserted points
    virtual
    void
    build(void) UTILS_OVERRIDE;

    virtual
    void
    setup( GenericContainer const & gc ) UTILS_OVERRIDE;

  };


  /*\
   |   ____        _ _            _ ____
   |  / ___| _ __ | (_)_ __   ___/ |  _ \
   |  \___ \| '_ \| | | '_ \ / _ \ | | | |
   |   ___) | |_) | | | | | |  __/ | |_| |
   |  |____/| .__/|_|_|_| |_|\___|_|____/
   |        |_|
  \*/
  //! Spline Management Class
  class Spline1D {
  protected:

    std::string m_name;

    Spline * m_pSpline;

    Spline1D( Spline1D const & ) = delete;
    Spline1D const & operator = ( Spline1D const & ) = delete;

  public:

    //! spline constructor
    Spline1D( std::string const & n )
    : m_name(n)
    , m_pSpline(nullptr)
    {}

    //! spline destructor
    ~Spline1D()
    {}

    string const & name() const { return m_pSpline->name(); }

    bool is_closed() const { return m_pSpline->is_closed(); }
    void make_closed() { m_pSpline->make_closed(); }
    void make_opened() { m_pSpline->make_opened(); }

    bool is_bounded() const { return m_pSpline->is_bounded(); }
    void make_unbounded()   { m_pSpline->make_unbounded(); }
    void make_bounded()     { m_pSpline->make_bounded(); }

    //! return the number of support points of the spline.
    integer numPoints(void) const { return m_pSpline->numPoints(); }

    //! return the i-th node of the spline (x component).
    real_type xNode( integer i ) const { return m_pSpline->xNode(i); }

    //! return the i-th node of the spline (y component).
    real_type yNode( integer i ) const { return m_pSpline->yNode(i); }

    //! return first node of the spline (x component).
    real_type xBegin() const { return m_pSpline->xBegin(); }

    //! return first node of the spline (y component).
    real_type yBegin() const { return m_pSpline->yBegin(); }

    //! return last node of the spline (x component).
    real_type xEnd() const { return m_pSpline->xEnd(); }

    //! return last node of the spline (y component).
    real_type yEnd() const { return m_pSpline->yEnd(); }

    //! Allocate memory for `npts` points
    void reserve( integer npts ) { return m_pSpline->reserve( npts ); }

    //! Add a support point (x,y) to the spline.
    void pushBack( real_type x, real_type y ) { return m_pSpline->pushBack( x, y ); }

    //! Drop a support point to the spline.
    void dropBack() { m_pSpline->dropBack(); }

    //! Build a spline.
    // must be defined in derived classes
    void build(void) { m_pSpline->build(); }

    void setup( GenericContainer const & gc );
    void build( GenericContainer const & gc ) { setup(gc); }

    //! Build a spline.
    /*!
     | \param x    vector of x-coordinates
     | \param incx access elements as x[0], x[incx], x[2*incx],...
     | \param y    vector of y-coordinates
     | \param incy access elements as y[0], y[incy], x[2*incy],...
     | \param n    total number of points
    \*/
    // must be defined in derived classes
    void
    build(
      SplineType1D tp,
      real_type const x[], integer incx,
      real_type const y[], integer incy,
      integer n
    );

    //! Build a spline.
    /*!
     | \param x vector of x-coordinates
     | \param y vector of y-coordinates
     | \param n total number of points
    \*/
    void
    build(
      SplineType1D    tp,
      real_type const x[],
      real_type const y[],
      integer         n
    ) {
      this->build( tp, x, 1, y, 1, n );
    }

    //! Build a spline.
    /*!
     | \param x vector of x-coordinates
     | \param y vector of y-coordinates
    \*/
    void
    build(
      SplineType1D              tp,
      vector<real_type> const & x,
      vector<real_type> const & y
    ) {
      integer n = integer(x.size());
      this->build( tp, &x.front(), &y.front(), n );
    }

    //! Cancel the support points, empty the spline.
    void
    clear(void) { m_pSpline->clear(); }

    //! return x-minumum spline value
    real_type xMin() const { return m_pSpline->xMin(); }

    //! return x-maximum spline value
    real_type xMax() const { return m_pSpline->xMax(); }

    //! return y-minumum spline value
    real_type yMin() const { return m_pSpline->yMin(); }

    //! return y-maximum spline value
    real_type yMax() const { return m_pSpline->yMax(); }

    ///////////////////////////////////////////////////////////////////////////
    //! change X-origin of the spline
    void
    setOrigin( real_type x0 )
    { return m_pSpline->setOrigin( x0 ); }

    //! change X-range of the spline
    void
    setRange( real_type xmin, real_type xmax )
    { return m_pSpline->setRange( xmin, xmax ); }

    ///////////////////////////////////////////////////////////////////////////
    //! dump a sample of the spline
    void
    dump(
      ostream_type & s,
      integer        nintervals,
      char const     header[] = "x\ty"
    ) const {
      m_pSpline->dump( s, nintervals, header );
    }

    void
    dump(
      char const fname[],
      integer    nintervals,
      char const header[] = "x\ty"
    ) const {
      m_pSpline->dump( fname, nintervals, header );
    }

    ///////////////////////////////////////////////////////////////////////////
    //! Evaluate spline value
    real_type
    operator () ( real_type x ) const { return (*m_pSpline)(x); }

    //! First derivative
    real_type
    D( real_type x ) const { return m_pSpline->D(x); }

    //! Second derivative
    real_type
    DD( real_type x ) const { return m_pSpline->DD(x); }

    //! Third derivative
    real_type
    DDD( real_type x ) const { return m_pSpline->DDD(x); }

    //! 4th derivative
    real_type
    DDDD( real_type x ) const { return m_pSpline->DDDD(x); }

    //! 5th derivative
    real_type
    DDDDD( real_type x ) const { return m_pSpline->DDDDD(x); }

    //! Some aliases
    real_type eval( real_type x ) const { return (*m_pSpline)(x); }
    real_type eval_D( real_type x ) const { return m_pSpline->D(x); }
    real_type eval_DD( real_type x ) const { return m_pSpline->DD(x); }
    real_type eval_DDD( real_type x ) const { return m_pSpline->DDD(x); }
    real_type eval_DDDD( real_type x ) const { return m_pSpline->DDDD(x); }
    real_type eval_DDDDD( real_type x ) const { return m_pSpline->DDDDD(x); }

    ///////////////////////////////////////////////////////////////////////////
    //! Evaluate spline value knowing interval
    real_type
    id_eval( integer ni, real_type x ) const { return m_pSpline->id_eval(ni,x); }

    //! First derivative
    real_type
    id_D( integer ni, real_type x ) const { return m_pSpline->id_D(ni,x); }

    //! Second derivative
    real_type
    id_DD( integer ni, real_type x ) const { return m_pSpline->id_DD(ni,x); }

    //! Third derivative
    real_type
    id_DDD( integer ni, real_type x ) const { return m_pSpline->id_DDD(ni,x); }

    //! 4th derivative
    real_type
    id_DDDD( integer ni, real_type x ) const { return m_pSpline->id_DDDD(ni,x); }

    //! 5th derivative
    real_type
    id_DDDDD( integer ni, real_type x ) const { return m_pSpline->id_DDDDD(ni,x); }

    //! get the piecewise polinomials of the spline
    integer // order
    coeffs(
      real_type cfs[],
      real_type nodes[],
      bool      transpose = false
    ) const {
      return m_pSpline->coeffs( cfs, nodes, transpose );
    }

    integer order() const { return m_pSpline->order(); }

    //! Print spline coefficients
    void
    writeToStream( ostream_type & s ) const
    { return m_pSpline->writeToStream( s ); }

    //! Return spline typename
    char const *
    type_name() const
    { return m_pSpline->type_name(); }

    //! Return spline type (as number)
    unsigned
    type() const
    { return m_pSpline->type(); }

    void
    info( ostream_type & s ) const
    { m_pSpline->info( s ); }
  };


  /*\
   |   ____        _ _          __     __
   |  / ___| _ __ | (_)_ __   __\ \   / /__  ___
   |  \___ \| '_ \| | | '_ \ / _ \ \ / / _ \/ __|
   |   ___) | |_) | | | | | |  __/\ V /  __/ (__
   |  |____/| .__/|_|_|_| |_|\___| \_/ \___|\___|
   |        |_|
  \*/

  //! Splines Management Class
  class SplineVec {

    SplineVec( SplineVec const & ) = delete;
    SplineVec const & operator = ( SplineVec const & ) = delete;

  protected:

    string const m_name;

    Utils::Malloc<real_type>  m_baseValue;
    Utils::Malloc<real_type*> m_basePointer;

    integer m_dim;
    integer m_npts;
    bool    m_curve_is_closed;
    bool    m_curve_can_extend;

    real_type *  m_X;
    real_type ** m_Y;
    real_type ** m_Yp;

    mutable Utils::BinarySearch<integer> m_bs;

    void initLastInterval();

    ///////////////////////////////////////////////////////////////////////////
    void
    allocate( integer dim, integer npts );

    void
    computeChords();

  public:

    //! spline constructor
    SplineVec( string const & name = "SplineVec" );

    //! spline destructor
    virtual
    ~SplineVec();

    integer search( real_type & x ) const;

    bool is_closed() const { return m_curve_is_closed; }
    void make_closed()     { m_curve_is_closed = true; }
    void make_open()       { m_curve_is_closed = false; }

    bool can_extend() const { return m_curve_can_extend; }
    void make_unbounded()   { m_curve_can_extend = true; }
    void make_buonded()     { m_curve_can_extend = false; }

    string const & name() const { return m_name; }

    //! return the number of support points of the splines
    integer
    numPoints(void) const
    { return m_npts; }

    //! return the number splines in the spline set
    integer
    dimension(void) const
    { return m_dim; }

    //! return the vector of values of x-nodes
    real_type const *
    xNodes() const
    { return m_X; }

    //! return the npt-th node of the spline (x component).
    real_type
    xNode( integer npt ) const
    { return m_X[size_t(npt)]; }

    //! return the npt-th node of the spline (y component).
    real_type
    yNode( integer npt, integer j ) const
    { return m_Y[size_t(j)][size_t(npt)]; }

    //! return x-minumum spline value
    real_type
    xMin() const
    { return m_X[0]; }

    //! return x-maximum spline value
    real_type
    xMax() const
    { return m_X[size_t(m_npts-1)]; }

    //! Evaluate spline value
    real_type
    operator () ( real_type x, integer i ) const;

    real_type
    eval( real_type x, integer i ) const
    { return operator() (x,i); }

    //! First derivative
    real_type
    D( real_type x, integer i ) const;

    real_type
    eval_D( real_type x, integer i ) const
    { return this->D(x,i); }

    //! Second derivative
    real_type
    DD( real_type x, integer i ) const;

    real_type
    eval_DD( real_type x, integer i ) const
    { return this->DD(x,i); }

    //! Third derivative
    real_type
    DDD( real_type x, integer i ) const;

    real_type
    eval_DDD( real_type x, integer i ) const
    { return this->DDD(x,i); }

    //! 4th derivative
    real_type
    DDDD( real_type x, integer i ) const;

    real_type
    eval_DDDD( real_type x, integer i ) const
    { return this->DDDD(x,i); }

    //! 5th derivative
    real_type
    DDDDD( real_type x, integer i ) const;

    real_type
    eval_DDDDD( real_type x, integer i ) const
    { return this->DDDDD(x,i); }

    //! Evaluate all the splines at `x`
    void
    eval(
      real_type x,
      real_type vals[],
      integer   inc
    ) const;

    //! Evaluate the fist derivative of all the splines at `x`
    void
    eval_D(
      real_type x,
      real_type vals[],
      integer   inc
    ) const;

    //! Evaluate the second derivative of all the splines at `x`
    void
    eval_DD(
      real_type x,
      real_type vals[],
      integer   inc
    ) const;

    //! Evaluate the third derivative of all the splines at `x`
    void
    eval_DDD(
      real_type x,
      real_type vals[],
      integer   inc
    ) const;

    //! Evaluate the 4th derivative of all the splines at `x`
    void
    eval_DDDD(
      real_type x,
      real_type vals[],
      integer   inc
    ) const;

    //! Evaluate the 5th derivative of all the splines at `x`
    void
    eval_DDDDD(
      real_type x,
      real_type vals[],
      integer   inc
    ) const;

    //! Evaluate all the splines at `x`
    void
    eval( real_type x, vector<real_type> & vals ) const;

    //! Evaluate the fist derivative of all the splines at `x`
    void
    eval_D( real_type x, vector<real_type> & vals ) const;

    //! Evaluate the second derivative of all the splines at `x`
    void
    eval_DD( real_type x, vector<real_type> & vals ) const;

    //! Evaluate the third derivative of all the splines at `x`
    void
    eval_DDD( real_type x, vector<real_type> & vals ) const;

    //! Evaluate the 4th derivative of all the splines at `x`
    void
    eval_DDDD( real_type x, vector<real_type> & vals ) const;

    //! Evaluate the 5th derivative of all the splines at `x`
    void
    eval_DDDDD( real_type x, vector<real_type> & vals ) const;

    /*!
     | Evaluate at `x` and fill a GenericContainer
    \*/
    void
    eval( real_type x, GenericContainer & vals ) const
    { eval( x, vals.set_vec_real() ); }

    void
    eval_D( real_type x, GenericContainer & vals ) const
    { eval_D( x, vals.set_vec_real() ); }

    void
    eval_DD( real_type x, GenericContainer & vals ) const
    { eval_DD( x, vals.set_vec_real() ); }

    void
    eval_DDD( real_type x, GenericContainer & vals ) const
    { eval_DDD( x, vals.set_vec_real() ); }

    void
    eval_DDDD( real_type x, GenericContainer & vals ) const
    { eval_DDDD( x, vals.set_vec_real() ); }

    void
    eval_DDDDD( real_type x, GenericContainer & vals ) const
    { eval_DDDDD( x, vals.set_vec_real() ); }

    /*!
     | Evaluate at `x` and fill a GenericContainer
    \*/
    void
    eval( vec_real_type const & x, GenericContainer & vals ) const;

    void
    eval_D( vec_real_type const & x, GenericContainer & vals ) const;

    void
    eval_DD( vec_real_type const & x, GenericContainer & vals ) const;

    void
    eval_DDD( vec_real_type const & x, GenericContainer & vals ) const;

    void
    eval_DDDD( vec_real_type const & x, GenericContainer & vals ) const;

    void
    eval_DDDDD( vec_real_type const & x, GenericContainer & vals ) const;

    void
    setup(
      integer           dim,
      integer           npts,
      real_type const * Y[]
    );

    void
    setup(
      integer         dim,
      integer         npts,
      real_type const Y[],
      integer         ldY
    );

    void
    setKnots( real_type const X[] );

    void
    setKnotsChordLength();

    void
    setKnotsCentripetal();

    void
    setKnotsFoley();

    void
    CatmullRom();

    real_type
    curvature( real_type x ) const;

    real_type
    curvature_D( real_type x ) const;

    virtual
    void
    setup( GenericContainer const & gc );

    void
    build( GenericContainer const & gc )
    { setup(gc); }

    //! Return spline type (as number)
    virtual
    unsigned
    type() const
    { return SPLINE_VEC_TYPE; }

    void
    info( ostream_type & s ) const;

    void
    dump_table( ostream_type & s, integer num_points ) const;

  };

  /*\
   |   ____        _ _            ____       _
   |  / ___| _ __ | (_)_ __   ___/ ___|  ___| |_
   |  \___ \| '_ \| | | '_ \ / _ \___ \ / _ \ __|
   |   ___) | |_) | | | | | |  __/___) |  __/ |_
   |  |____/| .__/|_|_|_| |_|\___|____/ \___|\__|
   |        |_|
  \*/

  //! Splines Management Class
  class SplineSet {

    SplineSet( SplineSet const & ) = delete;
    SplineSet const & operator = ( SplineSet const & ) = delete;

    class BinarySearch {
    public:
      typedef std::pair<std::string,integer> DATA_TYPE;
    private:
      mutable std::vector<DATA_TYPE> data;
    public:
      BinarySearch() { data.clear(); data.reserve(256); }
      ~BinarySearch() { data.clear(); }

      void clear() { data.clear(); data.reserve(256); }

      integer n_elem() const { return integer(data.size()); }

      DATA_TYPE const &
      get_elem( integer i ) const { return data[size_t(i)]; }

      integer search( std::string const & id ) const;
      void    insert( std::string const & id, integer position );
    };

  protected:

    string const m_name;

    Utils::Malloc<real_type>  m_baseValue;
    Utils::Malloc<real_type*> m_basePointer;
    Utils::Malloc<void*>      m_baseSplines;
    Utils::Malloc<int>        m_baseInt;

    integer m_npts;
    integer m_nspl;

    real_type *  m_X;
    real_type ** m_Y;
    real_type ** m_Yp;
    real_type ** m_Ypp;
    real_type *  m_Ymin;
    real_type *  m_Ymax;

    Spline ** m_splines;
    int     * m_is_monotone;

    BinarySearch m_header_to_position;

  private:

    /*!
     | find `x` value such that the monotone spline
     | `(spline[spl])(x)` intersect the value `zeta`
    \*/
    Spline const *
    intersect( integer spl, real_type zeta, real_type & x ) const;

  public:

    //! spline constructor
    SplineSet( string const & name = "SplineSet" );

    //! spline destructor
    virtual
    ~SplineSet();

    string const & name() const { return m_name; }

    string const &
    header( integer i ) const
    { return m_splines[i]->name(); }

    // vectorial values
    //! fill a vector of strings with the names of the splines
    void get_headers( std::vector<std::string> & names ) const;

    string name_list() const;

    // return +1 = strict monotone, 0 weak monotone, -1 non monotone
    int
    isMonotone( integer i ) const
    { return m_is_monotone[i]; }

    //! return the number of support points of the splines
    integer numPoints(void) const { return m_npts; }

    //! return the number splines in the spline set
    integer numSplines(void) const { return m_nspl; }

    //! return the column with header(i) == hdr, return -1 if not found
    integer getPosition( char const * hdr ) const;

    //! return the vector of values of x-nodes
    real_type const * xNodes() const { return m_X; }

    //! return the vector of values of x-nodes
    real_type const * yNodes( integer i ) const;

    //! return the i-th node of the spline (x component).
    real_type
    xNode( integer npt ) const
    { return m_X[npt]; }

    //! return the i-th node of the spline (y component).
    real_type
    yNode( integer npt, integer spl ) const
    { return m_Y[spl][npt]; }

    //! return x-minumum spline value
    real_type xMin() const { return m_X[0]; }

    //! return x-maximum spline value
    real_type xMax() const { return m_X[m_npts-1]; }

    //! return y-minumum spline value
    real_type yMin( integer spl ) const { return m_Ymin[size_t(spl)]; }

    //! return y-maximum spline value
    real_type yMax( integer spl ) const { return m_Ymax[size_t(spl)]; }

    //! return y-minumum spline value
    real_type
    yMin( char const spl[] ) const {
      integer idx = this->getPosition(spl);
      return m_Ymin[idx];
    }

    //! return y-maximum spline value
    real_type
    yMax( char const spl[] ) const {
      integer idx = this->getPosition(spl);
      return m_Ymax[idx];
    }

    //! Return pointer to the `i`-th spline
    Spline * getSpline( integer i ) const;

    //! Return pointer to the `i`-th spline
    Spline *
    getSpline( char const * hdr ) const {
      integer idx = this->getPosition(hdr);
      return m_splines[idx];
    }

    //! Evaluate spline value
    real_type
    operator () ( real_type x, integer spl ) const {
      Spline const * S = this->getSpline(spl);
      return (*S)(x);
    }

    real_type
    eval( real_type x, integer spl ) const {
      Spline const * S = this->getSpline(spl);
      return (*S)(x);
    }

    //! First derivative
    real_type
    D( real_type x, integer spl ) const {
      Spline const * S = this->getSpline(spl);
      return S->D(x);
    }

    real_type
    eval_D( real_type x, integer spl ) const {
      Spline const * S = this->getSpline(spl);
      return S->D(x);
    }

    //! Second derivative
    real_type
    DD( real_type x, integer spl ) const {
      Spline const * S = this->getSpline(spl);
      return S->DD(x);
    }

    real_type
    eval_DD( real_type x, integer spl ) const {
      Spline const * S = this->getSpline(spl);
      return S->DD(x);
    }

    //! Third derivative
    real_type
    DDD( real_type x, integer spl ) const {
      Spline const * S = this->getSpline(spl);
      return S->DDD(x);
    }

    real_type
    eval_DDD( real_type x, integer spl ) const {
      Spline const * S = this->getSpline(spl);
      return S->DDD(x);
    }

    //! 4th derivative
    real_type
    DDDD( real_type x, integer spl ) const {
      Spline const * S = this->getSpline(spl);
      return S->DDDD(x);
    }

    real_type
    eval_DDDD( real_type x, integer spl ) const {
      Spline const * S = this->getSpline(spl);
      return S->DDDD(x);
    }

    //! 5th derivative
    real_type
    DDDDD( real_type x, integer spl ) const {
      Spline const * S = this->getSpline(spl);
      return S->DDDDD(x);
    }

    real_type
    eval_DDDDD( real_type x, integer spl ) const {
      Spline const * S = this->getSpline(spl);
      return S->DDDDD(x);
    }

    //! Evaluate spline value
    real_type
    eval( real_type x, char const * name ) const {
      Spline const * S = this->getSpline(name);
      return (*S)(x);
    }

    //! First derivative
    real_type
    eval_D( real_type x, char const * name ) const {
      Spline const * S = this->getSpline(name);
      return S->D(x);
    }

    //! Second derivative
    real_type
    eval_DD( real_type x, char const * name ) const {
      Spline const * S = this->getSpline(name);
      return S->DD(x);
    }

    //! Third derivative
    real_type
    eval_DDD( real_type x, char const * name ) const {
      Spline const * S = this->getSpline(name);
      return S->DDD(x);
    }

    //! 4th derivative
    real_type
    eval_DDDD( real_type x, char const * name ) const {
      Spline const * S = this->getSpline(name);
      return S->DDDD(x);
    }

    //! 5th derivative
    real_type
    eval_DDDDD( real_type x, char const * name ) const {
      Spline const * S = this->getSpline(name);
      return S->DDDDD(x);
    }

    //! Evaluate all the splines at `x`
    void
    eval(
      real_type           x,
      vector<real_type> & vals
    ) const;

    //! Evaluate all the splines at `x`
    void
    eval(
      real_type x,
      real_type vals[],
      integer   incy = 1
    ) const;

    //! Evaluate the fist derivative of all the splines at `x`
    void
    eval_D(
      real_type           x,
      vector<real_type> & vals
    ) const;

    //! Evaluate the fist derivative of all the splines at `x`
    void
    eval_D(
      real_type x,
      real_type vals[],
      integer   incy = 1
    ) const;

    //! Evaluate the second derivative of all the splines at `x`
    void
    eval_DD(
      real_type           x,
      vector<real_type> & vals
    ) const;

    //! Evaluate the second derivative of all the splines at `x`
    void
    eval_DD(
      real_type x,
      real_type vals[],
      integer   incy = 1
    ) const;

    //! Evaluate the third derivative of all the splines at `x`
    void
    eval_DDD(
      real_type           x,
      vector<real_type> & vals
    ) const;

    //! Evaluate the third derivative of all the splines at `x`
    void
    eval_DDD(
      real_type x,
      real_type vals[],
      integer   incy = 1
    ) const;

    //! Evaluate the 4th derivative of all the splines at `x`
    void
    eval_DDDD(
      real_type           x,
      vector<real_type> & vals
    ) const;

    //! Evaluate the 4th derivative of all the splines at `x`
    void
    eval_DDDD(
      real_type x,
      real_type vals[],
      integer   incy = 1
    ) const;

    //! Evaluate the 5th derivative of all the splines at `x`
    void
    eval_DDDDD(
      real_type           x,
      vector<real_type> & vals
    ) const;

    //! Evaluate the 5th derivative of all the splines at `x`
    void
    eval_DDDDD(
      real_type x,
      real_type vals[],
      integer   incy = 1
    ) const;

    // change independent variable
    //! Evaluate all the splines at `zeta` using spline[spl] as independent
    void
    eval2(
      integer             spl,
      real_type           zeta,
      vector<real_type> & vals
    ) const;

    //! Evaluate all the splines at `zeta` using spline[spl] as independent
    void
    eval2(
      integer   spl,
      real_type zeta,
      real_type vals[],
      integer   incy = 1
    ) const;

    /*!
     | Evaluate the fist derivative of all the splines
     | at `zeta` using spline[spl] as independent
    \*/
    void
    eval2_D(
      integer             spl,
      real_type           zeta,
      vector<real_type> & vals
    ) const;

    /*!
     | Evaluate the fist derivative of all the splines
     | at `zeta` using spline[spl] as independent
    \*/
    void
    eval2_D(
      integer   spl,
      real_type zeta,
      real_type vals[],
      integer   incy = 1
    ) const;

    /*!
     | Evaluate the second derivative of all the splines
     | at `zeta` using spline[spl] as independent
    \*/
    void
    eval2_DD(
      integer             spl,
      real_type           zeta,
      vector<real_type> & vals
    ) const;

    /*!
     | Evaluate the second derivative of all the splines
     | at `zeta` using spline[spl] as independent
    \*/
    void
    eval2_DD(
      integer   spl,
      real_type zeta,
      real_type vals[],
      integer   incy = 1
    ) const;

    /*!
     | Evaluate the third derivative of all the splines
     | at `zeta` using spline[spl] as independent
    \*/
    void
    eval2_DDD(
      integer             spl,
      real_type           zeta,
      vector<real_type> & vals
    ) const;

    /*!
     | Evaluate the 3rd derivative of all the splines
     | at `zeta` using spline[spl] as independent
    \*/
    void
    eval2_DDD(
      integer   spl,
      real_type zeta,
      real_type vals[],
      integer   incy = 1
    ) const;

    /*!
     | Evaluate the spline `name` at `zeta` using
     | spline `indep` as independent
    \*/
    real_type
    eval2(
      real_type  zeta,
      char const indep[],
      char const name[]
    ) const;

    /*!
     | Evaluate first derivative of the spline `name`
     | at `zeta` using spline `indep` as independent
    \*/
    real_type
    eval2_D(
      real_type  zeta,
      char const indep[],
      char const name[]
    ) const;

    /*!
     | Evaluate second derivative of the spline `name`
     | at `zeta` using spline `indep` as independent
    \*/
    real_type
    eval2_DD(
      real_type  zeta,
      char const indep[],
      char const name[]
    ) const;

    /*!
     | Evaluate third derivative of the spline `name`
     | at `zeta` using spline `indep` as independent
    \*/
    real_type
    eval2_DDD(
      real_type  zeta,
      char const indep[],
      char const name[]
    ) const;

    /*!
     | Evaluate the spline `spl` at `zeta` using
     | spline `indep` as independent
    \*/
    real_type
    eval2(
      real_type zeta,
      integer   indep,
      integer   spl
    ) const;

    //! Evaluate first derivative of the spline `spl` at `zeta` using spline `indep` as independent
    real_type
    eval2_D(
      real_type zeta,
      integer   indep,
      integer   spl
    ) const;

    /*!
     | Evaluate second derivative of the spline `spl`
     | at `zeta` using spline `indep` as independent
    \*/
    real_type
    eval2_DD(
      real_type zeta,
      integer   indep,
      integer   spl
    ) const;

    /*!
     | Evaluate third derivative of the spline `spl`
     | at `zeta` using spline `indep` as independent
    \*/
    real_type
    eval2_DDD(
      real_type zeta,
      integer   indep,
      integer   spl
    ) const;

    /*!
     | Evaluate all the splines at `x`
     | and fill a map of values in a GenericContainer
    \*/
    void
    eval( real_type x, GenericContainer & vals ) const;

    /*!
     | Evaluate all the splines at `x` values contained
     | in vec and fill a map of vector in a GenericContainer
    \*/
    void
    eval( vec_real_type const & vec, GenericContainer & vals ) const;

    /*!
     | Evaluate all the splines at `x` and fill a map of
     | values in a GenericContainer with keys in `columns`
    \*/
    void
    eval(
      real_type               x,
      vec_string_type const & columns,
      GenericContainer      & vals
    ) const;

    /*!
     | Evaluate all the splines at `x` values contained
     | in vec and fill a map of vector in a GenericContainer
     | with keys in `columns`
    \*/
    void
    eval(
      vec_real_type   const & vec,
      vec_string_type const & columns,
      GenericContainer      & vals
    ) const;

    /*!
     | Evaluate all the splines at `zeta` and fill
     | a map of values in a GenericContainer and `indep`
     | as independent spline
    \*/
    void
    eval2(
      real_type          zeta,
      integer            indep,
      GenericContainer & vals
    ) const;

    /*!
     | Evaluate all the splines at `zeta` values contained in
     | vec and fill a map of vector in a GenericContainer and
     | `indep` as independent spline
    \*/
    void
    eval2(
      vec_real_type const & zetas,
      integer               indep,
      GenericContainer    & vals
    ) const;

    /*!
     | Evaluate all the splines at `zeta` and fill a map
     | of values in a GenericContainer with keys in `columns`
     | and `indep` as independent spline
    \*/
    void
    eval2(
      real_type               zeta,
      integer                 indep,
      vec_string_type const & columns,
      GenericContainer      & vals
    ) const;

    /*!
     | Evaluate all the splines at `zeta` values contained
     | in vec and fill a map of vector in a GenericContainer
     | with keys in `columns` and `indep` as independent spline
    \*/
    void
    eval2(
      vec_real_type   const & zetas,
      integer                 indep,
      vec_string_type const & columns,
      GenericContainer      & vals
    ) const;

    /*!
     | Evaluate all the splines at `zeta` and fill a map
     | of values in a GenericContainer and `indep` as independent spline
    \*/
    void
    eval2(
      real_type          zeta,
      char const       * indep,
      GenericContainer & vals
    ) const {
      this->eval2( zeta, this->getPosition(indep), vals );
    }

    /*!
     | Evaluate all the splines at `zeta` values contained
     | in vec and fill a map of vector in a GenericContainer
     | and `indep` as independent spline
    \*/
    void
    eval2(
      vec_real_type const & zetas,
      char const          * indep,
      GenericContainer    & vals
    ) const {
      this->eval2( zetas, this->getPosition(indep), vals );
    }

    /*!
     | Evaluate all the splines at `zeta` and fill a map of
     | values in a GenericContainer with keys in `columns`
     | and `indep` as independent spline
    \*/
    void
    eval2(
      real_type               zeta,
      char const            * indep,
      vec_string_type const & columns,
      GenericContainer      & vals
    ) const {
      this->eval2( zeta, this->getPosition(indep), columns, vals );
    }

    /*!
     | Evaluate all the splines at `zeta` values contained
     | in vec and fill a map of vector in a GenericContainer
     | with keys in `columns` and `indep` as independent spline
    \*/
    void
    eval2(
      vec_real_type const   & zetas,
      char const            * indep,
      vec_string_type const & columns,
      GenericContainer      & vals
    ) const {
      this->eval2( zetas, this->getPosition(indep), columns, vals );
    }

    /*!
     | Evaluate all the splines at `x` and fill a map of
     | values in a GenericContainer
    \*/
    void
    eval_D(
      real_type          x,
      GenericContainer & vals
    ) const;

    /*!
     | Evaluate all the splines at `x` values contained
     | in vec and fill a map of vector in a GenericContainer
    \*/
    void
    eval_D(
      vec_real_type const & vec,
      GenericContainer    & vals
    ) const;

    /*!
     | Evaluate all the splines at `x` and fill a map of
     | values in a GenericContainer with keys in `columns`
    \*/
    void
    eval_D(
      real_type               x,
      vec_string_type const & columns,
      GenericContainer      & vals
    ) const;

    /*!
     | Evaluate all the splines at `x` values contained
     | in vec and fill a map of vector in a GenericContainer
     | with keys in `columns`
    \*/
    void
    eval_D(
      vec_real_type const   & vec,
      vec_string_type const & columns,
      GenericContainer      & vals
    ) const;

    /*!
     | Evaluate all the splines at `zeta` and fill a map
     | of values in a GenericContainer and `indep` as independent spline
    \*/
    void
    eval2_D(
      real_type          zeta,
      integer            indep,
      GenericContainer & vals
    ) const;

    /*!
     | Evaluate all the splines at `zeta` values contained
     | in vec and fill a map of vector in a GenericContainer
     | and `indep` as independent spline
    \*/
    void
    eval2_D(
      vec_real_type const & zetas,
      integer               indep,
      GenericContainer    & vals
    ) const;

    /*!
     | Evaluate all the splines at `zeta` and fill a map
     | of values in a GenericContainer with keys in `columns`
     | and `indep` as independent spline
    \*/
    void
    eval2_D(
      real_type               zeta,
      integer                 indep,
      vec_string_type const & columns,
      GenericContainer      & vals
    ) const;

    /*!
     | Evaluate all the splines at `zeta` values contained
     | in vec and fill a map of vector in a GenericContainer
     | with keys in `columns` and `indep` as independent spline
    \*/
    void
    eval2_D(
      vec_real_type   const & zetas,
      integer                 indep,
      vec_string_type const & columns,
      GenericContainer      & vals
    ) const;

    /*!
     | Evaluate all the splines at `zeta` and fill a map
     | of values in a GenericContainer and `indep` as independent spline
    \*/
    void
    eval2_D(
      real_type          zeta,
      char const       * indep,
      GenericContainer & vals
    ) const {
      this->eval2_D( zeta, this->getPosition(indep), vals );
    }

    //! Evaluate all the splines at `zeta` values contained in vec and fill a map of vector in a GenericContainer and `indep` as independent spline
    void
    eval2_D(
      vec_real_type const & zetas,
      char const          * indep,
      GenericContainer    & vals
    ) const {
      this->eval2_D( zetas, this->getPosition(indep), vals );
    }

    /*!
     | Evaluate all the splines at `zeta` and fill a map of
     | values in a GenericContainer with keys in `columns`
     | and `indep` as independent spline
    \*/
    void
    eval2_D(
      real_type               zeta,
      char const            * indep,
      vec_string_type const & columns,
      GenericContainer      & vals
    ) const {
      this->eval2_D( zeta, this->getPosition(indep), columns, vals );
    }

    /*!
     | Evaluate all the splines at `zeta` values contained
     | in vec and fill a map of vector in a GenericContainer
     | with keys in `columns` and `indep` as independent spline
    \*/
    void
    eval2_D(
      vec_real_type const   & zetas,
      char const            * indep,
      vec_string_type const & columns,
      GenericContainer      & vals
    ) const {
      this->eval2_D( zetas, this->getPosition(indep), columns, vals );
    }

    /*!
     | Evaluate all the splines at `x` and fill a map of
     | values in a GenericContainer
    \*/
    void
    eval_DD(
      real_type          x,
      GenericContainer & vals
    ) const;

    /*!
     | Evaluate all the splines at `x` values contained in vec
     | and fill a map of vector in a GenericContainer
    \*/
    void
    eval_DD(
      vec_real_type const & vec,
      GenericContainer    & vals
    ) const;

    /*!
     | Evaluate all the splines at `x` and fill a map of
     | values in a GenericContainer with keys in `columns`
    \*/
    void
    eval_DD(
      real_type               x,
      vec_string_type const & columns,
      GenericContainer      & vals
    ) const;

    /*!
     | Evaluate all the splines at `x` values contained in vec and
     | fill a map of vector in a GenericContainer with keys in `columns`
    \*/
    void
    eval_DD(
      vec_real_type   const & vec,
      vec_string_type const & columns,
      GenericContainer      & vals
    ) const;

    /*!
     | Evaluate all the splines at `zeta` and fill a map of
     | values in a GenericContainer and `indep` as independent spline
    \*/
    void
    eval2_DD(
      real_type          zeta,
      integer            indep,
      GenericContainer & vals
    ) const;

    /*!
     | Evaluate all the splines at `zeta` values contained in vec
     | and fill a map of vector in a GenericContainer
     | and `indep` as independent spline
    \*/
    void
    eval2_DD(
      vec_real_type const & zetas,
      integer               indep,
      GenericContainer    & vals
    ) const;

    /*!
     | Evaluate all the splines at `zeta` and fill a map
     | of values in a GenericContainer with keys in `columns`
     | and `indep` as independent spline
    \*/
    void
    eval2_DD(
      real_type               zeta,
      integer                 indep,
      vec_string_type const & columns,
      GenericContainer      & vals
    ) const;

    /*!
     | Evaluate all the splines at `zeta` values contained
     | in vec and fill a map of vector in a GenericContainer
     | with keys in `columns` and `indep` as independent spline
    \*/
    void
    eval2_DD(
      vec_real_type   const & zetas,
      integer                 indep,
      vec_string_type const & columns,
      GenericContainer      & vals
    ) const;

    /*!
     | Evaluate all the splines at `zeta` and fill a map of
     | values in a GenericContainer and `indep` as independent spline
    \*/
    void
    eval2_DD(
      real_type          zeta,
      char const       * indep,
      GenericContainer & vals
    ) const {
      this->eval2_DD( zeta, this->getPosition(indep), vals );
    }

    /*!
     | Evaluate all the splines at `zeta` values contained in vec
     | and fill a map of vector in a GenericContainer
     | and `indep` as independent spline
    \*/
    void
    eval2_DD(
      vec_real_type const & zetas,
      char const          * indep,
      GenericContainer    & vals
    ) const {
      this->eval2_DD( zetas, this->getPosition(indep), vals );
    }

    /*!
     | Evaluate all the splines at `zeta` and fill a map
     | of values in a GenericContainer with keys in `columns`
     | and `indep` as independent spline
    \*/
    void
    eval2_DD(
      real_type               zeta,
      char const            * indep,
      vec_string_type const & columns,
      GenericContainer      & vals
    ) const {
      this->eval2_DD( zeta, this->getPosition(indep), columns, vals );
    }

    /*!
     | Evaluate all the splines at `zeta` values contained
     | in vec and fill a map of vector in a GenericContainer
     | with keys in `columns` and `indep` as independent spline
    \*/
    void
    eval2_DD(
      vec_real_type   const & zetas,
      char            const * indep,
      vec_string_type const & columns,
      GenericContainer      & vals
    ) const {
      this->eval2_DD( zetas, this->getPosition(indep), columns, vals );
    }

    /*!
     | Evaluate all the splines at `x` and fill a map
     | of values in a GenericContainer
    \*/
    void
    eval_DDD(
      real_type          x,
      GenericContainer & vals
    ) const;

    /*!
     | Evaluate all the splines at `x` values contained
     | in vec and fill a map of vector in a GenericContainer
    \*/
    void
    eval_DDD(
      vec_real_type const & vec,
      GenericContainer    & vals
    ) const;

    /*!
     | Evaluate all the splines at `x` and fill a map of
     | values in a GenericContainer with keys in `columns`
    \*/
    void
    eval_DDD(
      real_type               x,
      vec_string_type const & columns,
      GenericContainer      & vals
    ) const;

    /*!
     | Evaluate all the splines at `x` values contained in vec
     | and fill a map of vector in a GenericContainer with keys in `columns`
    \*/
    void
    eval_DDD(
      vec_real_type   const & vec,
      vec_string_type const & columns,
      GenericContainer      & vals
    ) const;

    /*!
     | Evaluate all the splines at `zeta` and fill a map of values
     | in a GenericContainer and `indep` as independent spline
    \*/
    void
    eval2_DDD(
      real_type          zeta,
      integer            indep,
      GenericContainer & vals
    ) const;

    /*!
     | Evaluate all the splines at `zeta` values contained in vec
     | and fill a map of vector in a GenericContainer
     | and `indep` as independent spline
    \*/
    void
    eval2_DDD(
      vec_real_type const & zetas,
      integer               indep,
      GenericContainer    & vals
    ) const;

    /*!
     | Evaluate all the splines at `zeta` and fill a map of values
     | in a GenericContainer with keys in `columns`
     | and `indep` as independent spline
    \*/
    void
    eval2_DDD(
      real_type               zeta,
      integer                 indep,
      vec_string_type const & columns,
      GenericContainer      & vals
    ) const;

    /*!
     | Evaluate all the splines at `zeta` values contained
     | in vec and fill a map of vector in a GenericContainer
     | with keys in `columns` and `indep` as independent spline
    \*/
    void
    eval2_DDD(
      vec_real_type   const & zetas,
      integer                 indep,
      vec_string_type const & columns,
      GenericContainer      & vals
    ) const;

    /*!
     | Evaluate all the splines at `zeta` and fill a map
     | of values in a GenericContainer and `indep` as independent spline
    \*/
    void
    eval2_DDD(
      real_type          zeta,
      char const       * indep,
      GenericContainer & vals
    ) const {
      this->eval2_DDD( zeta, this->getPosition(indep), vals );
    }

    /*!
     | Evaluate all the splines at `zeta` values contained
     | in vec and fill a map of vector in a GenericContainer
     | and `indep` as independent spline
    \*/
    void
    eval2_DDD(
      vec_real_type const & zetas,
      char          const * indep,
      GenericContainer    & vals
    ) const {
      this->eval2_DDD( zetas, this->getPosition(indep), vals );
    }

    /*!
     | Evaluate all the splines at `zeta` and fill a map of
     | values in a GenericContainer with keys in `columns`
     | and `indep` as independent spline
    \*/
    void
    eval2_DDD(
      real_type               zeta,
      char            const * indep,
      vec_string_type const & columns,
      GenericContainer      & vals
    ) const {
      this->eval2_DDD( zeta, this->getPosition(indep), columns, vals );
    }

    /*!
     | Evaluate all the splines at `zeta` values contained
     | in vec and fill a map of vector in a GenericContainer
     | with keys in `columns` and `indep` as independent spline
    \*/
    void
    eval2_DDD(
      vec_real_type   const & zetas,
      char            const * indep,
      vec_string_type const & columns,
      GenericContainer      & vals
    ) const {
      this->eval2_DDD( zetas, this->getPosition(indep), columns, vals );
    }

    ///////////////////////////////////////////////////////////////////////////
    /*! Build a set of splines
     | \param nspl       the number of splines
     | \param npts       the number of points of each splines
     | \param headers    the names of the splines
     | \param stype      the type of each spline
     | \param X          pointer to X independent values
     | \param Y          vector of `nspl` pointers to Y depentendent values.
    \*/

    void
    build(
      integer            nspl,
      integer            npts,
      char         const *headers[],
      SplineType1D const stype[],
      real_type    const X[],
      real_type    const *Y[],
      real_type    const *Yp[] = nullptr
    );

    void
    setup( GenericContainer const & gc );

    void
    build( GenericContainer const & gc )
    { this->setup(gc); }

    //! Return spline type (as number)
    unsigned
    type() const
    { return SPLINE_SET_TYPE; }

    void
    info( ostream_type & s ) const;

    void
    dump_table( ostream_type & s, integer num_points ) const;

  };

  /*\
   |   ____        _ _            ____              __
   |  / ___| _ __ | (_)_ __   ___/ ___| _   _ _ __ / _|
   |  \___ \| '_ \| | | '_ \ / _ \___ \| | | | '__| |_
   |   ___) | |_) | | | | | |  __/___) | |_| | |  |  _|
   |  |____/| .__/|_|_|_| |_|\___|____/ \__,_|_|  |_|
   |        |_|
  \*/
  //! Spline Management Class
  class SplineSurf {

    SplineSurf( SplineSurf const & ) = delete; // block copy constructor
    SplineSurf const & operator = ( SplineSurf const & ) = delete; // block copy method

    Utils::Malloc<real_type> m_mem;

  protected:

    string const m_name;
    bool         m_x_closed;
    bool         m_y_closed;
    bool         m_x_can_extend;
    bool         m_y_can_extend;

    integer      m_nx;
    integer      m_ny;

    real_type *  m_X;
    real_type *  m_Y;
    real_type *  m_Z;

    real_type m_Z_min, m_Z_max;

    mutable Utils::BinarySearch<integer> m_bs_x;
    mutable Utils::BinarySearch<integer> m_bs_y;

    integer search_x( real_type & x ) const;
    integer search_y( real_type & y ) const;

    void initLastInterval_x();
    void initLastInterval_y();

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

    virtual void makeSpline() UTILS_PURE_VIRTUAL;

  public:

    //! spline constructor
    SplineSurf( string const & name = "Spline" )
    : m_mem("SplineSurf")
    , m_name(name)
    , m_x_closed(false)
    , m_y_closed(false)
    , m_x_can_extend(true)
    , m_y_can_extend(true)
    , m_nx(0)
    , m_ny(0)
    , m_X(nullptr)
    , m_Y(nullptr)
    , m_Z(nullptr)
    , m_Z_min(0)
    , m_Z_max(0)
    {
      this->initLastInterval_x();
      this->initLastInterval_y();
    }

    //! spline destructor
    virtual
    ~SplineSurf();

    bool is_x_closed() const { return m_x_closed; }
    void make_x_closed()     { m_x_closed = true; }
    void make_x_opened()     { m_x_closed = false; }

    bool is_y_closed() const { return m_y_closed; }
    void make_y_closed()     { m_y_closed = true; }
    void make_y_opened()     { m_y_closed = false; }

    bool is_x_bounded() const { return m_x_can_extend; }
    void make_x_unbounded()   { m_x_can_extend = true; }
    void make_x_bounded()     { m_x_can_extend = false; }

    bool is_y_bounded() const { return m_y_can_extend; }
    void make_y_unbounded()   { m_y_can_extend = true; }
    void make_y_bounded()     { m_y_can_extend = false; }

    string const & name() const { return m_name; }

    //! Cancel the support points, empty the spline.
    void clear(void);

    //! return the number of support points of the spline along x direction
    integer numPointX(void) const { return m_nx; }

    //! return the number of support points of the spline along y direction
    integer numPointY(void) const { return m_ny; }

    //! return the i-th node of the spline (x component).
    real_type xNode( integer i ) const { return m_X[size_t(i)]; }

    //! return the i-th node of the spline (y component).
    real_type yNode( integer i ) const { return m_Y[size_t(i)]; }

    //! return the i-th node of the spline (y component).
    real_type
    zNode( integer i, integer j ) const
    { return m_Z[size_t(this->ipos_C(i,j))]; }

    //! return x-minumum spline value
    real_type xMin() const { return m_X[0]; }

    //! return x-maximum spline value
    real_type xMax() const { return m_X[m_nx-1]; }

    //! return y-minumum spline value
    real_type yMin() const { return m_Y[0]; }

    //! return y-maximum spline value
    real_type yMax() const { return m_Y[m_ny-1]; }

    //! return z-minumum spline value
    real_type zMin() const { return m_Z_min; }

    //! return z-maximum spline value
    real_type zMax() const { return m_Z_max; }

    ///////////////////////////////////////////////////////////////////////////
    /*! Build surface spline
     | \param x       vector of x-coordinates
     | \param incx    access elements as x[0], x[incx], x[2*incx],...
     | \param y       vector of y-coordinates
     | \param incy    access elements as y[0], y[incy], x[2*incy],...
     | \param z       matrix of z-values
     | \param ldZ     leading dimension of the matrix. Elements are stored
     |                by row Z(i,j) = z[i*ldZ+j] as C-matrix
     | \param nx      total number of points in direction x
     | \param ny      total number of points in direction y
     | \param fortran_storage if true elements are stored by column
     |                        i.e. Z(i,j) = z[i+j*ldZ] as Fortran-matrix
     | \param transposed      if true matrix Z is stored transposed
    \*/
    void
    build(
      real_type const x[], integer incx,
      real_type const y[], integer incy,
      real_type const z[], integer ldZ,
      integer nx, integer ny,
      bool fortran_storage = false,
      bool transposed      = false
    );

    /*! Build surface spline
     | \param x       vector of x-coordinates, nx = x.size()
     | \param y       vector of y-coordinates, ny = y.size()
     | \param z       matrix of z-values. Elements are stored
     |                by row Z(i,j) = z[i*ny+j] as C-matrix
     | \param fortran_storage if true elements are stored by column
     |                        i.e. Z(i,j) = z[i+j*nx] as Fortran-matrix
     | \param transposed      if true matrix Z is stored transposed
    \*/
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

    /*! Build surface spline
     | \param z               matrix of z-values. Elements are stored
     |                        by row Z(i,j) = z[i*ny+j] as C-matrix
     | \param ldZ             leading dimension of the matrix. Elements are stored
     |                        by row Z(i,j) = z[i*ldZ+j] as C-matrix
     | \param fortran_storage if true elements are stored by column
     |                        i.e. Z(i,j) = z[i+j*nx] as Fortran-matrix
     | \param transposed      if true matrix Z is stored transposed
    \*/
    void
    build(
      real_type const z[],
      integer         ldZ,
      integer         nx,
      integer         ny,
      bool fortran_storage = false,
      bool transposed      = false
    );

    /*! Build surface spline
     | \param z               matrix of z-values. Elements are stored
     |                        by row Z(i,j) = z[i*ny+j] as C-matrix.
     |                        ldZ leading dimension of the matrix is ny for C-storage
     |                        and nx for Fortran storage.
     | \param fortran_storage if true elements are stored by column
     |                        i.e. Z(i,j) = z[i+j*nx] as Fortran-matrix
     | \param transposed      if true matrix Z is stored transposed
    \*/
    void
    build(
      vector<real_type> const & z,
      integer                   nx,
      integer                   ny,
      bool fortran_storage = false,
      bool transposed      = false
    ) {
      if ( fortran_storage )
        this->build( &z.front(), nx, nx, ny, fortran_storage, transposed );
      else
        this->build( &z.front(), ny, nx, ny, fortran_storage, transposed );
    }

    void
    setup( GenericContainer const & gc );

    void
    build ( GenericContainer const & gc )
    { setup(gc); }

    //! Evaluate spline value
    virtual
    real_type
    operator () ( real_type x, real_type y ) const UTILS_PURE_VIRTUAL;

    //! First derivative
    virtual
    void
    D( real_type x, real_type y, real_type d[3] ) const UTILS_PURE_VIRTUAL;

    virtual
    real_type
    Dx( real_type x, real_type y ) const UTILS_PURE_VIRTUAL;

    virtual
    real_type
    Dy( real_type x, real_type y ) const UTILS_PURE_VIRTUAL;

    //! Second derivative
    virtual
    void
    DD( real_type x, real_type y, real_type dd[6] ) const UTILS_PURE_VIRTUAL;

    virtual
    real_type
    Dxx( real_type x, real_type y ) const UTILS_PURE_VIRTUAL;

    virtual
    real_type
    Dxy( real_type x, real_type y ) const UTILS_PURE_VIRTUAL;

    virtual
    real_type
    Dyy( real_type x, real_type y ) const UTILS_PURE_VIRTUAL;

    //! Evaluate spline value
    real_type
    eval( real_type x, real_type y ) const
    { return (*this)(x,y); }

    //! First derivative
    real_type
    eval_D_1( real_type x, real_type y ) const
    { return this->Dx(x,y); }

    real_type
    eval_D_2( real_type x, real_type y ) const
    { return this->Dy(x,y); }

    //! Second derivative
    real_type
    eval_D_1_1( real_type x, real_type y ) const
    { return this->Dxx(x,y); }

    real_type
    eval_D_1_2( real_type x, real_type y ) const
    { return this->Dxy(x,y); }

    real_type
    eval_D_2_2( real_type x, real_type y ) const
    { return this->Dyy(x,y); }

    //! Print spline coefficients
    virtual
    void
    writeToStream( ostream_type & s ) const UTILS_PURE_VIRTUAL;

    //! Return spline typename
    virtual
    char const *
    type_name() const UTILS_PURE_VIRTUAL;

    void
    info( ostream_type & s ) const;

  };

  /*\
   |   ____  _ _ _                       ____        _ _
   |  | __ )(_) (_)_ __   ___  __ _ _ __/ ___| _ __ | (_)_ __   ___
   |  |  _ \| | | | '_ \ / _ \/ _` | '__\___ \| '_ \| | | '_ \ / _ \
   |  | |_) | | | | | | |  __/ (_| | |   ___) | |_) | | | | | |  __/
   |  |____/|_|_|_|_| |_|\___|\__,_|_|  |____/| .__/|_|_|_| |_|\___|
   |                                          |_|
  \*/
  //! bilinear spline base class
  class BilinearSpline : public SplineSurf {
    virtual void makeSpline() UTILS_OVERRIDE {}
  public:

    using SplineSurf::m_nx;
    using SplineSurf::m_ny;

    using SplineSurf::m_X;
    using SplineSurf::m_Y;
    using SplineSurf::m_Z;

    //! spline constructor
    BilinearSpline( string const & name = "Spline" )
    : SplineSurf(name)
    {}

    virtual
    ~BilinearSpline() UTILS_OVERRIDE
    {}

    //! Evaluate spline value
    virtual
    real_type
    operator () ( real_type x, real_type y ) const UTILS_OVERRIDE;

    //! First derivative
    virtual
    void
    D( real_type x, real_type y, real_type d[3] ) const UTILS_OVERRIDE;

    virtual
    real_type
    Dx( real_type x, real_type y ) const UTILS_OVERRIDE;

    virtual
    real_type
    Dy( real_type x, real_type y ) const UTILS_OVERRIDE;

    //! Second derivative
    virtual
    void
    DD( real_type x, real_type y, real_type dd[6] ) const UTILS_OVERRIDE;

    virtual
    real_type
    Dxx( real_type , real_type ) const UTILS_OVERRIDE
    { return 0; }

    virtual
    real_type
    Dxy( real_type , real_type ) const UTILS_OVERRIDE
    { return 0; }

    virtual
    real_type
    Dyy( real_type , real_type ) const UTILS_OVERRIDE
    { return 0; }

    //! Print spline coefficients
    virtual
    void
    writeToStream( ostream_type & s ) const UTILS_OVERRIDE;

    //! Return spline typename
    virtual
    char const *
    type_name() const UTILS_OVERRIDE;

  };

  /*\
   |   ____  _  ____      _     _      ____        _ _            ____
   |  | __ )(_)/ ___|   _| |__ (_) ___/ ___| _ __ | (_)_ __   ___| __ )  __ _ ___  ___
   |  |  _ \| | |  | | | | '_ \| |/ __\___ \| '_ \| | | '_ \ / _ \  _ \ / _` / __|/ _ \
   |  | |_) | | |__| |_| | |_) | | (__ ___) | |_) | | | | | |  __/ |_) | (_| \__ \  __/
   |  |____/|_|\____\__,_|_.__/|_|\___|____/| .__/|_|_|_| |_|\___|____/ \__,_|___/\___|
   |                                        |_|
  \*/
  //! Bi-cubic spline base class
  class BiCubicSplineBase : public SplineSurf {
  protected:

    Utils::Malloc<real_type> m_mem_bicubic;

    real_type * m_DX;
    real_type * m_DY;
    real_type * m_DXY;

    void load( integer i, integer j, real_type bili3[4][4] ) const;

  public:

    using SplineSurf::m_nx;
    using SplineSurf::m_ny;

    using SplineSurf::m_X;
    using SplineSurf::m_Y;
    using SplineSurf::m_Z;

    //! spline constructor
    BiCubicSplineBase( string const & name = "Spline" )
    : SplineSurf( name )
    , m_mem_bicubic("BiCubicSplineBase")
    , m_DX(nullptr)
    , m_DY(nullptr)
    , m_DXY(nullptr)
    {}

    virtual
    ~BiCubicSplineBase() UTILS_OVERRIDE
    {}

    real_type
    DxNode ( integer i, integer j ) const
    { return m_DX[size_t(this->ipos_C(i,j))]; }

    real_type
    DyNode ( integer i, integer j ) const
    { return m_DY[size_t(this->ipos_C(i,j))]; }

    real_type
    DxyNode( integer i, integer j ) const
    { return m_DXY[size_t(this->ipos_C(i,j))]; }

    //! Evaluate spline value
    virtual
    real_type
    operator () ( real_type x, real_type y ) const UTILS_OVERRIDE;

    //! First derivative
    virtual
    void
    D( real_type x, real_type y, real_type d[3] ) const UTILS_OVERRIDE;

    virtual
    real_type
    Dx( real_type x, real_type y ) const UTILS_OVERRIDE;

    virtual
    real_type
    Dy( real_type x, real_type y ) const UTILS_OVERRIDE;

    //! Second derivative
    virtual
    void
    DD( real_type x, real_type y, real_type dd[6] ) const UTILS_OVERRIDE;

    virtual
    real_type
    Dxx( real_type x, real_type y ) const UTILS_OVERRIDE;

    virtual
    real_type
    Dxy( real_type x, real_type y ) const UTILS_OVERRIDE;

    virtual
    real_type
    Dyy( real_type x, real_type y ) const UTILS_OVERRIDE;
  };

  /*\
   |   ____  _  ____      _     _      ____        _ _
   |  | __ )(_)/ ___|   _| |__ (_) ___/ ___| _ __ | (_)_ __   ___
   |  |  _ \| | |  | | | | '_ \| |/ __\___ \| '_ \| | | '_ \ / _ \
   |  | |_) | | |__| |_| | |_) | | (__ ___) | |_) | | | | | |  __/
   |  |____/|_|\____\__,_|_.__/|_|\___|____/| .__/|_|_|_| |_|\___|
   |                                        |_|
  \*/
  //! cubic spline base class
  class BiCubicSpline : public BiCubicSplineBase {
    virtual void makeSpline() UTILS_OVERRIDE;

  public:

    using BiCubicSplineBase::m_mem_bicubic;
    using BiCubicSplineBase::m_DX;
    using BiCubicSplineBase::m_DY;
    using BiCubicSplineBase::m_DXY;

    //! spline constructor
    BiCubicSpline( string const & name = "Spline" )
    : BiCubicSplineBase( name )
    {}

    virtual
    ~BiCubicSpline() UTILS_OVERRIDE
    {}

    //! Print spline coefficients
    virtual
    void
    writeToStream( ostream_type & s ) const UTILS_OVERRIDE;

    //! Return spline typename
    virtual
    char const *
    type_name() const UTILS_OVERRIDE;

  };

  /*\
   |      _    _    _                 ____  ____            _ _
   |     / \  | | _(_)_ __ ___   __ _|___ \|  _ \ ___ _ __ | (_)_ __   ___
   |    / _ \ | |/ / | '_ ` _ \ / _` | __) | | | / __| '_ \| | | '_ \ / _ \
   |   / ___ \|   <| | | | | | | (_| |/ __/| |_| \__ \ |_) | | | | | |  __/
   |  /_/   \_\_|\_\_|_| |_| |_|\__,_|_____|____/|___/ .__/|_|_|_| |_|\___|
   |                                                 |_|
  \*/
  //! cubic spline base class
  class Akima2Dspline : public BiCubicSplineBase {
    virtual void makeSpline() UTILS_OVERRIDE;

  public:

    //! spline constructor
    Akima2Dspline( string const & name = "Spline" )
    : BiCubicSplineBase( name )
    {}

    virtual
    ~Akima2Dspline() UTILS_OVERRIDE
    {}

    //! Print spline coefficients
    virtual
    void
    writeToStream( ostream_type & s ) const UTILS_OVERRIDE;

    //! Return spline typename
    virtual
    char const *
    type_name() const UTILS_OVERRIDE;

  };

  /*\
   |   ____  _  ___        _       _   _      ____        _ _            ____
   |  | __ )(_)/ _ \ _   _(_)_ __ | |_(_) ___/ ___| _ __ | (_)_ __   ___| __ )  __ _ ___  ___
   |  |  _ \| | | | | | | | | '_ \| __| |/ __\___ \| '_ \| | | '_ \ / _ \  _ \ / _` / __|/ _ \
   |  | |_) | | |_| | |_| | | | | | |_| | (__ ___) | |_) | | | | | |  __/ |_) | (_| \__ \  __/
   |  |____/|_|\__\_\\__,_|_|_| |_|\__|_|\___|____/| .__/|_|_|_| |_|\___|____/ \__,_|___/\___|
   |                                               |_|
  \*/
  //! Bi-quintic spline base class
  class BiQuinticSplineBase : public SplineSurf {
  protected:

    Utils::Malloc<real_type> mem;

    real_type * m_DX;
    real_type * m_DXX;
    real_type * m_DY;
    real_type * m_DYY;
    real_type * m_DXY;
    real_type * m_DXYY;
    real_type * m_DXXY;
    real_type * m_DXXYY;
    void load( integer i, integer j, real_type bili5[6][6] ) const;

  public:

    //! spline constructor
    BiQuinticSplineBase( string const & name = "Spline" )
    : SplineSurf( name )
    , mem("BiQuinticSplineBase")
    , m_DX(nullptr)
    , m_DXX(nullptr)
    , m_DY(nullptr)
    , m_DYY(nullptr)
    , m_DXY(nullptr)
    , m_DXYY(nullptr)
    , m_DXXY(nullptr)
    {}

    virtual
    ~BiQuinticSplineBase() UTILS_OVERRIDE
    { mem.free(); }

    real_type
    DxNode( integer i, integer j ) const
    { return m_DX[size_t(this->ipos_C(i,j))]; }

    real_type
    DyNode( integer i, integer j ) const
    { return m_DY[size_t(this->ipos_C(i,j))]; }

    real_type
    DxxNode( integer i, integer j ) const
    { return m_DXX[size_t(this->ipos_C(i,j))]; }

    real_type
    DyyNode( integer i, integer j ) const
    { return m_DYY[size_t(this->ipos_C(i,j))]; }

    real_type
    DxyNode( integer i, integer j ) const
    { return m_DXY[size_t(this->ipos_C(i,j))]; }

    //! Evaluate spline value
    virtual
    real_type
    operator () ( real_type x, real_type y ) const UTILS_OVERRIDE;

    //! First derivative
    virtual
    void
    D( real_type x, real_type y, real_type d[3] ) const UTILS_OVERRIDE;

    virtual
    real_type
    Dx( real_type x, real_type y ) const UTILS_OVERRIDE;

    virtual
    real_type
    Dy( real_type x, real_type y ) const UTILS_OVERRIDE;

    //! Second derivative
    virtual
    void
    DD( real_type x, real_type y, real_type dd[6] ) const UTILS_OVERRIDE;

    virtual
    real_type
    Dxx( real_type x, real_type y ) const UTILS_OVERRIDE;

    virtual
    real_type
    Dxy( real_type x, real_type y ) const UTILS_OVERRIDE;

    virtual
    real_type
    Dyy( real_type x, real_type y ) const UTILS_OVERRIDE;
  };

  /*\
   |   ____  _  ___        _       _   _      ____        _ _
   |  | __ )(_)/ _ \ _   _(_)_ __ | |_(_) ___/ ___| _ __ | (_)_ __   ___
   |  |  _ \| | | | | | | | | '_ \| __| |/ __\___ \| '_ \| | | '_ \ / _ \
   |  | |_) | | |_| | |_| | | | | | |_| | (__ ___) | |_) | | | | | |  __/
   |  |____/|_|\__\_\\__,_|_|_| |_|\__|_|\___|____/| .__/|_|_|_| |_|\___|
   |                                               |_|
  \*/
  //! cubic spline base class
  class BiQuinticSpline : public BiQuinticSplineBase {
    virtual void makeSpline() UTILS_OVERRIDE;
  public:

    //! spline constructor
    BiQuinticSpline( string const & name = "Spline" )
    : BiQuinticSplineBase( name )
    {}

    virtual
    ~BiQuinticSpline() UTILS_OVERRIDE
    {}

    //! Print spline coefficients
    virtual
    void
    writeToStream( ostream_type & s ) const UTILS_OVERRIDE;

    //! Return spline typename
    virtual
    char const *
    type_name() const UTILS_OVERRIDE;

  };

  /*\
   |   ____        _ _            ____  ____
   |  / ___| _ __ | (_)_ __   ___|___ \|  _ \
   |  \___ \| '_ \| | | '_ \ / _ \ __) | | | |
   |   ___) | |_) | | | | | |  __// __/| |_| |
   |  |____/| .__/|_|_|_| |_|\___|_____|____/
   |        |_|
  \*/

  //! Bi-quintic spline base class
  class Spline2D {
  protected:
    std::string  m_name;
    SplineSurf * m_pSpline2D;
  public:

    //! spline constructor
    Spline2D( string const & name = "Spline2D" )
    : m_name(name)
    , m_pSpline2D( nullptr )
    {}

    ~Spline2D()
    {}

    bool is_x_closed() const { return m_pSpline2D->is_x_closed(); }
    void make_x_closed()     { m_pSpline2D->make_x_closed(); }
    void make_x_opened()     { m_pSpline2D->make_x_opened(); }

    bool is_y_closed() const { return m_pSpline2D->is_y_closed(); }
    void make_y_closed()     { m_pSpline2D->make_y_closed(); }
    void make_y_opened()     { m_pSpline2D->make_y_opened(); }

    bool is_x_bounded() const { return m_pSpline2D->is_x_bounded(); }
    void make_x_unbounded()   { m_pSpline2D->make_x_unbounded(); }
    void make_x_bounded()     { m_pSpline2D->make_x_bounded(); }

    bool is_y_bounded() const { return m_pSpline2D->is_y_bounded(); }
    void make_y_unbounded()   { m_pSpline2D->make_y_unbounded(); }
    void make_y_bounded()     { m_pSpline2D->make_y_bounded(); }

    string const & name() const { return m_pSpline2D->name(); }

    //! Cancel the support points, empty the spline.
    void clear(void) { m_pSpline2D->clear(); }

    //! return the number of support points of the spline along x direction
    integer
    numPointX(void) const { return m_pSpline2D->numPointX(); }

    //! return the number of support points of the spline along y direction
    integer
    numPointY(void) const { return m_pSpline2D->numPointY(); }

    //! return the i-th node of the spline (x component).
    real_type
    xNode( integer i ) const { return m_pSpline2D->xNode(i); }

    //! return the i-th node of the spline (y component).
    real_type
    yNode( integer i ) const { return m_pSpline2D->yNode(i); }

    //! return the i-th node of the spline (y component).
    real_type
    zNode( integer i, integer j ) const { return m_pSpline2D->zNode(i,j); }

    //! return x-minumum spline value
    real_type
    xMin() const { return m_pSpline2D->xMin(); }

    //! return x-maximum spline value
    real_type
    xMax() const { return m_pSpline2D->xMax(); }

    //! return y-minumum spline value
    real_type
    yMin() const { return m_pSpline2D->yMin(); }

    //! return y-maximum spline value
    real_type
    yMax() const { return m_pSpline2D->yMax(); }

    //! return z-minumum spline value
    real_type
    zMin() const { return m_pSpline2D->zMin(); }

    //! return z-maximum spline value
    real_type
    zMax() const { return m_pSpline2D->zMax(); }

    void
    build(
      SplineType2D    tp,
      real_type const x[], integer incx,
      real_type const y[], integer incy,
      real_type const z[], integer ldZ,
      integer         nx,
      integer         ny,
      bool            fortran_storage = false,
      bool            transposed      = false
    );

    /*! Build surface spline
     | \param x       vector of x-coordinates, nx = x.size()
     | \param y       vector of y-coordinates, ny = y.size()
     | \param z       matrix of z-values. Elements are stored
     |                by row Z(i,j) = z[i*ny+j] as C-matrix
     | \param fortran_storage if true elements are stored by column
     |                        i.e. Z(i,j) = z[i+j*nx] as Fortran-matrix
     | \param transposed      if true matrix Z is stored transposed
    \*/
    void
    build(
      SplineType2D              tp,
      vector<real_type> const & x,
      vector<real_type> const & y,
      vector<real_type> const & z,
      bool fortran_storage = false,
      bool transposed      = false
    );

    /*! Build surface spline
     | \param z               matrix of z-values. Elements are stored
     |                        by row Z(i,j) = z[i*ny+j] as C-matrix
     | \param ldZ             leading dimension of the matrix. Elements are stored
     |                        by row Z(i,j) = z[i*ldZ+j] as C-matrix
     | \param fortran_storage if true elements are stored by column
     |                        i.e. Z(i,j) = z[i+j*nx] as Fortran-matrix
     | \param transposed      if true matrix Z is stored transposed
    \*/
    void
    build(
      SplineType2D    tp,
      real_type const z[],
      integer         ldZ,
      integer         nx,
      integer         ny,
      bool fortran_storage = false,
      bool transposed      = false
    );

    /*! Build surface spline
     | \param z               matrix of z-values. Elements are stored
     |                        by row Z(i,j) = z[i*ny+j] as C-matrix.
     |                        ldZ leading dimension of the matrix is ny for C-storage
     |                        and nx for Fortran storage.
     | \param fortran_storage if true elements are stored by column
     |                        i.e. Z(i,j) = z[i+j*nx] as Fortran-matrix
     | \param transposed      if true matrix Z is stored transposed
    \*/
    void
    build(
      SplineType2D              tp,
      vector<real_type> const & z,
      integer                   nx,
      integer                   ny,
      bool fortran_storage = false,
      bool transposed      = false
    );

    void
    setup( GenericContainer const & gc )
    { build(gc); }

    void
    build( GenericContainer const & gc );

    //! Evaluate spline value
    real_type
    operator () ( real_type x, real_type y ) const
    { return (*m_pSpline2D)( x, y ); }

    //! First derivative
    void
    D( real_type x, real_type y, real_type d[3] ) const
    { return m_pSpline2D->D( x, y, d ); }

    real_type
    Dx( real_type x, real_type y ) const
    { return m_pSpline2D->Dx( x, y ); }

    real_type
    Dy( real_type x, real_type y ) const
    { return m_pSpline2D->Dy( x, y ); }

    //! Second derivative
    void
    DD( real_type x, real_type y, real_type dd[6] ) const
    { return m_pSpline2D->DD( x, y, dd ); }

    real_type
    Dxx( real_type x, real_type y ) const
    { return m_pSpline2D->Dxx( x, y ); }

    real_type
    Dxy( real_type x, real_type y ) const
    { return m_pSpline2D->Dxy( x, y ); }

    real_type
    Dyy( real_type x, real_type y ) const
    { return m_pSpline2D->Dyy( x, y ); }

    //! Evaluate spline value
    real_type
    eval( real_type x, real_type y ) const
    { return (*this)(x,y); }

    //! First derivative
    real_type
    eval_D_1( real_type x, real_type y ) const
    { return this->Dx(x,y); }

    real_type
    eval_D_2( real_type x, real_type y ) const
    { return this->Dy(x,y); }

    //! Second derivative
    real_type
    eval_D_1_1( real_type x, real_type y ) const
    { return this->Dxx(x,y); }

    real_type
    eval_D_1_2( real_type x, real_type y ) const
    { return this->Dxy(x,y); }

    real_type
    eval_D_2_2( real_type x, real_type y ) const
    { return this->Dyy(x,y); }

    //! Print spline coefficients
    void
    writeToStream( ostream_type & s ) const
    { return m_pSpline2D->writeToStream( s ); }

    //! Return spline typename
    char const * type_name() const { return m_pSpline2D->type_name(); }

    void info( ostream_type & s ) const { m_pSpline2D->info( s ); }

  };

}

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
