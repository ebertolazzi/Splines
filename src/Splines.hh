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

/*!
  \mainpage  Splines
  \author    Enrico Bertolazzi (enrico.bertolazzi@unitn.it), homepage: http://www.ing.unitn.it/~bertolaz
  \version   1.0.7
  \note      first release Jan 12, 1998
  \date      2020
  \copyright see License.txt.
  \par       Affiliation:
             Department of Industrial Engineering<br>
             University of Trento<br>
             Via Sommarive 9, I -- 38123 Trento, Italy <br>
             enrico.bertolazzi@unitn.it
*/

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

  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  extern SplineType1D string_to_splineType( string const & n );
  #endif

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

  /*!
   * Compute nodes for the spline using uniform distribution
   *
   * \param[in]  dim     dimension of the points
   * \param[in]  npts    number of points
   * \param[in]  pnts    matrix whose columns are the points
   * \param[in]  ld_pnts leading dimension of the matrix (fortran storage)  
   * \param[out] t       vector of the computed nodes
   * 
   */
  void
  uniform(
    integer           dim,
    integer           npts,
    real_type const * pnts,
    integer           ld_pnts,
    real_type       * t
  );

  /*!
   * Compute nodes for the spline using chordal distribution
   *
   * \param[in]  dim     dimension of the points
   * \param[in]  npts    number of points
   * \param[in]  pnts    matrix whose columns are the points
   * \param[in]  ld_pnts leading dimension of the matrix (fortran storage)  
   * \param[out] t       vector of the computed nodes
   * 
   */
  void
  chordal(
    integer           dim,
    integer           npts,
    real_type const * pnts,
    integer           ld_pnts,
    real_type       * t
  );

  /*!
   * Compute nodes for the spline using centripetal distribution
   *
   * \param[in]  dim     dimension of the points
   * \param[in]  npts    number of points
   * \param[in]  pnts    matrix whose columns are the points
   * \param[in]  ld_pnts leading dimension of the matrix (fortran storage)  
   * \param[in]  alpha   power factor
   * \param[out] t       vector of the computed nodes
   * 
   */
  void
  centripetal(
    integer           dim,
    integer           npts,
    real_type const * pnts,
    integer           ld_pnts,
    real_type         alpha,
    real_type       * t
  );

  /*!
   * Compute nodes for the spline using universal distribution
   *
   * \param[in]  dim     dimension of the points
   * \param[in]  npts    number of points
   * \param[in]  pnts    matrix whose columns are the points
   * \param[in]  ld_pnts leading dimension of the matrix (fortran storage)  
   * \param[out] t       vector of the computed nodes
   * 
   */
  void
  universal(
    integer           dim,
    integer           npts,
    real_type const * pnts,
    integer           ld_pnts,
    real_type       * t
  );

  /*!
   * Compute nodes for the spline using ``FoleyNielsen`` distribution
   *
   * \param[in]  dim     dimension of the points
   * \param[in]  npts    number of points
   * \param[in]  pnts    matrix whose columns are the points
   * \param[in]  ld_pnts leading dimension of the matrix (fortran storage)  
   * \param[out] t       vector of the computed nodes
   * 
   */
  void
  FoleyNielsen(
    integer           dim,
    integer           npts,
    real_type const * pnts,
    integer           ld_pnts,
    real_type       * t
  );

  /*!
   * Compute nodes for the spline using ``FangHung`` distribution
   *
   * \param[in]  dim     dimension of the points
   * \param[in]  npts    number of points
   * \param[in]  pnts    matrix whose columns are the points
   * \param[in]  ld_pnts leading dimension of the matrix (fortran storage)  
   * \param[out] t       vector of the computed nodes
   * 
   */
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
  //! Spline Management Class
  class Spline {
    friend class SplineSet;
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

    //! find interval containing ``x`` using binary search.
    integer search( real_type & x ) const;

    /*!
     * \name Spline Data Info
     */
    //@{

    //! return the name of the spline used in the constructor
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
    integer numPoints() const { return m_npts; }

    //! return the i-th node of the spline (x component).
    real_type xNode( integer i ) const { return m_X[size_t(i)]; }

    //! return the i-th node of the spline (y component).
    real_type yNode( integer i ) const { return m_Y[size_t(i)]; }

    //! return first node of the spline (x component).
    real_type xBegin() const { return m_X[0]; }

    //! return first node of the spline (y component).
    real_type yBegin() const { return m_Y[0]; }

    //! return last node of the spline (x component).
    real_type xEnd() const { return m_X[size_t(m_npts-1)]; }

    //! return last node of the spline (y component).
    real_type yEnd() const { return m_Y[size_t(m_npts-1)]; }

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

    //@}

    /*!
     * \name Build
     */
    //@{

    /*!
     * Build a spline using a generic container.
     */
    void
    build( GenericContainer const & gc )
    { setup(gc); }

    /*!
     * \brief Build a spline.
     * 
     * \param x    vector of x-coordinates
     * \param incx access elements as x[0], x[incx], x[2*incx],...
     * \param y    vector of y-coordinates
     * \param incy access elements as y[0], y[incy], y[2*incy],...
     * \param n    total number of points
     */
    virtual
    void
    build(
      real_type const * x, integer incx,
      real_type const * y, integer incy,
      integer n
    );

    /*!
     * \brief Build a spline.
     * 
     * \param x vector of x-coordinates
     * \param y vector of y-coordinates
     * \param n total number of points
     */
    inline
    void
    build(
      real_type const * x,
      real_type const * y,
      integer           n
    )
    { this->build( x, 1, y, 1, n ); }

    /*!
     * \brief Build a spline.
     * 
     * \param x vector of x-coordinates
     * \param y vector of y-coordinates
     */
    inline
    void
    build( vector<real_type> const & x, vector<real_type> const & y ) {
      integer N = integer(x.size());
      if ( N > integer(y.size()) ) N = integer(y.size());
      this->build( &x.front(), 1, &y.front(), 1, N );
    }

    /*!
     * Build a spline using internal stored data 
     */
    virtual
    void build() = 0;

    /*!
     * Build a spline using a generic container.
     */
    virtual
    void setup( GenericContainer const & gc );

    //@}

    /*!
     * \name Incremental Build
     * 
     * Various constructors for the spline class.
     * 
     */
    //@{

    //! Allocate memory for `npts` points
    virtual
    void reserve( integer npts ) = 0;

    //! Add a support point (x,y) to the spline.
    void pushBack( real_type x, real_type y );

    //! Drop last inserted point of the spline.
    void dropBack() { if ( m_npts > 0 ) --m_npts; }

    //! Cancel the support points, empty the spline.
    virtual
    void clear() = 0;

    //@}

    /*!
     * \name Dump Data
     */
    //@{

    //! dump a sample of the spline
    void
    dump(
      ostream_type & s,
      integer        nintervals,
      char const *   header = "x\ty"
    ) const;

    //! dump a sample of the spline
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

    //! Print spline coefficients
    virtual
    void writeToStream( ostream_type & s ) const = 0;

    //@}

    /*!
     * \name Evaluation
     */
    //@{

    //! Evaluate spline value
    virtual
    real_type operator () ( real_type x ) const = 0;

    //! First derivative
    virtual
    real_type D( real_type x ) const = 0;

    //! Second derivative
    virtual
    real_type DD( real_type x ) const = 0;

    //! Third derivative
    virtual
    real_type DDD( real_type x ) const = 0;

    //! 4th derivative
    virtual
    real_type DDDD( real_type ) const { return real_type(0); }

    //! 4th derivative
    virtual
    real_type DDDDD( real_type ) const { return real_type(0); }

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
    real_type id_eval( integer ni, real_type x ) const = 0;

    //! First derivative
    virtual
    real_type id_D( integer ni, real_type x ) const = 0;

    //! Second derivative
    virtual
    real_type id_DD( integer ni, real_type x ) const = 0;

    //! Third derivative
    virtual
    real_type id_DDD( integer ni, real_type x ) const = 0;

    //! 4th derivative
    virtual
    real_type id_DDDD( integer, real_type ) const { return real_type(0); }

    //! 4th derivative
    virtual
    real_type id_DDDDD( integer, real_type ) const { return real_type(0); }

    //@}

    /*!
     * \name Get Info
     */
    //@{

    //! get the piecewise polinomials of the spline
    virtual
    integer // order
    coeffs(
      real_type * const cfs,
      real_type * const nodes,
      bool              transpose = false
    ) const = 0;

    virtual
    integer
    order() const = 0;
    //! Return spline typename
    char const *
    type_name() const
    { return Splines::spline_type_1D[type()]; }

    //! Return spline type (as number)
    virtual
    unsigned type() const = 0;

    //! Return a string information of the kind and order of the spline
    string info() const;

    //! Print information of the kind and order of the spline
    void info( ostream_type & stream ) const { stream << this->info() << '\n'; }

    //@}

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

    #ifndef DOXYGEN_SHOULD_SKIP_THIS
    using Spline::build;
    #endif

    //! spline constructor
    CubicSplineBase( string const & name = "CubicSplineBase" )
    : Spline(name)
    , m_baseValue(name+"_memory")
    , m_Yp(nullptr)
    , m_external_alloc(false)
    {}

    ~CubicSplineBase() override {}

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

    real_type operator () ( real_type x ) const override;
    real_type D( real_type x ) const override;
    real_type DD( real_type x ) const override;
    real_type DDD( real_type x ) const override;
    real_type id_eval( integer ni, real_type x ) const override;
    real_type id_D( integer ni, real_type x ) const override;
    real_type id_DD( integer ni, real_type x ) const override;
    real_type id_DDD( integer ni, real_type x ) const override;

    void writeToStream( ostream_type & s ) const override;

    // --------------------------- VIRTUALS -----------------------------------

    void reserve( integer npts ) override;

    /*!
     * \brief Build a spline.
     *
     * \param x     vector of x-coordinates
     * \param incx  access elements as x[0], x[incx], x[2*incx],...
     * \param y     vector of y-coordinates
     * \param incy  access elements as y[0], y[incy], y[2*incy],...
     * \param yp    vector of y-defivative
     * \param incyp access elements as yp[0], yp[incyp], yp[2*incyy],...
     * \param n     total number of points
     */
    // must be defined in derived classes
    void
    build(
      real_type const * x,  integer incx,
      real_type const * y,  integer incy,
      real_type const * yp, integer incyp,
      integer n
    );

    /*!
     * \brief Build a spline.
     *
     * \param x  vector of x-coordinates
     * \param y  vector of y-coordinates
     * \param yp vector of y'-coordinates
     * \param n  total number of points
     */
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

    /*!
     * \brief Build a spline.
     *
     * \param x  vector of x-coordinates
     * \param y  vector of y-coordinates
     * \param yp vector of y'-coordinates
     */
    void
    build(
      vector<real_type> const & x,
      vector<real_type> const & y,
      vector<real_type> const & yp
    );

    //! Cancel the support points, empty the spline.
    void clear() override;

    //! get the piecewise polinomials of the spline
    integer // order
    coeffs(
      real_type * const cfs,
      real_type * const nodes,
      bool              transpose = false
    ) const override;

    integer order() const override;

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

    virtual void makeSpline() = 0;

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
    void clear();

    //! return the number of support points of the spline along x direction
    integer numPointX() const { return m_nx; }

    //! return the number of support points of the spline along y direction
    integer numPointY() const { return m_ny; }

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
     * \param x       vector of x-coordinates
     * \param incx    access elements as x[0], x[incx], x[2*incx],...
     * \param y       vector of y-coordinates
     * \param incy    access elements as y[0], y[incy], x[2*incy],...
     * \param z       matrix of z-values
     * \param ldZ     leading dimension of the matrix. Elements are stored
     *                by row Z(i,j) = z[i*ldZ+j] as C-matrix
     * \param nx      total number of points in direction x
     * \param ny      total number of points in direction y
     * \param fortran_storage if true elements are stored by column
     *                        i.e. Z(i,j) = z[i+j*ldZ] as Fortran-matrix
     * \param transposed      if true matrix Z is stored transposed
     */
    void
    build(
      real_type const * x, integer incx,
      real_type const * y, integer incy,
      real_type const * z, integer ldZ,
      integer nx, integer ny,
      bool fortran_storage = false,
      bool transposed      = false
    );

    /*! Build surface spline
     * \param x       vector of x-coordinates, nx = x.size()
     * \param y       vector of y-coordinates, ny = y.size()
     * \param z       matrix of z-values. Elements are stored
     *                by row Z(i,j) = z[i*ny+j] as C-matrix
     * \param fortran_storage if true elements are stored by column
     *                        i.e. Z(i,j) = z[i+j*nx] as Fortran-matrix
     * \param transposed      if true matrix Z is stored transposed
     */
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

    /*!
     * \brief Build surface spline
     *
     * \param z               matrix of z-values. Elements are stored
     *                        by row Z(i,j) = z[i*ny+j] as C-matrix
     * \param ldZ             leading dimension of the matrix. Elements are stored
     *                        by row Z(i,j) = z[i*ldZ+j] as C-matrix
     * \param nx              x-dimension
     * \param ny              y-dimension
     * \param fortran_storage if true elements are stored by column
     *                        i.e. Z(i,j) = z[i+j*nx] as Fortran-matrix
     * \param transposed      if true matrix Z is stored transposed
     */
    void
    build(
      real_type const * z,
      integer           ldZ,
      integer           nx,
      integer           ny,
      bool fortran_storage = false,
      bool transposed      = false
    );

    /*!
     * \brief Build surface spline
     * 
     * \param z               matrix of z-values. Elements are stored
     *                        by row Z(i,j) = z[i*ny+j] as C-matrix.
     *                        ldZ leading dimension of the matrix is ny for C-storage
     *                        and nx for Fortran storage.
     * \param nx              x-dimension
     * \param ny              y-dimension
     * \param fortran_storage if true elements are stored by column
     *                        i.e. Z(i,j) = z[i+j*nx] as Fortran-matrix
     * \param transposed      if true matrix Z is stored transposed
     */
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

    //! Evaluate spline value at point \f$ (x,y) \f$
    virtual
    real_type
    operator () ( real_type x, real_type y ) const = 0;

    /*!
     * value and first derivatives at point \f$ (x,y) \f$
     * - d[0] value of the spline \f$ S(x,y) \f$
     * - d[1] derivative respect to \f$ x \f$ of the spline: \f$ S_x(x,y) \f$
     * - d[2] derivative respect to \f$ y \f$ of the spline: \f$ S_y(x,y) \f$
     */
    virtual
    void
    D( real_type x, real_type y, real_type d[3] ) const = 0;

    /*!
     * first derivatives respect to \f$ x \f$ at point \f$ (x,y) \f$ of the spline: \f$ S_x(x,y) \f$
     */
    virtual
    real_type
    Dx( real_type x, real_type y ) const = 0;

    /*!
     * first derivatives respect to \f$ y \f$ at point \f$ (x,y) \f$ of the spline: \f$ S_y(x,y) \f$
     */
    virtual
    real_type
    Dy( real_type x, real_type y ) const = 0;

    /*!
     * value, first and second derivatives at point \f$ (x,y) \f$
     * - d[0] value of the spline \f$ S(x,y) \f$
     * - d[1] derivative respect to \f$ x \f$ of the spline: \f$ S_x(x,y) \f$
     * - d[2] derivative respect to \f$ y \f$ of the spline: \f$ S_y(x,y) \f$
     * - d[3] second derivative respect to \f$ x \f$ of the spline: \f$ S_{xx}(x,y) \f$
     * - d[4] mixed second derivative: \f$ S_{xy}(x,y) \f$
     * - d[5] second derivative respect to \f$ y \f$ of the spline: \f$ S_{yy}(x,y) \f$
     */
    virtual
    void
    DD( real_type x, real_type y, real_type dd[6] ) const = 0;

    /*!
     * second derivatives respect to \f$ x \f$ at point \f$ (x,y) \f$ of the spline: \f$ S_{xx}(x,y) \f$
     */
    virtual
    real_type
    Dxx( real_type x, real_type y ) const = 0;

    /*!
     * mixed second derivatives: \f$ S_{xy}(x,y) \f$
     */
    virtual
    real_type
    Dxy( real_type x, real_type y ) const = 0;

    /*!
     * second derivatives respect to \f$ y \f$ at point \f$ (x,y) \f$ of the spline: \f$ S_{yy}(x,y) \f$
     */
    virtual
    real_type
    Dyy( real_type x, real_type y ) const = 0;

    //! Evaluate spline value at point \f$ (x,y) \f$
    real_type
    eval( real_type x, real_type y ) const
    { return (*this)(x,y); }

    //! Alias for ``Dx(x,y)``
    real_type
    eval_D_1( real_type x, real_type y ) const
    { return this->Dx(x,y); }

    //! Alias for ``Dy(x,y)``
    real_type
    eval_D_2( real_type x, real_type y ) const
    { return this->Dy(x,y); }

    //! Alias for ``Dxx(x,y)``
    real_type
    eval_D_1_1( real_type x, real_type y ) const
    { return this->Dxx(x,y); }

    //! Alias for ``Dxy(x,y)``
    real_type
    eval_D_1_2( real_type x, real_type y ) const
    { return this->Dxy(x,y); }

    //! Alias for ``Dyy(x,y)``
    real_type
    eval_D_2_2( real_type x, real_type y ) const
    { return this->Dyy(x,y); }

    //! Print spline coefficients
    virtual void writeToStream( ostream_type & s ) const = 0;

    //! Return spline typename
    virtual char const * type_name() const = 0;

    //! return a string with information about the spline.
    virtual string info() const;

    //! write a string with information about the spline.
    void
    info( ostream_type & stream ) const
    { stream << this->info() << '\n'; }

  };

}

#include "SplineAkima.hxx"
#include "SplineBessel.hxx"
#include "SplineConstant.hxx"
#include "SplineLinear.hxx"
#include "SplineCubic.hxx"
#include "SplineHermite.hxx"
#include "SplinePchip.hxx"
#include "SplineQuinticBase.hxx"
#include "SplineQuintic.hxx"
#include "SplineBilinear.hxx"
#include "SplineBiCubic.hxx"
#include "SplineAkima2D.hxx"
#include "SplineBiQuintic.hxx"

#include "SplineVec.hxx"
#include "SplineSet.hxx"
#include "Splines1D.hxx"
#include "Splines2D.hxx"

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
