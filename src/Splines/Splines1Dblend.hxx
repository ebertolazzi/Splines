/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2025                                                      |
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

/*\
 |   ____        _ _            _ ____
 |  / ___| _ __ | (_)_ __   ___/ |  _ \
 |  \___ \| '_ \| | | '_ \ / _ \ | | | |
 |   ___) | |_) | | | | | |  __/ | |_| |
 |  |____/| .__/|_|_|_| |_|\___|_|____/
 |        |_|
\*/

namespace Splines {

  //! Spline Management Class
  class Spline1Dblend {
  protected:

    Spline1D m_spline0;
    Spline1D m_spline1;

    void check_compatibility() const;

  public:

    Spline1Dblend( Spline1Dblend const & ) = delete;
    Spline1Dblend const & operator = ( Spline1Dblend const & ) = delete;

    //! \name Constructors
    ///@{

    //!
    //! Build an empty spline of `Spline1D` type
    //!
    //! \param n the name of the spline
    //!
    Spline1Dblend( string_view n )
    : m_spline0(fmt::format("{}_0",n))
    , m_spline1(fmt::format("{}_1",n))
    {}

    //!
    //! Spline destructor.
    //!
    ~Spline1Dblend()
    {}

    ///@}

    Spline1D const & get_spline0() const { return m_spline0; }
    Spline1D const & get_spline1() const { return m_spline1; }

    //!
    //! Return the number of support points of the spline.
    //!
    integer num_points0() const { return m_spline0.num_points(); }
    integer num_points1() const { return m_spline1.num_points(); }

    //!
    //! Build a spline using data in `GenericContainer`
    //!
    void setup( GenericContainer const & gc );
    //!
    //! Build a spline using data in `GenericContainer`
    //!
    void build( GenericContainer const & gc ) { setup(gc); }

    //!
    //! Build a spline.
    //!
    //! \param tp0   spline type
    //! \param x0    vector of `x`-coordinates
    //! \param incx0 access elements as `x[0]`, `x[incx]`, `x[2*incx]`,...
    //! \param y0    vector of `y`-coordinates
    //! \param incy0 access elements as `y[0]`, `y[incy]`, `x[2*incy]`,...
    //! \param n0    total number of points
    //!
    //! \param tp1   spline type
    //! \param x1    vector of `x`-coordinates
    //! \param incx1 access elements as `x[0]`, `x[incx]`, `x[2*incx]`,...
    //! \param y1    vector of `y`-coordinates
    //! \param incy1 access elements as `y[0]`, `y[incy]`, `x[2*incy]`,...
    //! \param n1    total number of points
    //!
    // must be defined in derived classes
    void
    build(
      SplineType1D    tp0,
      real_type const x0[], integer const incx0,
      real_type const y0[], integer const incy0,
      integer   const n0,
      SplineType1D    tp1,
      real_type const x1[], integer const incx1,
      real_type const y1[], integer const incy1,
      integer   const n1
    ) {
      m_spline0.build( tp0, x0, incx0, y0, incy0, n0 );
      m_spline1.build( tp1, x1, incx1, y1, incy1, n1 );
      check_compatibility();
    }

    //!
    //! Build a spline.
    //!
    //! \param tp0   spline type
    //! \param x0    vector of `x`-coordinates
    //! \param y0    vector of `y`-coordinates
    //! \param n0    total number of points
    //!
    //! \param tp1   spline type
    //! \param x1    vector of `x`-coordinates
    //! \param y1    vector of `y`-coordinates
    //! \param n1    total number of points
    //!
    void
    build(
      SplineType1D    tp0,
      real_type const x0[],
      real_type const y0[],
      integer   const n0,
      SplineType1D    tp1,
      real_type const x1[],
      real_type const y1[],
      integer   const n1
    ) {
      m_spline0.build( tp0, x0, y0, n0 );
      m_spline1.build( tp1, x1, y1, n1 );
      check_compatibility();
    }

    //!
    //! Build a spline.
    //!
    //! \param tp0   spline type
    //! \param x0    vector of `x`-coordinates
    //! \param y0    vector of `y`-coordinates
    //! \param n0    total number of points
    //!
    //! \param tp1   spline type
    //! \param x1    vector of `x`-coordinates
    //! \param y1    vector of `y`-coordinates
    //! \param n1    total number of points
    //!
    void
    build(
      SplineType1D              tp0,
      vector<real_type> const & x0,
      vector<real_type> const & y0,
      SplineType1D              tp1,
      vector<real_type> const & x1,
      vector<real_type> const & y1
    ) {
      m_spline0.build( tp0, x0, y0 );
      m_spline1.build( tp1, x1, y1 );
      check_compatibility();
    }

    //!
    //! Return x-initial spline value.
    //!
    real_type x_begin( real_type const s ) const { return (1-s) * m_spline0.x_begin() + s * m_spline1.x_begin(); }

    //!
    //! Return x-final spline value.
    //!
    real_type x_end( real_type const s ) const { return (1-s) * m_spline0.x_end() + s * m_spline1.x_end(); }

    //!
    //! Return y-initial spline value
    //!
    real_type y_begin( real_type const s ) const { return (1-s) * m_spline0.y_begin() + s * m_spline1.y_begin(); }

    //!
    //! Return y-final spline value
    //!
    real_type y_end( real_type const s ) const { return (1-s) * m_spline0.y_end() + s * m_spline1.y_end(); }

    ///////////////////////////////////////////////////////////////////////////
    //!
    //! Change X-origin of the spline.
    //!
    void
    set_origin( real_type const x0 ) const {
      m_spline0.set_origin( x0 );
      m_spline1.set_origin( x0 );
    }

    //!
    //! Change X-range of the spline.
    //!
    void
    set_range( real_type const xmin, real_type const xmax ) {
      m_spline0.set_range( xmin, xmax );
      m_spline1.set_range( xmin, xmax );
    }

    ///////////////////////////////////////////////////////////////////////////
    //! \name Evaluation
    ///@{

    //!
    //! Evaluate spline value at `x`.
    //!
    real_type
    eval( real_type const x, real_type const s ) const {
      return (1-s) * m_spline0.eval(x) + s * m_spline1.eval(x);
    }

    //!
    //! First derivative at `x`.
    //!
    real_type
    D( real_type const x, real_type const s ) const {
      return (1-s) * m_spline0.D(x) + s * m_spline1.D(x);
    }

    //!
    //! Second derivative at `x`.
    //!
    real_type
    DD( real_type const x, real_type const s ) const {
      return (1-s) * m_spline0.DD(x) + s * m_spline1.DD(x);
    }

    //!
    //! Third derivative at `x`.
    //!
    real_type
    DDD( real_type const x, real_type const s ) const {
      return (1-s) * m_spline0.DDD(x) + s * m_spline1.DDD(x);
    }

    //!
    //! 4th derivative at `x`.
    //!
    real_type
    DDDD( real_type const x, real_type const s ) const {
      return (1-s) * m_spline0.DDDD(x) + s * m_spline1.DDDD(x);
    }

    //!
    //! 5th derivative at `x`.
    //!
    real_type
    DDDDD( real_type const x, real_type const s ) const {
      return (1-s) * m_spline0.DDDDD(x) + s * m_spline1.DDDDD(x);
    }

    void
    D( real_type const x, real_type const s, real_type dd[2] ) const {
      real_type d0[2], d1[2];
      m_spline0.D(x,d0);
      m_spline1.D(x,d1);
      dd[0] = (1-s)*d0[0] + s*d1[0];
      dd[1] = (1-s)*d0[1] + s*d1[1];
    }
    
    void
    DD( real_type const x, real_type const s, real_type dd[3] ) const {
      real_type d0[3], d1[3];
      m_spline0.DD(x,d0);
      m_spline1.DD(x,d1);
      dd[0] = (1-s)*d0[0] + s*d1[0];
      dd[1] = (1-s)*d0[1] + s*d1[1];
      dd[2] = (1-s)*d0[2] + s*d1[2];
    }

    ///@}

    #ifdef AUTODIFF_SUPPORT
    //!
    //! \name Autodiff
    //!
    ///@{
    autodiff::dual1st eval( autodiff::dual1st const & x, real_type const s ) const;
    autodiff::dual2nd eval( autodiff::dual2nd const & x, real_type const s ) const;

    template <typename T>
    autodiff::HigherOrderDual<autodiff::detail::DualOrder<T>::value,real_type>
    eval( T const & x, real_type const s ) const {
      autodiff::HigherOrderDual<autodiff::detail::DualOrder<T>::value,real_type> X{x};
      return eval( X, s );
    }
    ///@}
    #endif

    //!
    //! \name Evaluation Aliases
    //!
    ///@{
    //! the value of the spline at `x`
    real_type
    operator () ( real_type const x, real_type const s ) const
    { return (1-s)*m_spline0.eval(x) + s*m_spline1.eval(x); }

    //! the value of the first derivative of the spline at `x`
    real_type
    eval_D( real_type const x, real_type const s ) const
    { return (1-s)*m_spline0.D(x) + s*m_spline1.D(x); }

    //! the value of the second derivative of the spline at `x`
    real_type
    eval_DD( real_type const x, real_type const s ) const
    { return (1-s)*m_spline0.DD(x) + s*m_spline1.DD(x); }

    //! the value of the third derivative of the spline at `x`
    real_type
    eval_DDD( real_type const x, real_type const s ) const
    { return (1-s)*m_spline0.DDD(x) + s*m_spline1.DDD(x); }

    //! the value of the 4-th derivative of the spline at `x`
    real_type
    eval_DDDD( real_type const x, real_type const s ) const
    { return (1-s)*m_spline0.DDDD(x) + s*m_spline1.DDDD(x); }

    //! the value of the 5-th derivative of the spline at `x`
    real_type
    eval_DDDDD( real_type const x, real_type const s ) const
    { return (1-s)*m_spline0.DDDDD(x) + s*m_spline1.DDDDD(x); }
    ///@}

    //!
    //! \return the order of the spline (`degree+1`)
    //!
    integer order0() const { return m_spline0.order(); }
    integer order1() const { return m_spline1.order(); }

    //!
    //! Return spline typename.
    //!
    char const * type_name0() const { return m_spline0.type_name(); }
    char const * type_name1() const { return m_spline1.type_name(); }

    //!
    //! Return spline type (as number).
    //!
    SplineType1D type0() const { return m_spline0.type(); }
    SplineType1D type1() const { return m_spline1.type(); }

  };

}

// EOF: Spline1Dblend.hxx
