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
  class Spline1D {
  protected:

    std::string m_name;

    Spline * m_pSpline;

    Spline1D( Spline1D const & ) = delete;
    Spline1D const & operator = ( Spline1D const & ) = delete;

  public:

    //! \name Constructors
    ///@{

    //! spline constructor
    Spline1D( std::string const & n )
    : m_name(n)
    , m_pSpline(nullptr)
    {}

    //! spline destructor
    ~Spline1D()
    {}

    ///@}

    string const & name() const { return m_pSpline->name(); }

    bool is_closed() const { return m_pSpline->is_closed(); }
    void make_closed() { m_pSpline->make_closed(); }
    void make_opened() { m_pSpline->make_opened(); }

    bool is_bounded() const { return m_pSpline->is_bounded(); }
    void make_unbounded()   { m_pSpline->make_unbounded(); }
    void make_bounded()     { m_pSpline->make_bounded(); }

    //! return the number of support points of the spline.
    integer numPoints() const { return m_pSpline->numPoints(); }

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

    //! Drop last inserted point of the spline.
    void dropBack() { m_pSpline->dropBack(); }

    //! Build a spline.
    // must be defined in derived classes
    void build() { m_pSpline->build(); }

    void setup( GenericContainer const & gc );
    void build( GenericContainer const & gc ) { setup(gc); }

    //! 
    //! Build a spline.
    //! 
    //! \param tp   spline type
    //! \param x    vector of x-coordinates
    //! \param incx access elements as x[0], x[incx], x[2*incx],...
    //! \param y    vector of y-coordinates
    //! \param incy access elements as y[0], y[incy], x[2*incy],...
    //! \param n    total number of points
    //! 
    // must be defined in derived classes
    void
    build(
      SplineType1D tp,
      real_type const * x, integer incx,
      real_type const * y, integer incy,
      integer n
    );

    //! 
    //! Build a spline.
    //! 
    //! \param tp spline type
    //! \param x  vector of x-coordinates
    //! \param y  vector of y-coordinates
    //! \param n  total number of points
    //! 
    void
    build(
      SplineType1D      tp,
      real_type const * x,
      real_type const * y,
      integer           n
    ) {
      this->build( tp, x, 1, y, 1, n );
    }

    //! 
    //! Build a spline.
    //! 
    //! \param tp spline type
    //! \param x  vector of x-coordinates
    //! \param y  vector of y-coordinates
    //! 
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
    clear() { m_pSpline->clear(); }

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
    void
    dump(
      ostream_type & s,
      integer        nintervals,
      char const *   header = "x\ty"
    ) const {
      m_pSpline->dump( s, nintervals, header );
    }

    void
    dump(
      char const * fname,
      integer      nintervals,
      char const * header = "x\ty"
    ) const {
      m_pSpline->dump( fname, nintervals, header );
    }

    ///////////////////////////////////////////////////////////////////////////
    //! Evaluate spline value
    real_type operator () ( real_type x ) const { return (*m_pSpline)(x); }

    //! First derivative
    real_type D( real_type x ) const { return m_pSpline->D(x); }

    //! Second derivative
    real_type DD( real_type x ) const { return m_pSpline->DD(x); }

    //! Third derivative
    real_type DDD( real_type x ) const { return m_pSpline->DDD(x); }

    //! 4th derivative
    real_type DDDD( real_type x ) const { return m_pSpline->DDDD(x); }

    //! 5th derivative
    real_type DDDDD( real_type x ) const { return m_pSpline->DDDDD(x); }

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
      real_type * const cfs,
      real_type * const nodes,
      bool              transpose = false
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

    string
    info() const
    { return m_pSpline->info(); }

    void
    info( ostream_type & stream ) const
    { m_pSpline->info( stream ); }
  };

}

// EOF: Spline1D.hxx
