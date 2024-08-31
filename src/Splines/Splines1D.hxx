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

    #ifndef DOXYGEN_SHOULD_SKIP_THIS

    std::string m_name;
    Spline * m_pSpline{nullptr};

    #endif

  public:

    Spline1D( Spline1D const & ) = delete;
    Spline1D const & operator = ( Spline1D const & ) = delete;

    //! \name Constructors
    ///@{

    //!
    //! Build an empty spline of `Spline1D` type
    //!
    //! \param n the name of the spline
    //!
    Spline1D( std::string const & n )
    : m_name(n)
    {}

    //!
    //! Spline destructor.
    //!
    ~Spline1D()
    {}

    ///@}


    //!
    //! \name Open/Close
    //!
    ///@{

    //!
    //! \return string with the name of the spline
    //!
    string const & name() const { return m_pSpline->name(); }

    //! \return `true` if spline is a closed spline
    bool is_closed() const { return m_pSpline->is_closed(); }

    //!
    //! Set spline as a closed spline.
    //! When evaluated if parameter is outside the domain
    //! is wrapped cyclically before evalation.
    //!
    void make_closed() { m_pSpline->make_closed(); }

    //!
    //! Set spline as an opened spline.
    //! When evaluated if parameter is outside the domain
    //! an error is produced.
    //!
    void make_opened() { m_pSpline->make_opened(); }

    //!
    //! \return `true` if spline cannot extend outside interval of definition
    //!
    bool is_bounded() const { return m_pSpline->is_bounded(); }

    //!
    //! Set spline as unbounded.
    //! When evaluated if parameter is outside the domain
    //! an extrapolated value is used.
    //!
    void make_unbounded() { m_pSpline->make_unbounded();  }
    //!
    //! Set spline as bounded.
    //! When evaluated if parameter is outside the domain
    //! an error is issued.
    //!
    void make_bounded() { m_pSpline->make_bounded(); }

    //!
    //! \return `true` if the spline extend with a constant value
    //!
    bool is_extended_constant() const { return m_pSpline->is_extended_constant(); }
    //!
    //! Set spline to extend constant.
    //! When evaluated if parameter is outside the domain
    //! the value returned is the value of the closed border.
    //!
    void make_extended_constant() { m_pSpline->make_extended_constant(); }
    //!
    //! Set spline to extend NOT constant.
    //! When evaluated if parameter is outside the domain
    //! teh value returned is extrapolated using the last spline polynomial.
    //!
    void make_extended_not_constant() { m_pSpline->make_extended_not_constant(); }

    ///@}


    //!
    //! Return the number of support points of the spline.
    //!
    integer num_points() const { return m_pSpline->num_points(); }

    //!
    //! Return the i-th node of the spline (x component).
    //!
    real_type x_node( integer i ) const { return m_pSpline->x_node(i); }

    //!
    //! Return the i-th node of the spline (y component).
    //!
    real_type y_node( integer i ) const { return m_pSpline->y_node(i); }

    //!
    //! Return first node of the spline (x component).
    //!
    real_type x_begin() const { return m_pSpline->x_begin(); }

    //!
    //! Return first node of the spline (y component).
    //!
    real_type y_begin() const { return m_pSpline->y_begin(); }

    //!
    //! Return last node of the spline (x component).
    //!
    real_type x_end() const { return m_pSpline->x_end(); }

    //!
    //! Return last node of the spline (y component).
    //!
    real_type y_end() const { return m_pSpline->y_end(); }

    //!
    //! Allocate memory for `npts` points.
    //!
    void reserve( integer npts ) { return m_pSpline->reserve( npts ); }

    //!
    //! Add a support point (x,y) to the spline.
    //!
    void push_back( real_type x, real_type y ) { return m_pSpline->push_back( x, y ); }
    //!
    //! Drop last inserted point of the spline.
    //!
    void drop_back() { m_pSpline->drop_back(); }

    //!
    //! Build a spline.
    //!
    // must be defined in derived classes
    void build() { m_pSpline->build(); }

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
    //! \param tp   spline type
    //! \param x    vector of `x`-coordinates
    //! \param incx access elements as `x[0]`, `x[incx]`, `x[2*incx]`,...
    //! \param y    vector of `y`-coordinates
    //! \param incy access elements as `y[0]`, `y[incy]`, `x[2*incy]`,...
    //! \param n    total number of points
    //!
    // must be defined in derived classes
    void
    build(
      SplineType1D tp,
      real_type const x[], integer incx,
      real_type const y[], integer incy,
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
      SplineType1D    tp,
      real_type const x[],
      real_type const y[],
      integer         n
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
      integer n{integer(x.size())};
      this->build( tp, &x.front(), &y.front(), n );
    }

    //!
    //! Cancel the support points, empty the spline.
    //!
    void clear() { m_pSpline->clear(); }

    //!
    //! Return x-minumum spline value.
    //!
    real_type x_min() const { return m_pSpline->x_min(); }

    //!
    //! Return x-maximum spline value.
    //!
    real_type x_max() const { return m_pSpline->x_max(); }

    //!
    //! Return y-minumum spline value (on the support point of the spline).
    //!
    real_type y_min() const { return m_pSpline->y_min(); }

    //!
    //! Return y-maximum spline value (on the support point of the spline).
    //!
    real_type y_max() const { return m_pSpline->y_max(); }

    ///////////////////////////////////////////////////////////////////////////
    //!
    //! Change X-origin of the spline.
    //!
    void set_origin( real_type x0 ) { return m_pSpline->set_origin( x0 ); }

    //!
    //! Change X-range of the spline.
    //!
    void set_range( real_type xmin, real_type xmax ) { return m_pSpline->set_range( xmin, xmax ); }

    ///////////////////////////////////////////////////////////////////////////
    //!
    //! dump a sample of the spline
    //!
    void
    dump(
      ostream_type & s,
      integer        nintervals,
      char const     header[] = "x\ty"
    ) const {
      m_pSpline->dump( s, nintervals, header );
    }

    //!
    //! dump a sample of the spline
    //!
    void
    dump(
      char const fname[],
      integer    nintervals,
      char const header[] = "x\ty"
    ) const {
      m_pSpline->dump( fname, nintervals, header );
    }

    ///////////////////////////////////////////////////////////////////////////
    //! \name Evaluation
    ///@{

    //!
    //! Evaluate spline value at `x`.
    //!
    real_type eval( real_type x ) const { return m_pSpline->eval(x); }

    //!
    //! First derivative at `x`.
    //!
    real_type D( real_type x ) const { return m_pSpline->D(x); }

    //!
    //! Second derivative at `x`.
    //!
    real_type DD( real_type x ) const { return m_pSpline->DD(x); }

    //!
    //! Third derivative at `x`.
    //!
    real_type DDD( real_type x ) const { return m_pSpline->DDD(x); }

    //!
    //! 4th derivative at `x`.
    //!
    real_type DDDD( real_type x ) const { return m_pSpline->DDDD(x); }

    //!
    //! 5th derivative at `x`.
    //!
    real_type DDDDD( real_type x ) const { return m_pSpline->DDDDD(x); }

    ///@}

    //!
    //! \name Evaluation Aliases
    //!
    ///@{
    //! the value of the spline at `x`
    real_type operator () ( real_type x ) const { return m_pSpline->eval(x); }

    //! the value of the first derivative of the spline at `x`
    real_type eval_D( real_type x ) const { return m_pSpline->D(x); }

    //! the value of the second derivative of the spline at `x`
    real_type eval_DD( real_type x ) const { return m_pSpline->DD(x); }

    //! the value of the third derivative of the spline at `x`
    real_type eval_DDD( real_type x ) const { return m_pSpline->DDD(x); }

    //! the value of the 4-th derivative of the spline at `x`
    real_type eval_DDDD( real_type x ) const { return m_pSpline->DDDD(x); }

    //! the value of the 5-th derivative of the spline at `x`
    real_type eval_DDDDD( real_type x ) const { return m_pSpline->DDDDD(x); }
    ///@}

    ///////////////////////////////////////////////////////////////////////////
    //!
    //! \name Evaluation when segment is known
    ///@{
    //!
    //! Evaluate spline at `x`.
    //! \param x  value at which spline is evaluated
    //! \param ni select the component
    //!
    real_type
    id_eval( integer ni, real_type x ) const { return m_pSpline->id_eval(ni,x); }

    //!
    //! First derivative at `x`.
    //! \param x  value at which spline is evaluated
    //! \param ni select the component
    //!
    real_type
    id_D( integer ni, real_type x ) const { return m_pSpline->id_D(ni,x); }

    //!
    //! Second derivative at `x`.
    //! \param x  value at which spline is evaluated
    //! \param ni select the component
    //!
    real_type
    id_DD( integer ni, real_type x ) const { return m_pSpline->id_DD(ni,x); }

    //!
    //! Third derivative at `x`.
    //! \param x  value at which spline is evaluated
    //! \param ni select the component
    //!
    real_type
    id_DDD( integer ni, real_type x ) const { return m_pSpline->id_DDD(ni,x); }

    //!
    //! 4th derivative at `x`.
    //! \param x  value at which spline is evaluated
    //! \param ni select the component
    //!
    real_type
    id_DDDD( integer ni, real_type x ) const { return m_pSpline->id_DDDD(ni,x); }

    //!
    //! 5th derivative at `x`.
    //! \param x  value at which spline is evaluated
    //! \param ni select the component
    //!
    real_type
    id_DDDDD( integer ni, real_type x ) const { return m_pSpline->id_DDDDD(ni,x); }

    ///@}

    //!
    //! Get the piecewise polinomials of the spline
    //!
    integer // order
    coeffs(
      real_type cfs[],
      real_type nodes[],
      bool      transpose = false
    ) const {
      return m_pSpline->coeffs( cfs, nodes, transpose );
    }

    //!
    //! \return the order of the spline (`degree+1`)
    //!
    integer order() const { return m_pSpline->order(); }

    //!
    //! Print spline coefficients.
    //!
    void
    write_to_stream( ostream_type & s ) const
    { return m_pSpline->write_to_stream( s ); }

    //!
    //! Return spline typename.
    //!
    char const *
    type_name() const
    { return m_pSpline->type_name(); }

    //!
    //! Return spline type (as number).
    //!
    SplineType1D
    type() const
    { return m_pSpline->type(); }

    //!
    //! String information of the kind and order of the spline
    //!
    string
    info() const
    { return m_pSpline->info(); }

    //!
    //! Print information of the kind and order of the spline
    //!
    void
    info( ostream_type & stream ) const
    { m_pSpline->info( stream ); }

    #ifdef SPLINES_BACK_COMPATIBILITY
    integer numPoints() const { return m_pSpline->num_points(); }
    real_type xNode( integer i ) const { return this->x_node(i); }
    real_type yNode( integer i ) const { return this->y_node(i); }
    real_type xBegin() const { return this->x_begin(); }
    real_type yBegin() const { return this->y_begin(); }
    real_type xEnd() const { return this->x_end(); }
    real_type yEnd() const { return this->y_end(); }
    void pushBack( real_type x, real_type y ) { return this->push_back( x, y ); }
    void dropBack() { this->drop_back(); }
    real_type xMin() const { return this->x_min(); }
    real_type xMax() const { return this->x_max(); }
    real_type yMin() const { return this->y_min(); }
    real_type yMax() const { return this->y_max(); }
    void setOrigin( real_type x0 ) { return this->set_origin( x0 ); }
    void setRange( real_type xmin, real_type xmax ) { return this->set_range( xmin, xmax ); }
    #endif

  };

}

// EOF: Spline1D.hxx
