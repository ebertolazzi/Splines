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
    std::unique_ptr<Spline> m_spline;

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
    Spline1D( string_view n )
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
    string_view name() const { return m_spline->name(); }

    //! \return `true` if spline is a closed spline
    bool is_closed() const { return m_spline->is_closed(); }

    //!
    //! Set spline as a closed spline.
    //! When evaluated if parameter is outside the domain
    //! is wrapped cyclically before evalation.
    //!
    void make_closed() { m_spline->make_closed(); }

    //!
    //! Set spline as an opened spline.
    //! When evaluated if parameter is outside the domain
    //! an error is produced.
    //!
    void make_opened() { m_spline->make_opened(); }

    //!
    //! \return `true` if spline cannot extend outside interval of definition
    //!
    bool is_bounded() const { return m_spline->is_bounded(); }

    //!
    //! Set spline as unbounded.
    //! When evaluated if parameter is outside the domain
    //! an extrapolated value is used.
    //!
    void make_unbounded() { m_spline->make_unbounded();  }
    //!
    //! Set spline as bounded.
    //! When evaluated if parameter is outside the domain
    //! an error is issued.
    //!
    void make_bounded() { m_spline->make_bounded(); }

    //!
    //! \return `true` if the spline extend with a constant value
    //!
    bool is_extended_constant() const { return m_spline->is_extended_constant(); }
    //!
    //! Set spline to extend constant.
    //! When evaluated if parameter is outside the domain
    //! the value returned is the value of the closed border.
    //!
    void make_extended_constant() { m_spline->make_extended_constant(); }
    //!
    //! Set spline to extend NOT constant.
    //! When evaluated if parameter is outside the domain
    //! teh value returned is extrapolated using the last spline polynomial.
    //!
    void make_extended_not_constant() { m_spline->make_extended_not_constant(); }

    ///@}


    //!
    //! Return the number of support points of the spline.
    //!
    integer num_points() const { return m_spline->num_points(); }

    //!
    //! Return the i-th node of the spline (x component).
    //!
    real_type x_node( integer i ) const { return m_spline->x_node(i); }

    //!
    //! Return the i-th node of the spline (y component).
    //!
    real_type y_node( integer i ) const { return m_spline->y_node(i); }

    //!
    //! Return first node of the spline (x component).
    //!
    real_type x_begin() const { return m_spline->x_begin(); }

    //!
    //! Return first node of the spline (y component).
    //!
    real_type y_begin() const { return m_spline->y_begin(); }

    //!
    //! Return last node of the spline (x component).
    //!
    real_type x_end() const { return m_spline->x_end(); }

    //!
    //! Return last node of the spline (y component).
    //!
    real_type y_end() const { return m_spline->y_end(); }

    //!
    //! Allocate memory for `npts` points.
    //!
    void reserve( integer npts ) { return m_spline->reserve( npts ); }

    //!
    //! Add a support point (x,y) to the spline.
    //!
    void push_back( real_type x, real_type y ) { return m_spline->push_back( x, y ); }
    //!
    //! Drop last inserted point of the spline.
    //!
    void drop_back() { m_spline->drop_back(); }

    //!
    //! Build a spline.
    //!
    // must be defined in derived classes
    void build() { m_spline->build(); }

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
      this->build( tp, x.data(), y.data(), n );
    }

    //!
    //! Cancel the support points, empty the spline.
    //!
    void clear() { m_spline->clear(); }

    //!
    //! Return x-minumum spline value.
    //!
    real_type x_min() const { return m_spline->x_min(); }

    //!
    //! Return x-maximum spline value.
    //!
    real_type x_max() const { return m_spline->x_max(); }

    //!
    //! Return y-minumum spline value (on the support point of the spline).
    //!
    real_type y_min() const { return m_spline->y_min(); }

    //!
    //! Return y-maximum spline value (on the support point of the spline).
    //!
    real_type y_max() const { return m_spline->y_max(); }

    ///////////////////////////////////////////////////////////////////////////
    //!
    //! Change X-origin of the spline.
    //!
    void set_origin( real_type x0 ) const { return m_spline->set_origin( x0 ); }

    //!
    //! Change X-range of the spline.
    //!
    void set_range( real_type xmin, real_type xmax ) { return m_spline->set_range( xmin, xmax ); }

    ///////////////////////////////////////////////////////////////////////////
    //!
    //! dump a sample of the spline
    //!
    void
    dump(
      ostream_type & s,
      integer        nintervals,
      string_view    header = "x\ty"
    ) const {
      m_spline->dump( s, nintervals, header );
    }

    //!
    //! dump a sample of the spline
    //!
    void
    dump(
      string_view fname,
      integer     nintervals,
      string_view header = "x\ty"
    ) const {
      m_spline->dump( fname, nintervals, header );
    }

    ///////////////////////////////////////////////////////////////////////////
    //! \name Evaluation
    ///@{

    //!
    //! Evaluate spline value at `x`.
    //!
    real_type eval( real_type const x ) const { return m_spline->eval(x); }

    //!
    //! First derivative at `x`.
    //!
    real_type D( real_type const x ) const { return m_spline->D(x); }

    //!
    //! Second derivative at `x`.
    //!
    real_type DD( real_type const x ) const { return m_spline->DD(x); }

    //!
    //! Third derivative at `x`.
    //!
    real_type DDD( real_type const x ) const { return m_spline->DDD(x); }

    //!
    //! 4th derivative at `x`.
    //!
    real_type DDDD( real_type const x ) const { return m_spline->DDDD(x); }

    //!
    //! 5th derivative at `x`.
    //!
    real_type DDDDD( real_type const x ) const { return m_spline->DDDDD(x); }

    void D  ( real_type const x, real_type dd[2] ) const { m_spline->D(x,dd); }
    void DD ( real_type const x, real_type dd[2] ) const { m_spline->DD(x,dd); }

    ///@}

    #ifdef AUTIDIFF_SUPPORT
    //!
    //! \name Autodiff
    //!
    ///@{
    autodiff::dual1st eval( autodiff::dual1st const & x ) const;
    autodiff::dual2nd eval( autodiff::dual2nd const & x ) const;

    template <typename T>
    autodiff::HigherOrderDual<autodiff::detail::DualOrder<T>::value,real_type>
    eval( T const & x ) const {
      autodiff::HigherOrderDual<autodiff::detail::DualOrder<T>::value,real_type> X{x};
      return eval( X );
    }
    ///@}
    #endif

    //!
    //! \name Evaluation Aliases
    //!
    ///@{
    //! the value of the spline at `x`
    real_type operator () ( real_type const x ) const { return m_spline->eval(x); }

    //! the value of the first derivative of the spline at `x`
    real_type eval_D( real_type const x ) const { return m_spline->D(x); }

    //! the value of the second derivative of the spline at `x`
    real_type eval_DD( real_type const x ) const { return m_spline->DD(x); }

    //! the value of the third derivative of the spline at `x`
    real_type eval_DDD( real_type const x ) const { return m_spline->DDD(x); }

    //! the value of the 4-th derivative of the spline at `x`
    real_type eval_DDDD( real_type const x ) const { return m_spline->DDDD(x); }

    //! the value of the 5-th derivative of the spline at `x`
    real_type eval_DDDDD( real_type const x ) const { return m_spline->DDDDD(x); }
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
    real_type id_eval( integer const ni, real_type const x ) const { return m_spline->id_eval(ni,x); }

    //!
    //! First derivative at `x`.
    //! \param x  value at which spline is evaluated
    //! \param ni select the component
    //!
    real_type id_D( integer const ni, real_type const x ) const { return m_spline->id_D(ni,x); }

    //!
    //! Second derivative at `x`.
    //! \param x  value at which spline is evaluated
    //! \param ni select the component
    //!
    real_type id_DD( integer const ni, real_type const x ) const { return m_spline->id_DD(ni,x); }

    //!
    //! Third derivative at `x`.
    //! \param x  value at which spline is evaluated
    //! \param ni select the component
    //!
    real_type id_DDD( integer const ni, real_type const x ) const { return m_spline->id_DDD(ni,x); }

    //!
    //! 4th derivative at `x`.
    //! \param x  value at which spline is evaluated
    //! \param ni select the component
    //!
    real_type id_DDDD( integer const ni, real_type const x ) const { return m_spline->id_DDDD(ni,x); }

    //!
    //! 5th derivative at `x`.
    //! \param x  value at which spline is evaluated
    //! \param ni select the component
    //!
    real_type id_DDDDD( integer const ni, real_type const x ) const { return m_spline->id_DDDDD(ni,x); }

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
      return m_spline->coeffs( cfs, nodes, transpose );
    }

    //!
    //! \return the order of the spline (`degree+1`)
    //!
    integer order() const { return m_spline->order(); }

    //!
    //! Print spline coefficients.
    //!
    void
    write_to_stream( ostream_type & s ) const
    { return m_spline->write_to_stream( s ); }

    //!
    //! Return spline typename.
    //!
    char const * type_name() const { return m_spline->type_name(); }

    //!
    //! Return spline type (as number).
    //!
    SplineType1D type() const { return m_spline->type(); }

    //!
    //! String information of the kind and order of the spline
    //!
    string info() const { return m_spline->info(); }

    //!
    //! Print information of the kind and order of the spline
    //!
    void info( ostream_type & stream ) const { m_spline->info( stream ); }

    #ifdef SPLINES_BACK_COMPATIBILITY
    integer numPoints() const { return m_spline->num_points(); }
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
