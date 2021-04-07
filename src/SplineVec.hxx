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
 |   ____        _ _          __     __
 |  / ___| _ __ | (_)_ __   __\ \   / /__  ___
 |  \___ \| '_ \| | | '_ \ / _ \ \ / / _ \/ __|
 |   ___) | |_) | | | | | |  __/\ V /  __/ (__
 |  |____/| .__/|_|_|_| |_|\___| \_/ \___|\___|
 |        |_|
\*/

namespace Splines {

  //!
  //! Splines Management Class
  //!
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
    void allocate( integer dim, integer npts );
    void computeChords();

  public:

    //! \name Constructors
    ///@{

    //! spline constructor
    SplineVec( string const & name = "SplineVec" );

    //! spline destructor
    virtual
    ~SplineVec();

    ///@}

    //!
    //! Search the segment containing `x`
    //!
    integer search( real_type & x ) const;

    //!
    //! Spline name usd in the constructor
    //!
    string const & name() const { return m_name; }

    //! \name Open/Close
    ///@{
    bool is_closed() const { return m_curve_is_closed; }
    void make_closed()     { m_curve_is_closed = true; }
    void make_open()       { m_curve_is_closed = false; }

    bool can_extend() const { return m_curve_can_extend; }
    void make_unbounded()   { m_curve_can_extend = true; }
    void make_buonded()     { m_curve_can_extend = false; }
    ///@}

    //! \name Info
    ///@{

    //!
    //! return the number of support points of the splines
    //!
    integer numPoints() const { return m_npts; }

    //!
    //! return the number splines in the spline set
    //!
    integer dimension() const { return m_dim; }

    //!
    //! return the vector of values of x-nodes
    //!
    real_type const * xNodes() const { return m_X; }

    //!
    //! return the npt-th node of the spline (x component).
    //!
    real_type xNode( integer npt ) const { return m_X[size_t(npt)]; }

    //!
    //! return the npt-th node of the spline (y component).
    //!
    real_type
    yNode( integer npt, integer j ) const
    { return m_Y[size_t(j)][size_t(npt)]; }

    //!
    //! return x-minumum spline value
    //!
    real_type xMin() const { return m_X[0]; }

    //!
    //! return x-maximum spline value
    //!
    real_type xMax() const { return m_X[size_t(m_npts-1)]; }

    ///@}

    //! \name Evaluate
    ///@{

    //!
    //! Evaluate spline value
    //!
    real_type
    operator () ( real_type x, integer i ) const;

    //!
    //! Evaluate spline value at segment `i`
    //!
    real_type
    eval( real_type x, integer i ) const
    { return operator() (x,i); }

    //!
    //! First derivative
    //!
    real_type
    D( real_type x, integer i ) const;

    real_type
    eval_D( real_type x, integer i ) const
    { return this->D(x,i); }

    //!
    //! Second derivative
    //!
    real_type
    DD( real_type x, integer i ) const;

    real_type
    eval_DD( real_type x, integer i ) const
    { return this->DD(x,i); }

    //!
    //! Third derivative
    //!
    real_type
    DDD( real_type x, integer i ) const;

    real_type
    eval_DDD( real_type x, integer i ) const
    { return this->DDD(x,i); }

    //!
    //! 4th derivative
    //!
    real_type
    DDDD( real_type x, integer i ) const;

    real_type
    eval_DDDD( real_type x, integer i ) const
    { return this->DDDD(x,i); }

    //!
    //! 5th derivative
    //!
    real_type
    DDDDD( real_type x, integer i ) const;

    real_type
    eval_DDDDD( real_type x, integer i ) const
    { return this->DDDDD(x,i); }

    //!
    //! Evaluate all the splines at `x`
    //!
    void
    eval(
      real_type         x,
      real_type * const vals,
      integer           inc
    ) const;

    //!
    //! Evaluate the fist derivative of all the splines at `x`
    //!
    void
    eval_D(
      real_type         x,
      real_type * const vals,
      integer           inc
    ) const;

    //!
    //! Evaluate the second derivative of all the splines at `x`
    //!
    void
    eval_DD(
      real_type         x,
      real_type * const vals,
      integer           inc
    ) const;

    //!
    //! Evaluate the third derivative of all the splines at `x`
    //!
    void
    eval_DDD(
      real_type         x,
      real_type * const vals,
      integer           inc
    ) const;

    //!
    //! Evaluate the 4th derivative of all the splines at `x`
    //!
    void
    eval_DDDD(
      real_type         x,
      real_type * const vals,
      integer           inc
    ) const;

    //!
    //! Evaluate the 5th derivative of all the splines at `x`
    //!
    void
    eval_DDDDD(
      real_type         x,
      real_type * const vals,
      integer           inc
    ) const;

    //!
    //! Evaluate all the splines at `x`
    //!
    void eval( real_type x, vector<real_type> & vals ) const;

    //!
    //! Evaluate the fist derivative of all the splines at `x`
    //!
    void eval_D( real_type x, vector<real_type> & vals ) const;

    //!
    //! Evaluate the second derivative of all the splines at `x`
    //!
    void eval_DD( real_type x, vector<real_type> & vals ) const;

    //!
    //! Evaluate the third derivative of all the splines at `x`
    //!
    void eval_DDD( real_type x, vector<real_type> & vals ) const;

    //!
    //! Evaluate the 4th derivative of all the splines at `x`
    //!
    void eval_DDDD( real_type x, vector<real_type> & vals ) const;

    //!
    //! Evaluate the 5th derivative of all the splines at `x`
    //!
    void eval_DDDDD( real_type x, vector<real_type> & vals ) const;

    //!
    //! Evaluate at `x` and fill a GenericContainer
    //!
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

    //!
    //! Evaluate at `x` (x is a vector with many values) and fill a GenericContainer
    //!
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
    ///@}

    //! \name Setup Spline
    ///@{

    //!
    //! Initialize the interpolation point of the splines
    //! \param[in] dim  the dimension of the points
    //! \param[in] npts the numeber of interpolation points
    //! \param[in] Y    the matrix of points values, `Y[i]` is the 
    //!                 pointer to the i-th components 
    void
    setup(
      integer            dim,
      integer            npts,
      real_type const ** Y
    );

    //!
    //! Initialize the interpolation point of the splines
    //! \param[in] dim  the dimension of the points
    //! \param[in] npts the numeber of interpolation points
    //! \param[in] Y    the matrix of points values
    //! \param[in] ldY  leading dimension of the matrix
    //!
    void
    setup(
      integer           dim,
      integer           npts,
      real_type const * Y,
      integer           ldY
    );

    //!
    //! Set the knots of the spline
    //!
    void setKnots( real_type const * X );

    //!
    //! Computes the knots of the spline using chordal length
    //!
    void setKnotsChordLength();

    //!
    //! Computes the knots of the spline using centripetal approach
    //!
    void setKnotsCentripetal();

    //!
    //! Computes the knots of the spline using Foley algorithm
    //!
    void setKnotsFoley();

    //!
    //! Computes the spline using Catmull Rom algorithm.
    //! Points and nodes must be previously stored.
    //!
    void CatmullRom();

    void
    build( GenericContainer const & gc )
    { setup(gc); }

    ///@}

    //!
    //! Compute spline curvature at `x`
    //!
    real_type curvature( real_type x ) const;

    //!
    //! Compute spline curvature derivative at `x`
    //!
    real_type curvature_D( real_type x ) const;

    void setup( GenericContainer const & gc );

    //! Return spline type (as number)
    unsigned type() const { return SPLINE_VEC_TYPE; }

    string info() const;

    void
    info( ostream_type & stream ) const
    { stream << this->info() << '\n'; }

    void
    dump_table( ostream_type & s, integer num_points ) const;

  };

}

// EOF: SplineVec.hxx
