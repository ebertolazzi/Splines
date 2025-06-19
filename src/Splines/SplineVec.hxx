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

  protected:

    string const m_name;

    Utils::Malloc<real_type>  m_mem;
    Utils::Malloc<real_type*> m_mem_p;

    integer m_dim{0};
    integer m_npts{0};
    bool    m_curve_is_closed{false};
    bool    m_curve_can_extend{true};

    real_type *  m_X{nullptr};
    real_type ** m_Y{nullptr};
    real_type ** m_Yp{nullptr};

    SearchInterval m_search;

    void allocate( integer dim, integer npts );
    void compute_chords();

  public:

    SplineVec( SplineVec const & ) = delete;
    SplineVec const & operator = ( SplineVec const & ) = delete;

    //!
    //! \name Constructors
    //!
    ///@{

    //!
    //! Build an empty spline of `SplineVec` type
    //!
    //! \param name the name of the spline
    //!
    explicit
    SplineVec( string_view name = "SplineVec" );

    //!
    //! Spline destructor.
    //!
    virtual
    ~SplineVec();

    ///@}

    //!
    //! Spline name usd in the constructor
    //!
    string_view name() const { return m_name; }

    //!
    //! \name Open/Close
    //!
    ///@{

    //! \return `true` if spline is a closed spline
    bool is_closed() const { return m_curve_is_closed; }

    //!
    //! Set spline as a closed spline.
    //! When evaluated if parameter is outside the domain
    //! is wrapped cyclically before evalation.
    //!
    void make_closed() { m_curve_is_closed = true; }

    //!
    //! Set spline as an opened spline.
    //! When evaluated if parameter is outside the domain
    //! an error is produced.
    //!
    void make_open() { m_curve_is_closed = false; }

    //!
    //! \return `true` if spline can extend outside interval of definition
    //!
    bool can_extend() const { return m_curve_can_extend; }

    //!
    //! Set spline as unbounded.
    //! When evaluated if parameter is outside the domain
    //! an extrapolated value is used.
    //!
    void make_unbounded() { m_curve_can_extend = true; }

    //!
    //! Set spline as bounded.
    //! When evaluated if parameter is outside the domain
    //! an error is issued.
    //!
    void make_buonded() { m_curve_can_extend = false; }
    ///@}

    //!
    //! \name Info
    //!
    ///@{

    //!
    //! Return the number of support points of the splines.
    //!
    integer num_points() const { return m_npts; }

    //!
    //! Return the number splines in the spline set.
    //!
    integer dimension() const { return m_dim; }

    //!
    //! Return the vector of values of x-nodes.
    //!
    real_type const * x_nodes() const { return m_X; }

    //!
    //! Return the npt-th node of the spline (x component).
    //!
    real_type x_node( integer npt ) const { return m_X[npt]; }

    //!
    //! Return the vector of values of y-nodes, component `j`
    //!
    real_type const * y_nodes( integer j ) const { return m_Y[j]; }

    //!
    //! Return the npt-th node of the spline (`j` component of y).
    //!
    real_type
    y_node( integer npt, integer j ) const
    { return m_Y[j][npt]; }

    //!
    //! Return x-minumum spline value.
    //!
    real_type x_min() const { return m_X[0]; }

    //!
    //! Return x-maximum spline value.
    //!
    real_type x_max() const { return m_X[m_npts-1]; }

    ///@}

    //!
    //! \name Evaluation when segment is known.
    //!
    ///@{

    //!
    //! Evaluate spline value at `x` component `i`-th.
    //!
    real_type
    eval( real_type const x, integer const i ) const;

    //!
    //! Evaluate spline value at `x` component `i`-th.
    //!
    real_type
    operator () ( real_type const x, integer const i ) const
    { return this->eval( x, i ); }

    //!
    //! First derivative value at `x` component `i`-th.
    //!
    real_type
    D( real_type const x, integer const i ) const;

    //!
    //! First derivative value at `x` component `i`-th.
    //!
    real_type
    eval_D( real_type const x, integer const i ) const
    { return this->D(x,i); }

    //!
    //! Second derivative value at `x` component `i`-th.
    //!
    real_type
    DD( real_type const x, integer const i ) const;

    //!
    //! Second derivative value at `x` component `i`-th.
    //!
    real_type
    eval_DD( real_type const x, integer const i ) const
    { return this->DD(x,i); }

    //!
    //! Third derivative value at `x` component `i`-th.
    //!
    real_type
    DDD( real_type const x, integer const i ) const;

    //!
    //! Third derivative value at `x` component `i`-th.
    //!
    real_type
    eval_DDD( real_type const x, integer const i ) const
    { return this->DDD(x,i); }

    //!
    //! 4th derivative value at `x` component `i`-th.
    //!
    real_type
    DDDD( real_type const x, integer const i ) const;

    //!
    //! 4th derivative value at `x` component `i`-th.
    //!
    real_type
    eval_DDDD( real_type const x, integer const i ) const
    { return this->DDDD(x,i); }

    //!
    //! 5th derivative value at `x` component `i`-th.
    //!
    real_type
    DDDDD( real_type x, integer i ) const;

    //!
    //! 5th derivative value at `x` component `i`-th.
    //!
    real_type
    eval_DDDDD( real_type const x, integer const i ) const
    { return this->DDDDD(x,i); }

    ///@}

    #ifdef AUTIDIFF_SUPPORT
    //!
    //! \name Autodiff
    //!
    ///@{
    autodiff::dual1st
    eval( autodiff::dual1st const & x, integer const i ) const {
      using autodiff::dual1st;
      using autodiff::derivative;
      real_type xv  { val(x) };
      dual1st   res { eval(xv,i) };
      res.grad = eval_D(xv,i) * x.grad;
      return res;
    }

    autodiff::dual2nd
    eval( autodiff::dual2nd const & x, integer const i ) const {
      using autodiff::dual2nd;
      using autodiff::derivative;

      real_type xv  { val(x) };
      real_type xg  { val(x.grad) };
      real_type dfx { eval_D(xv,i) };
      real_type dxx { eval_DD(xv,i) };
      dual2nd   res { eval(xv,i) };

      res.grad      = dfx * xg;
      res.grad.grad = dfx * x.grad.grad + dxx * (xg*xg);
      return res;
    }

    template <typename T>
    autodiff::HigherOrderDual<autodiff::detail::DualOrder<T>::value,real_type>
    eval( T const & x, integer const i ) const {
      autodiff::HigherOrderDual<autodiff::detail::DualOrder<T>::value,real_type> X{x};
      return eval( X, i );
    }
    ///@}
    #endif
    
    //!
    //! \name Evaluate all the splines in a vector.
    //!
    ///@{

    //!
    //! Evaluate all the splines at `x` and
    //! store values in `vals` with stride `inc`.
    //!
    void eval( real_type const x, real_type vals[], integer const inc ) const;

    //!
    //! Evaluate the fist derivative of all the splines at `x` and
    //! store values in `vals` with stride `inc`.
    //!
    void eval_D( real_type const x, real_type vals[], integer const inc ) const;

    //!
    //! Evaluate the second derivative of all the splines at `x` and
    //! store values in `vals` with stride `inc`.
    //!
    void eval_DD( real_type const x, real_type vals[], integer const inc ) const;

    //!
    //! Evaluate the third derivative of all the splines at `x` and
    //! store values in `vals` with stride `inc`.
    //!
    void eval_DDD( real_type const x, real_type vals[], integer const inc ) const;

    //!
    //! Evaluate the 4th derivative of all the splines at `x` and
    //! store values in `vals` with stride `inc`.
    //!
    void eval_DDDD( real_type const x, real_type vals[], integer const inc ) const;

    //!
    //! Evaluate the 5th derivative of all the splines at `x` and
    //! store values in `vals` with stride `inc`.
    //!
    void eval_DDDDD( real_type const x, real_type vals[], integer const inc ) const;
    ///@}

    //!
    //! \name Evaluate all the splines in an STL vector
    //!
    ///@{

    //!
    //! Evaluate all the splines at `x` and
    //! store values in `vals`.
    //!
    void eval( real_type const x, vector<real_type> & vals ) const;

    //!
    //! Evaluate the fist derivative of all the splines at `x` and
    //! store values in `vals`.
    //!
    void eval_D( real_type const x, vector<real_type> & vals ) const;

    //!
    //! Evaluate the second derivative of all the splines at `x` and
    //! store values in `vals`.
    //!
    void eval_DD( real_type const x, vector<real_type> & vals ) const;

    //!
    //! Evaluate the third derivative of all the splines at `x` and
    //! store values in `vals`.
    //!
    void eval_DDD( real_type const x, vector<real_type> & vals ) const;

    //!
    //! Evaluate the 4th derivative of all the splines at `x` and
    //! store values in `vals`.
    //!
    void eval_DDDD( real_type const x, vector<real_type> & vals ) const;

    //!
    //! Evaluate the 5th derivative of all the splines at `x` and
    //! store values in `vals`.
    //!
    void eval_DDDDD( real_type const x, vector<real_type> & vals ) const;
    ///@}

    //!
    //! \name Evaluate all the splines in a `GenericContainer`.
    //!
    ///@{

    //!
    //! Evaluate at `x` and fill a `GenericContainer`
    //!
    void
    eval( real_type const x, GenericContainer & vals ) const
    { eval( x, vals.set_vec_real(m_dim) ); }

    //!
    //! Evaluate first derivatives at `x` and fill a `GenericContainer`
    //!
    void
    eval_D( real_type const x, GenericContainer & vals ) const
    { eval_D( x, vals.set_vec_real(m_dim) ); }

    //!
    //! Evaluate second derivatives at `x` and fill a `GenericContainer`
    //!
    void
    eval_DD( real_type const x, GenericContainer & vals ) const
    { eval_DD( x, vals.set_vec_real(m_dim) ); }

    //!
    //! Evaluate third derivatives at `x` and fill a `GenericContainer`
    //!
    void
    eval_DDD( real_type const x, GenericContainer & vals ) const
    { eval_DDD( x, vals.set_vec_real(m_dim) ); }

    //!
    //! Evaluate 4th derivatives at `x` and fill a `GenericContainer`
    //!
    void
    eval_DDDD( real_type const x, GenericContainer & vals ) const
    { eval_DDDD( x, vals.set_vec_real(m_dim) ); }

    //!
    //! Evaluate 5th derivatives at `x` and fill a `GenericContainer`
    //!
    void
    eval_DDDDD( real_type const x, GenericContainer & vals ) const
    { eval_DDDDD( x, vals.set_vec_real(m_dim) ); }
    ///@}

    //!
    //! \name Evaluate all the splines in a point set in a `GenericContainer`
    //!
    ///@{

    //!
    //! Evaluate at `x` (x is a vector with many values)
    //! and fill a GenericContainer
    //!
    void
    eval( vec_real_type const & x, GenericContainer & vals ) const;

    //!
    //! Evaluate first derivative at `x` (x is a vector with many values)
    //! and fill a GenericContainer
    //!
    void
    eval_D( vec_real_type const & x, GenericContainer & vals ) const;

    //!
    //! Evaluate second derivative at `x` (x is a vector with many values)
    //! and fill a GenericContainer
    //!
    void
    eval_DD( vec_real_type const & x, GenericContainer & vals ) const;

    //!
    //! Evaluate third derivative at `x` (x is a vector with many values)
    //! and fill a GenericContainer
    //!
    void
    eval_DDD( vec_real_type const & x, GenericContainer & vals ) const;

    //!
    //! Evaluate 4th derivative at `x` (x is a vector with many values)
    //! and fill a GenericContainer
    //!
    void
    eval_DDDD( vec_real_type const & x, GenericContainer & vals ) const;

    //!
    //! Evaluate 5th derivative at `x` (x is a vector with many values)
    //! and fill a GenericContainer
    //!
    void
    eval_DDDDD( vec_real_type const & x, GenericContainer & vals ) const;
    ///@}

    //!
    //! \name Setup Splines.
    //!
    ///@{

    //!
    //! Initialize the interpolation point of the splines.
    //!
    //! \param[in] dim  the dimension of the points
    //! \param[in] npts the number of interpolation points
    //! \param[in] Y    the matrix of points values, `Y[i]` is the
    //!                 pointer to the i-th components
    //!
    void
    setup(
      integer                 dim,
      integer                 npts,
      real_type const * const Y[]
    );

    //!
    //! Initialize the interpolation point of the splines.
    //!
    //! \param[in] dim  the dimension of the points
    //! \param[in] npts the number of interpolation points
    //! \param[in] Y    the matrix of points values
    //! \param[in] ldY  leading dimension of the matrix
    //!
    void
    setup(
      integer         dim,
      integer         npts,
      real_type const Y[],
      integer         ldY
    );

    //!
    //! Copy to SplineVec `S`
    //!
    void deep_copy_to( SplineVec & S ) const;
  
    //!
    //! Set the knots of the spline.
    //!
    void set_knots( real_type const X[] );

    //!
    //! Computes the knots of the spline using chordal length.
    //!
    void set_knots_chord_length();

    //!
    //! Computes the knots of the spline using centripetal approach.
    //!
    void set_knots_centripetal();

    //!
    //! Computes the knots of the spline using Foley algorithm.
    //!
    void set_knots_foley();

    //!
    //! Computes the spline using Catmull Rom algorithm.
    //! Points and nodes must be previously stored.
    //!
    void catmull_rom();

    //!
    //! Build a spline using data in `GenericContainer`
    //!
    void
    build( GenericContainer const & gc )
    { setup( gc ); }

    //!
    //! Build a spline using data in `GenericContainer`
    //!
    void setup( GenericContainer const & gc );
  
    ///@}

    //!
    //! Compute spline curvature at `x`.
    //!
    real_type curvature( real_type x ) const;

    //!
    //! Compute spline curvature derivative at `x`.
    //!
    real_type curvature_D( real_type x ) const;

    //!
    //! Return spline type (as number).
    //!
    static SplineType1D type() { return SplineType1D::SPLINE_VEC; }

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

    //!
    //! Dump values of the spline on a stream for plotting
    //!
    void
    dump_table( ostream_type & s, integer num_points ) const;

    #ifdef SPLINES_BACK_COMPATIBILITY
    integer numPoints() const { return m_npts; }
    real_type const * xNodes() const { return m_X; }
    real_type xNode( integer npt ) const { return m_X[npt]; }
    real_type const * yNodes( integer j ) const { return m_Y[j]; }
    real_type yNode( integer npt, integer j ) const { return y_node(npt,j); }
    real_type xMin() const { return m_X[0]; }
    real_type xMax() const { return m_X[m_npts-1]; }
    void setKnots( real_type const X[] ) { this->set_knots( X ); }
    void setKnotsChordLength() { this->set_knots_chord_length(); }
    void setKnotsCentripetal() { this->set_knots_centripetal(); }
    void setKnotsFoley() { this->set_knots_foley(); }
    void CatmullRom() { this->catmull_rom(); }
    #endif

  };

}

// EOF: SplineVec.hxx
