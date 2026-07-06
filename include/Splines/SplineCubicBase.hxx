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
 |    ____      _     _        ____        _ _              ____
 |   / ___|   _| |__ (_) ___  / ___| _ __ | (_)_ __   ___  | __ )  __ _ ___  ___
 |  | |  | | | | '_ \| |/ __| \___ \| '_ \| | | '_ \ / _ \ |  _ \ / _` / __|/ _ \
 |  | |__| |_| | |_) | | (__   ___) | |_) | | | | | |  __/ | |_) | (_| \__ \  __/
 |   \____\__,_|_.__/|_|\___| |____/| .__/|_|_|_| |_|\___| |____/ \__,_|___/\___|
 |                                  |_|
\*/

namespace Splines
{

  //!
  //! cubic spline base class
  //!
  class CubicSplineBase : public Spline
  {
  protected:
    Malloc_real m_mem_cubic;
    real_type * m_Yp             = nullptr;
    bool        m_external_alloc = false;

  public:
    using Spline::build;

    //!
    //! \name Contructors/Destructors
    ///@{
    //!
    //! Spline constructor.
    //!
    explicit CubicSplineBase( string_view name = "CubicSplineBase" )
      : Spline( name ), m_mem_cubic( fmt::format( "CubicSplineBase[{}]::m_mem_cubic", name ) )
    {
    }

    ~CubicSplineBase() override {}
    ///@}

    //!
    //! Build a copy of spline `S`
    //!
    void copy_spline( CubicSplineBase const & S );

    //!
    //! Return the pointer of values of yp-nodes.
    //!
    [[nodiscard]] real_type const * yp_nodes() const noexcept { return m_Yp; }

    //!
    //! Return the i-th node of the spline (y' component).
    //!
    [[nodiscard]] real_type yp_node( integer i ) const noexcept { return m_Yp[i]; }

    //!
    //! Change X-range of the spline.
    //!
    void set_range( real_type xmin, real_type xmax );

    //!
    //! Use externally allocated memory for `npts` points.
    //!
    void reserve_external( integer n, real_type *& p_x, real_type *& p_y, real_type *& p_dy )
    {
      m_npts_reserved  = n;
      m_X              = p_x;
      m_Y              = p_y;
      m_Yp             = p_dy;
      m_external_alloc = true;
      m_npts           = 0;
    }

    void y_min_max(
      integer &   i_min_pos,
      real_type & x_min_pos,
      real_type & y_min,
      integer &   i_max_pos,
      real_type & x_max_pos,
      real_type & y_max ) const override;

    void y_min_max(
      vector<integer> &   i_min_pos,
      vector<real_type> & x_min_pos,
      vector<real_type> & y_min,
      vector<integer> &   i_max_pos,
      vector<real_type> & x_max_pos,
      vector<real_type> & y_max ) const override;

    // --------------------------- VIRTUALS -----------------------------------

    //!
    //! \name Evaluation
    //!
    ///@{

    [[nodiscard]] real_type eval( real_type const x ) const override
    {
      std::pair<integer, real_type> res( 0, x );
      m_search.find( res );
      return this->id_eval( res.first, res.second );
    }

    [[nodiscard]] real_type D( real_type const x ) const override
    {
      std::pair<integer, real_type> res( 0, x );
      m_search.find( res );
      return this->id_D( res.first, res.second );
    }

    [[nodiscard]] real_type DD( real_type const x ) const override
    {
      std::pair<integer, real_type> res( 0, x );
      m_search.find( res );
      return this->id_DD( res.first, res.second );
    }

    [[nodiscard]] real_type DDD( real_type const x ) const override
    {
      std::pair<integer, real_type> res( 0, x );
      m_search.find( res );
      return this->id_DDD( res.first, res.second );
    }

    [[nodiscard]] real_type DDDD( real_type const ) const override { return 0; }
    [[nodiscard]] real_type DDDDD( real_type const ) const override { return 0; }

    void D( real_type const x, real_type dd[2] ) const override;

    void DD( real_type const x, real_type dd[3] ) const override;

    ///@}

    //!
    //! \name Evaluation when segment is known
    ///@{
    [[nodiscard]] real_type id_eval( integer const ni, real_type const x ) const override;

    [[nodiscard]] real_type id_D( integer const ni, real_type const x ) const override;

    [[nodiscard]] real_type id_DD( integer const ni, real_type const x ) const override;

    [[nodiscard]] real_type id_DDD( integer const ni, real_type const x ) const override;

    [[nodiscard]] real_type id_DDDD( integer const, real_type const ) const override { return 0; }
    [[nodiscard]] real_type id_DDDDD( integer const, real_type const ) const override { return 0; }
    ///@}

#ifdef AUTODIFF_SUPPORT
    //!
    //! \name Autodiff
    //!
    ///@{
    [[nodiscard]] autodiff::dual1st eval( autodiff::dual1st const & x ) const override;
    [[nodiscard]] autodiff::dual2nd eval( autodiff::dual2nd const & x ) const override;
    ///@}
#endif

    void write_to_stream( ostream_type & s ) const override;
    void export_csv( ostream_type & s ) const;
    void export_csv( string_view fname ) const
    {
      std::ofstream file{ std::string{ fname } };
      this->export_csv( file );
    }
    void export_json( ostream_type & s ) const;
    void export_json( string_view fname ) const
    {
      std::ofstream file{ std::string{ fname } };
      this->export_json( file );
    }
    void export_yaml( ostream_type & s ) const;
    void export_yaml( string_view fname ) const
    {
      std::ofstream file{ std::string{ fname } };
      this->export_yaml( file );
    }

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
    void build(
      real_type const x[],
      integer const   incx,
      real_type const y[],
      integer const   incy,
      real_type const yp[],
      integer const   incyp,
      integer         n );

    //!
    //! Build a spline.
    //!
    //! \param[in] x  vector of x-coordinates
    //! \param[in] y  vector of y-coordinates
    //! \param[in] yp vector of y'-coordinates
    //! \param[in] n  total number of points
    //!
    inline void build( real_type const x[], real_type const y[], real_type const yp[], integer const n )
    { this->build( x, 1, y, 1, yp, 1, n ); }

    //!
    //! Build a spline.
    //!
    //! \param[in] x  vector of x-coordinates
    //! \param[in] y  vector of y-coordinates
    //! \param[in] yp vector of y'-coordinates
    //!
    void build( vector<real_type> const & x, vector<real_type> const & y, vector<real_type> const & yp )
    {
      integer N = static_cast<integer>( x.size() );
      if ( N > static_cast<integer>( y.size() ) ) N = static_cast<integer>( y.size() );
      if ( N > static_cast<integer>( yp.size() ) ) N = static_cast<integer>( yp.size() );
      this->build( x.data(), 1, y.data(), 1, yp.data(), 1, N );
    }

    void reserve( integer npts ) override;

    ///@}

    void clear() override;

    integer  // order
    coeffs( real_type cfs[], real_type nodes[], bool transpose = false ) const override;

    [[nodiscard]] integer order() const override { return 4; }

    [[nodiscard]] bool is_monotone() const { return check_cubic_spline_monotonicity( m_X, m_Y, m_Yp, m_npts ); }

#ifdef SPLINES_BACK_COMPATIBILITY
    void                       copySpline( CubicSplineBase const & S ) { this->copy_spline( S ); }
    [[nodiscard]] integer      numPoints() const noexcept { return m_npts; }
    [[nodiscard]] real_type    xNode( integer i ) const noexcept { return m_X[i]; }
    [[nodiscard]] real_type    yNode( integer i ) const noexcept { return m_Y[i]; }
    [[nodiscard]] real_type    ypNode( integer i ) const noexcept { return this->yp_node( i ); }
    [[nodiscard]] real_type    xBegin() const noexcept { return m_X[0]; }
    [[nodiscard]] real_type    yBegin() const noexcept { return m_Y[0]; }
    [[nodiscard]] real_type    xEnd() const noexcept { return m_X[m_npts - 1]; }
    [[nodiscard]] real_type    yEnd() const noexcept { return m_Y[m_npts - 1]; }
    [[nodiscard]] real_type    xMin() const noexcept { return m_X[0]; }
    [[nodiscard]] real_type    xMax() const noexcept { return m_X[m_npts - 1]; }
    [[nodiscard]] real_type    yMin() const { return y_min(); }
    [[nodiscard]] real_type    yMax() const { return y_max(); }
#endif
  };

}  // namespace Splines
