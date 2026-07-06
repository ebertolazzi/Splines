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
 |    ___        _       _   _      ____        _ _            ____
 |   / _ \ _   _(_)_ __ | |_(_) ___/ ___| _ __ | (_)_ __   ___| __ )  __ _ ___  ___
 |  | | | | | | | | '_ \| __| |/ __\___ \| '_ \| | | '_ \ / _ \  _ \ / _` / __|/ _ \
 |  | |_| | |_| | | | | | |_| | (__ ___) | |_) | | | | | |  __/ |_) | (_| \__ \  __/
 |   \__\_\\__,_|_|_| |_|\__|_|\___|____/| .__/|_|_|_| |_|\___|____/ \__,_|___/\___|
 |                                       |_|
\*/

#ifndef SPLINE_QUINTIC_BASE_HH
#define SPLINE_QUINTIC_BASE_HH

namespace Splines
{

  //!
  //! Quintic spline base class
  //!
  class QuinticSplineBase : public Spline
  {
  protected:
    Malloc_real     m_base_quintic;
    real_type *     m_Yp             = nullptr;
    real_type *     m_Ypp            = nullptr;
    bool            m_external_alloc = false;
    Spline_sub_type m_sub_type;

  public:
    //!
    //! \name Constructors
    //!
    ///@{

    using Spline::build;

    //!
    //! Build an empty spline of `QuinticSplineBase` type
    //!
    //! \param name the name of the spline
    //!
    explicit QuinticSplineBase( Spline_sub_type type = Spline_sub_type::PCHIP, string_view name = "QuinticSplineBase" )
      : Spline( name ), m_base_quintic( fmt::format( "QuinticSplineBase[{}]", name ) ), m_sub_type( type )
    {
    }

    //!
    //! Spline destructor.
    //!
    ~QuinticSplineBase() override {}

    ///@}

    //!
    //! Build a copy of spline `S`
    //!
    void copy_spline( QuinticSplineBase const & S );

    //!
    //! \name Info
    //!
    ///@{

    //!
    //! Return the pointer of values of yp-nodes.
    //!
    [[nodiscard]] real_type const * yp_nodes() const noexcept { return m_Yp; }

    //!
    //! Return the pointer of values of ypp-nodes.
    //!
    [[nodiscard]] real_type const * ypp_nodes() const noexcept { return m_Ypp; }

    //!
    //! Return the i-th node of the spline (y' component).
    //!
    [[nodiscard]] real_type yp_node( integer const i ) const noexcept { return m_Yp[i]; }

    //!
    //! Return the i-th node of the spline (y'' component).
    //!
    [[nodiscard]] real_type ypp_node( integer const i ) const noexcept { return m_Ypp[i]; }

    void y_min_max(
      integer &   i_min_pos,
      real_type & x_min_pos,
      real_type & y_min,
      integer &   i_max_pos,
      real_type & x_max_pos,
      real_type & y_max ) const override;

    void y_min_max(
      std::vector<integer> &   i_min_pos,
      std::vector<real_type> & x_min_pos,
      std::vector<real_type> & y_min,
      std::vector<integer> &   i_max_pos,
      std::vector<real_type> & x_max_pos,
      std::vector<real_type> & y_max ) const override;

    void write_to_stream( std::ostream & s ) const override;
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

    [[nodiscard]] SplineType1D type() const override
    {
      switch ( m_sub_type )
      {
        case Spline_sub_type::CUBIC: return SplineType1D::QUINTIC_CUBIC;
        case Spline_sub_type::AKIMA: return SplineType1D::QUINTIC_AKIMA;
        case Spline_sub_type::VANLEER: return SplineType1D::QUINTIC_VANLEER;
        case Spline_sub_type::PCHIP: return SplineType1D::QUINTIC_PCHIP;
      }
    }

    ///@}

    //!
    //! Change X-range of the spline
    //!
    void set_range( real_type const xmin, real_type const xmax );

    //!
    //! Use externally allocated memory for `npts` points
    //!
    void reserve_external( integer const n, real_type *& p_x, real_type *& p_y, real_type *& p_Yp, real_type *& p_Ypp )
    {
      if ( !m_external_alloc ) m_base_quintic.free();
      m_npts           = 0;
      m_npts_reserved  = n;
      m_external_alloc = true;
      m_X              = p_x;
      m_Y              = p_y;
      m_Yp             = p_Yp;
      m_Ypp            = p_Ypp;
    }

    // --------------------------- VIRTUALS -----------------------------------

    //!
    //! \name Evaluation Aliases
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

    [[nodiscard]] real_type DDDD( real_type const x ) const override
    {
      std::pair<integer, real_type> res( 0, x );
      m_search.find( res );
      return this->id_DDDD( res.first, res.second );
    }

    [[nodiscard]] real_type DDDDD( real_type const x ) const override
    {
      std::pair<integer, real_type> res( 0, x );
      m_search.find( res );
      return this->id_DDDDD( res.first, res.second );
    }

    void D( real_type const x, real_type dd[2] ) const override;

    void DD( real_type const x, real_type dd[3] ) const override;

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

    //!
    //! \name Evaluation when segment is known
    ///@{
    [[nodiscard]] real_type id_eval( integer const ni, real_type const x ) const override;
    [[nodiscard]] real_type id_D( integer const ni, real_type const x ) const override;
    [[nodiscard]] real_type id_DD( integer const ni, real_type const x ) const override;
    [[nodiscard]] real_type id_DDD( integer const ni, real_type const x ) const override;
    [[nodiscard]] real_type id_DDDD( integer const ni, real_type const x ) const override;
    [[nodiscard]] real_type id_DDDDD( integer const ni, real_type const x ) const override;
    ///@}

    void reserve( integer npts ) override
    {
      if ( m_external_alloc && npts <= m_npts_reserved )
      {
        // nothing to do!, already allocated
      }
      else
      {
        m_base_quintic.reallocate( 4 * npts );
        m_npts_reserved  = npts;
        m_external_alloc = false;
        m_X              = m_base_quintic( npts );
        m_Y              = m_base_quintic( npts );
        m_Yp             = m_base_quintic( npts );
        m_Ypp            = m_base_quintic( npts );
      }
      m_npts = 0;
    }

    void clear() override;

    //!
    //! Get the piecewise polinomials of the spline
    //!
    integer  // order
    coeffs( real_type cfs[], real_type nodes[], bool transpose = false ) const override;

    [[nodiscard]] integer order() const override { return 6; }

    [[nodiscard]] bool is_monotone() const { return check_quintic_spline_monotonicity( m_X, m_Y, m_Yp, m_Ypp, m_npts ); }

#ifdef SPLINES_BACK_COMPATIBILITY
    void                    copySpline( QuinticSplineBase const & S ) { this->copy_spline( S ); }
    [[nodiscard]] real_type ypNode( integer i ) const noexcept { return this->yp_node( i ); }
    [[nodiscard]] real_type yppNode( integer i ) const noexcept { return this->ypp_node( i ); }
    void                    setRange( real_type xmin, real_type xmax ) { this->set_range( xmin, xmax ); }
#endif
  };

}  // namespace Splines

#endif

//
// EOF: SplineQuinticBase.hxx
//
