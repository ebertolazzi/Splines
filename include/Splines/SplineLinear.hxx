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
 |   _     _                       ____        _ _
 |  | |   (_)_ __   ___  __ _ _ __/ ___| _ __ | (_)_ __   ___
 |  | |   | | '_ \ / _ \/ _` | '__\___ \| '_ \| | | '_ \ / _ \
 |  | |___| | | | |  __/ (_| | |   ___) | |_) | | | | | |  __/
 |  |_____|_|_| |_|\___|\__,_|_|  |____/| .__/|_|_|_| |_|\___|
 |                                      |_|
\*/

namespace Splines
{

  //! Linear spline class
  class LinearSpline final : public Spline
  {
    Malloc_real m_mem_linear;
    bool        m_external_alloc = false;

  public:
    using Spline::build;

    //!
    //! Build an empty spline of `LinearSpline` type
    //!
    //! \param name the name of the spline
    //!
    explicit LinearSpline( string_view name = "LinearSpline" )
      : Spline( name ), m_mem_linear( fmt::format( "LinearSpline[{}]", name ) )
    {
      m_curve_extended_constant = true;  // by default linear spline extend constant
    }

    //!
    //! Spline destructor.
    //!
    ~LinearSpline() override {}

    //! Use externally allocated memory for `npts` points
    void reserve_external( integer const n, real_type *& p_x, real_type *& p_y )
    {
      if ( !m_external_alloc ) m_mem_linear.free();
      m_npts           = 0;
      m_npts_reserved  = n;
      m_external_alloc = true;
      m_X              = p_x;
      m_Y              = p_y;
    }

    // --------------------------- VIRTUALS -----------------------------------

    [[nodiscard]] real_type eval( real_type const x ) const override
    {
      std::pair<integer, real_type> res( 0, x );
      m_search.find( res );
      return this->id_eval( res.first, res.second );
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    [[nodiscard]] real_type D( real_type const x ) const override;
    [[nodiscard]] real_type DD( real_type const ) const override { return 0; }
    [[nodiscard]] real_type DDD( real_type const ) const override { return 0; }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    void D( real_type const x, real_type dd[2] ) const override;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    void DD( real_type const x, real_type dd[3] ) const override;

    [[nodiscard]] real_type id_eval( integer const ni, real_type const x ) const override;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    [[nodiscard]] real_type id_D( integer const i, real_type const x ) const override;
    [[nodiscard]] real_type id_DD( integer const, real_type const ) const override { return 0; }
    [[nodiscard]] real_type id_DDD( integer const, real_type const ) const override { return 0; }

    void write_to_stream( ostream_type & s ) const override;

    [[nodiscard]] SplineType1D type() const override { return SplineType1D::LINEAR; }

    // --------------------------- VIRTUALS -----------------------------------

#ifdef AUTODIFF_SUPPORT
    //!
    //! \name Autodiff
    //!
    ///@{
    [[nodiscard]] autodiff::dual1st eval( autodiff::dual1st const & x ) const override;
    [[nodiscard]] autodiff::dual2nd eval( autodiff::dual2nd const & x ) const override;
    ///@}
#endif

    void reserve( integer const npts ) override;

    void build() override { m_search.must_reset(); }

    void clear() override;

    integer  // order
    coeffs( real_type cfs[], real_type nodes[], bool const transpose ) const override;

    [[nodiscard]] integer order() const override { return 2; }

    void setup( GenericContainer const & gc ) override;

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

    void copy_spline( LinearSpline const & S );
  };

}  // namespace Splines

// EOF: SplineLinear.hxx
