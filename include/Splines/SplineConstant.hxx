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
 |    ____                _              _       ____        _ _
 |   / ___|___  _ __  ___| |_ __ _ _ __ | |_ ___/ ___| _ __ | (_)_ __   ___
 |  | |   / _ \| '_ \/ __| __/ _` | '_ \| __/ __\___ \| '_ \| | | '_ \ / _ \
 |  | |__| (_) | | | \__ \ || (_| | | | | |_\__ \___) | |_) | | | | | |  __/
 |   \____\___/|_| |_|___/\__\__,_|_| |_|\__|___/____/| .__/|_|_|_| |_|\___|
 |                                                    |_|
\*/

namespace Splines
{

  using std::copy_n;

  //! Picewise constants spline class
  class ConstantSpline final : public Spline
  {
    Malloc_real m_mem_constant;
    bool        m_external_alloc = false;

  public:
    using Spline::build;

    //!
    //! Build an empty spline of `ConstantSpline` type
    //!
    //! \param name the name of the spline
    //!
    explicit ConstantSpline( string_view name = "ConstantSpline" )
      : Spline( name ), m_mem_constant( fmt::format( "ConstantSpline[{}]", name ) )
    {
    }

    //!
    //! Spline destructor.
    //!
    ~ConstantSpline() override {}

    //! Use externally allocated memory for `npts` points
    void reserve_external( integer const n, real_type *& p_x, real_type *& p_y )
    {
      if ( !m_external_alloc ) m_mem_constant.free();
      m_npts           = 0;
      m_npts_reserved  = n;
      m_external_alloc = true;
      m_X              = p_x;
      m_Y              = p_y;
    }

    // --------------------------- VIRTUALS -----------------------------------
    //!
    //! \name Build
    //!
    ///@{

    //!
    //! Build the spline with the data stored
    //!
    void build() override { m_search.must_reset(); }

    ///@}

    //!
    //! \name Evaluate
    //!
    ///@{
    [[nodiscard]] real_type eval( real_type const x ) const override
    {
      std::pair<integer, real_type> res( 0, x );
      m_search.find( res );
      return m_Y[res.first];
    }

    [[nodiscard]] real_type D( real_type const ) const override { return 0; }
    [[nodiscard]] real_type DD( real_type const ) const override { return 0; }
    [[nodiscard]] real_type DDD( real_type const ) const override { return 0; }

    void D( real_type const x, real_type dd[2] ) const override
    {
      dd[0] = eval( x );
      dd[1] = 0;
    }

    void DD( real_type const x, real_type dd[3] ) const override
    {
      dd[0] = eval( x );
      dd[1] = 0;
      dd[2] = 0;
    }

    ///@}

    //!
    //! \name Evaluation when segment is known
    //!
    ///@{
    [[nodiscard]] real_type id_eval( integer const ni, [[maybe_unused]] real_type const x ) const override { return m_Y[ni]; }
    [[nodiscard]] real_type id_D( [[maybe_unused]] integer const ni, [[maybe_unused]] real_type const x ) const override { return 0; }
    [[nodiscard]] real_type id_DD( [[maybe_unused]] integer const ni, [[maybe_unused]] real_type const x ) const override
    { return 0; }
    [[nodiscard]] real_type id_DDD( [[maybe_unused]] integer const ni, [[maybe_unused]] real_type const x ) const override
    { return 0; }
    ///@}

#ifdef AUTODIFF_SUPPORT
    [[nodiscard]] autodiff::dual1st eval( autodiff::dual1st const & x ) const override
    {
      autodiff::dual1st res;
      res.val  = eval( x.val );
      res.grad = 0;
      return res;
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    [[nodiscard]] autodiff::dual2nd eval( autodiff::dual2nd const & x ) const override
    {
      autodiff::dual2nd res;
      res.val.val   = eval( x.val.val );
      res.val.grad  = 0;
      res.grad.val  = 0;
      res.grad.grad = 0;
      return res;
    }
#endif

    void write_to_stream( ostream_type & s ) const override;

    [[nodiscard]] SplineType1D type() const override { return SplineType1D::CONSTANT; }

    // --------------------------- VIRTUALS -----------------------------------

    void reserve( integer const npts ) override;

    void clear() override;

    integer coeffs( real_type cfs[], real_type nodes[], [[maybe_unused]] bool transpose = false ) const override
    {
      // 1. Early exit: Se non ci sono punti, non fare nulla.
      if ( m_npts <= 0 ) return 1;

      // 2. Ottimizzazione Memoria: Usa memcpy per copia diretta di byte.
      // È la via più veloce per copiare array di numeri (bypassando gli iteratori).
      std::memcpy( nodes, m_X, m_npts * sizeof( real_type ) );

      // 3. Gestione Coeffs: Copiamo m_Y solo se ci sono segmenti (punti > 1).
      // Nota: Se m_npts è 1, nseg è 0, quindi non copiamo nulla in cfs.
      if ( m_npts > 1 ) std::memcpy( cfs, m_Y, ( m_npts - 1 ) * sizeof( real_type ) );

      return 1;
    }

    [[nodiscard]] integer order() const override { return 1; }

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

    void copy_spline( ConstantSpline const & S );
  };
}  // namespace Splines

//
// EOF: SplineConstant.cc
//
