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

namespace Splines {

  using std::copy_n;

  //! Linear spline class
  class LinearSpline : public Spline {
    Malloc_real m_mem_linear;
    bool        m_external_alloc{false};

  public:

    using Spline::build;

    //!
    //! Build an empty spline of `LinearSpline` type
    //!
    //! \param name the name of the spline
    //!
    explicit
    LinearSpline( string_view name = "LinearSpline" );

    //!
    //! Spline destructor.
    //!
    ~LinearSpline() override {}

    //! Use externally allocated memory for `npts` points
    void
    reserve_external(
      integer const n,
      real_type * & p_x,
      real_type * & p_y
    );

    // --------------------------- VIRTUALS -----------------------------------

    real_type eval ( real_type const x ) const override;
    real_type D    ( real_type const x ) const override;
    real_type DD   ( real_type const   ) const override { return 0; }
    real_type DDD  ( real_type const   ) const override { return 0; }

    void D  ( real_type const x, real_type dd[2] ) const override;
    void DD ( real_type const x, real_type dd[3] ) const override;

    real_type id_eval ( integer const ni, real_type const x ) const override;
    real_type id_D    ( integer const   , real_type const   ) const override;
    real_type id_DD   ( integer const   , real_type const   ) const override { return 0; }
    real_type id_DDD  ( integer const   , real_type const   ) const override { return 0; }

    void write_to_stream( ostream_type & s ) const override;
    SplineType1D type() const override { return SplineType1D::LINEAR; }

    // --------------------------- VIRTUALS -----------------------------------

    #ifdef AUTODIFF_SUPPORT
    autodiff::dual1st eval( autodiff::dual1st const & x ) const override;
    autodiff::dual2nd eval( autodiff::dual2nd const & x ) const override;

    template <typename T>
    autodiff::HigherOrderDual<autodiff::detail::DualOrder<T>::value,real_type>
    eval( T const & x ) const { return eval( autodiff::detail::to_dual(x) ); }

    template <typename T>
    autodiff::HigherOrderDual<autodiff::detail::DualOrder<T>::value,real_type>
    operator () ( T const & x ) const { return eval( autodiff::detail::to_dual(x) ); }
    #endif

    void reserve( integer const npts ) override;
    void build() override { m_search.must_reset(); }
    void clear() override;

    integer // order
    coeffs(
      real_type cfs[],
      real_type nodes[],
      bool      transpose = false
    ) const override;

    integer order() const override;
    void setup( GenericContainer const & gc ) override;

    void
    y_min_max(
      integer   & i_min_pos,
      real_type & x_min_pos,
      real_type & y_min,
      integer   & i_max_pos,
      real_type & x_max_pos,
      real_type & y_max
    ) const override;

    void
    y_min_max(
      vector<integer>   & i_min_pos,
      vector<real_type> & x_min_pos,
      vector<real_type> & y_min,
      vector<integer>   & i_max_pos,
      vector<real_type> & x_max_pos,
      vector<real_type> & y_max
    ) const override;

    void
    copy_spline( LinearSpline const & S ) {
      LinearSpline::reserve(S.m_npts);
      m_npts = S.m_npts;
      copy_n( S.m_X, m_npts, m_X );
      copy_n( S.m_Y, m_npts, m_Y );
      copy_flags( S );
    }

  };

}

// EOF: SplineLinear.hxx
