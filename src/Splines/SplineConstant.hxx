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

namespace Splines {

  using std::copy_n;

  //! Picewise constants spline class
  class ConstantSpline : public Spline {
    Malloc_real m_mem_constant;
    bool        m_external_alloc{false};

  public:

    #ifndef DOXYGEN_SHOULD_SKIP_THIS
    using Spline::build;
    #endif

    //!
    //! Build an empty spline of `ConstantSpline` type
    //!
    //! \param name the name of the spline
    //!
    explicit
    ConstantSpline( string_view name = "ConstantSpline" );

    //!
    //! Spline destructor.
    //!
    ~ConstantSpline() override {}

    //! Use externally allocated memory for `npts` points
    void
    reserve_external(
      integer       n,
      real_type * & p_x,
      real_type * & p_y
    );

    // --------------------------- VIRTUALS -----------------------------------
    //!
    //! \name Build
    //!
    ///@{

    //!
    //! Build the spline with the data stored
    //!
    void build() override { m_search.reset(); }

    //!
    //! Build the spline with the data passed as arguments
    //!
    //! \param x    \f$ x \f$ coordinates of the points
    //! \param incx access elements as `x[0]`, `x[incx]`, `x[2*incx]`,...
    //! \param y    \f$ y \f$ coordinates of the points
    //! \param incy access elements as `y[0]`, `y[incx]`, `y[2*incx]`,...
    //! \param n    the number of the points
    //!
    void
    build(
      real_type const x[], integer incx,
      real_type const y[], integer incy,
      integer n
    ) override;
    ///@}

    //!
    //! \name Evaluate
    //!
    ///@{
    real_type eval( real_type x ) const override;
    real_type D( real_type ) const override { return 0; }
    real_type DD( real_type ) const override { return 0; }
    real_type DDD( real_type ) const override { return 0; }
    ///@}

    //!
    //! \name Evaluation when segment is known
    //!
    ///@{
    real_type id_eval( integer ni, real_type x ) const override;
    real_type id_D( integer, real_type ) const override { return 0; }
    real_type id_DD( integer, real_type ) const override { return 0; }
    real_type id_DDD( integer, real_type ) const override { return 0; }
    ///@}

    void write_to_stream( ostream_type & ) const override;
    SplineType1D type() const override { return SplineType1D::CONSTANT; }

    // --------------------------- VIRTUALS -----------------------------------

    void reserve( integer npts ) override;

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
    copy_spline( ConstantSpline const & S ) {
      ConstantSpline::reserve(S.m_npts);
      m_npts = S.m_npts;
      copy_n( S.m_X,  m_npts,   m_X );
      copy_n( S.m_Y,  m_npts-1, m_Y );
      copy_flags( S );
    }

  };
}

// EOF: SplineConstant.hxx

