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
 |    ____                _              _       ____        _ _
 |   / ___|___  _ __  ___| |_ __ _ _ __ | |_ ___/ ___| _ __ | (_)_ __   ___
 |  | |   / _ \| '_ \/ __| __/ _` | '_ \| __/ __\___ \| '_ \| | | '_ \ / _ \
 |  | |__| (_) | | | \__ \ || (_| | | | | |_\__ \___) | |_) | | | | | |  __/
 |   \____\___/|_| |_|___/\__\__,_|_| |_|\__|___/____/| .__/|_|_|_| |_|\___|
 |                                                    |_|
\*/

namespace Splines {

  //! Picewise constants spline class
  class ConstantSpline : public Spline {
    Utils::Malloc<real_type> m_baseValue;
    bool                     m_external_alloc;

  public:

    #ifndef DOXYGEN_SHOULD_SKIP_THIS
    using Spline::build;
    #endif

    ConstantSpline( string const & name = "ConstantSpline" )
    : Spline(name)
    , m_baseValue(name+"_memory")
    , m_external_alloc(false)
    {}

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
    void build() override {} // nothing to do

    void
    build(
      real_type const * x, integer incx,
      real_type const * y, integer incy,
      integer n
    ) override;
    ///@}

    //!
    //! \name Evaluate
    //!
    ///@{
    real_type operator () ( real_type x ) const override;
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
    unsigned type() const override { return CONSTANT_TYPE; }

    // --------------------------- VIRTUALS -----------------------------------

    void reserve( integer npts ) override;

    void clear() override;

    integer // order
    coeffs(
      real_type * const cfs,
      real_type * const nodes,
      bool              transpose = false
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

  };
}

// EOF: SplineConstant.hxx

