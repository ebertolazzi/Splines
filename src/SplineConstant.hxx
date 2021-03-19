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

    using Spline::build;

    ConstantSpline( string const & name = "ConstantSpline" )
    : Spline(name)
    , m_baseValue(name+"_memory")
    , m_external_alloc(false)
    {}

    ~ConstantSpline() override
    {}

    //! Use externally allocated memory for `npts` points
    void
    reserve_external(
      integer       n,
      real_type * & p_x,
      real_type * & p_y
    );

    // --------------------------- VIRTUALS -----------------------------------
    //! Build a spline.
    virtual
    void
    internal_build() override
    {} // nothing to do

    virtual
    void
    build(
      real_type const x[], integer incx,
      real_type const y[], integer incy,
      integer n
    ) override;

    //! Evaluate spline value at `x`
    virtual
    real_type
    operator () ( real_type x ) const override;

    //! First derivative
    virtual
    real_type
    D( real_type ) const override
    { return 0; }

    //! Second derivative
    virtual
    real_type
    DD( real_type ) const override
    { return 0; }

    //! Third derivative
    virtual
    real_type
    DDD( real_type ) const override
    { return 0; }

    //! Evaluate spline value at `x` knowing interval
    virtual
    real_type
    id_eval( integer ni, real_type x ) const override;

    //! First derivative
    virtual
    real_type
    id_D( integer, real_type ) const override
    { return 0; }

    //! Second derivative
    virtual
    real_type
    id_DD( integer, real_type ) const override
    { return 0; }

    //! Third derivative
    virtual
    real_type
    id_DDD( integer, real_type ) const override
    { return 0; }

    //! Print spline coefficients
    virtual
    void
    writeToStream( ostream_type & ) const override;

    //! Return spline type (as number)
    virtual
    unsigned
    type() const override
    { return CONSTANT_TYPE; }

    // --------------------------- VIRTUALS -----------------------------------

    //! Allocate memory for `npts` points
    virtual
    void
    reserve( integer npts ) override;

    //! Cancel the support points, empty the spline.
    virtual
    void
    clear() override;

    //! get the piecewise polinomials of the spline
    virtual
    integer // order
    coeffs(
      real_type cfs[],
      real_type nodes[],
      bool      transpose = false
    ) const override;

    virtual
    integer // order
    order() const override;

    virtual
    void
    setup( GenericContainer const & gc ) override;

  };
}

// EOF: SplineConstant.hxx

