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
 |   _     _                       ____        _ _
 |  | |   (_)_ __   ___  __ _ _ __/ ___| _ __ | (_)_ __   ___
 |  | |   | | '_ \ / _ \/ _` | '__\___ \| '_ \| | | '_ \ / _ \
 |  | |___| | | | |  __/ (_| | |   ___) | |_) | | | | | |  __/
 |  |_____|_|_| |_|\___|\__,_|_|  |____/| .__/|_|_|_| |_|\___|
 |                                      |_|
\*/

namespace Splines {

  //! Linear spline class
  class LinearSpline : public Spline {
    Utils::Malloc<real_type> m_baseValue;
    bool                     m_external_alloc;

  public:

    #ifndef DOXYGEN_SHOULD_SKIP_THIS
    using Spline::build;
    #endif

    LinearSpline( string const & name = "LinearSpline" )
    : Spline(name)
    , m_baseValue( name+"_memory")
    , m_external_alloc(false)
    {
      m_curve_extended_constant = true; // by default linear spline extend constant
    }

    virtual
    ~LinearSpline() override
    {}

    //! Use externally allocated memory for `npts` points
    void
    reserve_external(
      integer      n,
      real_type *& p_x,
      real_type *& p_y
    );

    // --------------------------- VIRTUALS -----------------------------------

    //! Evalute spline value at `x`
    virtual
    real_type
    operator () ( real_type x ) const override;

    //! First derivative
    virtual
    real_type
    D( real_type x ) const override;

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

    //! Evaluate spline value knowing interval
    virtual
    real_type
    id_eval( integer ni, real_type x ) const override;

    //! First derivative
    virtual
    real_type
    id_D( integer, real_type ) const override;

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
    writeToStream( ostream_type & s ) const override;

    //! Return spline type (as number)
    virtual
    unsigned
    type() const override
    { return LINEAR_TYPE; }

    // --------------------------- VIRTUALS -----------------------------------

    //! Allocate memory for `npts` points
    virtual
    void
    reserve( integer npts ) override;

    //! added for compatibility with cubic splines
    virtual
    void
    build() override
    {}

    //! Cancel the support points, empty the spline.
    virtual
    void
    clear() override;

    //! get the piecewise polinomials of the spline
    virtual
    integer // order
    coeffs(
      real_type * const cfs,
      real_type * const nodes,
      bool              transpose = false
    ) const override;

    virtual
    integer // order
    order() const override;

    virtual
    void
    setup( GenericContainer const & gc ) override;

  };

}

// EOF: SplineLinbear.hxx
