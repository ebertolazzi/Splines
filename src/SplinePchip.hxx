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
 |   ____      _     _      ____        _ _
 |  |  _ \ ___| |__ (_)_ __/ ___| _ __ | (_)_ __   ___
 |  | |_) / __| '_ \| | '_ \___ \| '_ \| | | '_ \ / _ \
 |  |  __/ (__| | | | | |_) |__) | |_) | | | | | |  __/
 |  |_|   \___|_| |_|_| .__/____/| .__/|_|_|_| |_|\___|
 |                    |_|        |_|
\*/

namespace Splines {
  void
  Pchip_build(
    real_type const * X,
    real_type const * Y,
    real_type       * Yp,
    integer           npts
  );

  //! Pchip (Piecewise Cubic Hermite Interpolating Polynomial) spline class
  class PchipSpline : public CubicSplineBase {
  public:

    #ifndef DOXYGEN_SHOULD_SKIP_THIS
    using CubicSplineBase::build;
    using CubicSplineBase::reserve;
    #endif

    //! spline constructor
    PchipSpline( string const & name = "PchipSpline" )
    : CubicSplineBase( name )
    {}

    //! spline destructor
    virtual
    ~PchipSpline() override
    {}

    //! Return spline type (as number)
    virtual
    unsigned
    type() const override
    { return PCHIP_TYPE; }

    // --------------------------- VIRTUALS -----------------------------------

    //! Build a Monotone spline from previously inserted points
    virtual
    void
    build() override;

    virtual
    void
    setup( GenericContainer const & gc ) override;

  };

}

// EOF: SplinePchip.hxx
