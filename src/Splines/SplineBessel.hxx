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
 |   ____                     _ ____        _ _
 |  | __ )  ___  ___ ___  ___| / ___| _ __ | (_)_ __   ___
 |  |  _ \ / _ \/ __/ __|/ _ \ \___ \| '_ \| | | '_ \ / _ \
 |  | |_) |  __/\__ \__ \  __/ |___) | |_) | | | | | |  __/
 |  |____/ \___||___/___/\___|_|____/| .__/|_|_|_| |_|\___|
 |                                   |_|
\*/

namespace Splines {

  void
  Bessel_build(
    real_type const X[],
    real_type const Y[],
    real_type       Yp[],
    integer   const npts
  );

  //!
  //! Bessel spline class
  //!
  class BesselSpline : public CubicSplineBase {
  public:

    using CubicSplineBase::build;
    using CubicSplineBase::reserve;

    //!
    //! Build an empty spline of `BesselSpline` type
    //!
    //! \param name the name of the spline
    //!
    explicit
    BesselSpline( string_view name = "BesselSpline" )
    : CubicSplineBase( name )
    {}

    //!
    //! spline destructor
    //!
    ~BesselSpline() override {}

    //!
    //! Return spline type (as number)
    //!
    SplineType1D type() const override { return SplineType1D::BESSEL; }

    // --------------------------- VIRTUALS -----------------------------------

    void build() override;
    void setup( GenericContainer const & gc ) override;
  };

}
