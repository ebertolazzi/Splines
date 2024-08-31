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
 |    _   _                     _ _       ____        _ _
 |   | | | | ___ _ __ _ __ ___ (_) |_ ___/ ___| _ __ | (_)_ __   ___
 |   | |_| |/ _ \ '__| '_ ` _ \| | __/ _ \___ \| '_ \| | | '_ \ / _ \
 |   |  _  |  __/ |  | | | | | | | ||  __/___) | |_) | | | | | |  __/
 |   |_| |_|\___|_|  |_| |_| |_|_|\__\___|____/| .__/|_|_|_| |_|\___|
 |                                             |_|
\*/

namespace Splines {

  //! Hermite Spline Management Class
  class HermiteSpline : public CubicSplineBase {
  public:

    #ifndef DOXYGEN_SHOULD_SKIP_THIS
    using CubicSplineBase::build;
    using CubicSplineBase::reserve;
    #endif

    //!
    //! Build an empty spline of `HermiteSpline` type
    //!
    //! \param name the name of the spline
    //!
    HermiteSpline( string const & name = "HermiteSpline" )
    : CubicSplineBase( name )
    {}

    //!
    //! Spline destructor.
    //!
    ~HermiteSpline() override {}

    //! Return spline type (as number)
    SplineType1D type() const override { return SplineType1D::HERMITE; }

    // --------------------------- VIRTUALS -----------------------------------

    void build() override {} // nothing to do

    // block method!
    void
    build(
      real_type const [], integer,
      real_type const [], integer,
      integer
    ) override;

    void setup( GenericContainer const & gc ) override;

  };

}

// EOF: SplineHermite.hxx
