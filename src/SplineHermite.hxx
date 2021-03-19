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

    using CubicSplineBase::build;
    using CubicSplineBase::reserve;

    //! spline constructor
    HermiteSpline( string const & name = "HermiteSpline" )
    : CubicSplineBase( name )
    {}

    //! spline destructor
    virtual
    ~HermiteSpline() override
    {}

    //! Return spline type (as number)
    virtual
    unsigned
    type() const override
    { return HERMITE_TYPE; }

    // --------------------------- VIRTUALS -----------------------------------

    virtual
    void
    internal_build() override
    {} // nothing to do

    // block method!
    virtual
    void
    build(
      real_type const [], integer,
      real_type const [], integer,
      integer
    ) override;

    virtual
    void
    setup( GenericContainer const & gc ) override;

  };

}

// EOF: SplineHermite.hxx