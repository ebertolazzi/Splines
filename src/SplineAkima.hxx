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
 |      _    _    _                   ____        _ _
 |     / \  | | _(_)_ __ ___   __ _  / ___| _ __ | (_)_ __   ___
 |    / _ \ | |/ / | '_ ` _ \ / _` | \___ \| '_ \| | | '_ \ / _ \
 |   / ___ \|   <| | | | | | | (_| |  ___) | |_) | | | | | |  __/
 |  /_/   \_\_|\_\_|_| |_| |_|\__,_| |____/| .__/|_|_|_| |_|\___|
 |                                         |_|
\*/

namespace Splines {

  #ifndef DOXYGEN_SHOULD_SKIP_THIS

  void
  Akima_build(
    real_type const * X,
    real_type const * Y,
    real_type       * Yp,
    integer           npts
  );

  #endif

  //! Akima spline class
  /*!
   |  Reference
   |  =========
   |  Hiroshi Akima, Journal of the ACM, Vol. 17, No. 4, October 1970, pages 589-602.
   */
  class AkimaSpline : public CubicSplineBase {
  public:

    #ifndef DOXYGEN_SHOULD_SKIP_THIS
    using CubicSplineBase::reserve;
    using CubicSplineBase::build;
    #endif

    //! spline constructor
    AkimaSpline( string const & name = "AkimaSpline" )
    : CubicSplineBase( name )
    {}

    //! spline destructor
    virtual
    ~AkimaSpline() override
    {}

    //! Return spline type (as number)
    virtual
    unsigned
    type() const override
    { return AKIMA_TYPE; }

    // --------------------------- VIRTUALS -----------------------------------

    //! Build an Akima spline from previously inserted points
    virtual
    void
    build() override;

    virtual
    void
    setup( GenericContainer const & gc ) override;

  };

}

// EOF: SplineAkima.hxx
