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
    real_type const X[],
    real_type const Y[],
    real_type       Yp[],
    integer         npts
  );

  #endif

  //!
  //! Smooth Curve Fitting Based on Local Procedures
  //!
  //! *Reference*
  //!
  //! - *Hiroshi Akima*, Journal of the ACM, Vol.17, No. 4, 589-602, 1970.
  //!
  class AkimaSpline : public CubicSplineBase {
  public:

    #ifndef DOXYGEN_SHOULD_SKIP_THIS
    using CubicSplineBase::reserve;
    using CubicSplineBase::build;
    #endif

    //!
    //! Build an empty spline of `AkimaSpline` type
    //!
    //! \param name the name of the spline
    //!
    AkimaSpline( string_view name = "AkimaSpline" )
    : CubicSplineBase( name )
    {}

    //!
    //! Spline destructor.
    //!
    ~AkimaSpline() override {}

    //!
    //! Return spline type (as number).
    //!
    SplineType1D type() const override { return SplineType1D::AKIMA; }

    // --------------------------- VIRTUALS -----------------------------------

    void build() override;
    void setup( GenericContainer const & gc ) override;

  };

}

// EOF: SplineAkima.hxx
