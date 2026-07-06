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

#ifndef SPLINE_HERMITE_HXX
#define SPLINE_HERMITE_HXX

namespace Splines
{

  /*\
   |    _   _                     _ _       ____        _ _
   |   | | | | ___ _ __ _ __ ___ (_) |_ ___/ ___| _ __ | (_)_ __   ___
   |   | |_| |/ _ \ '__| '_ ` _ \| | __/ _ \___ \| '_ \| | | '_ \ / _ \
   |   |  _  |  __/ |  | | | | | | | ||  __/___) | |_) | | | | | |  __/
   |   |_| |_|\___|_|  |_| |_| |_|_|\__\___|____/| .__/|_|_|_| |_|\___|
   |                                             |_|
  \*/

  //! Hermite Spline Management Class
  class HermiteSpline final : public CubicSplineBase
  {
  public:
    using CubicSplineBase::build;
    using CubicSplineBase::reserve;

    //!
    //! Build an empty spline of `HermiteSpline` type
    //!
    //! \param name the name of the spline
    //!
    HermiteSpline( string_view name = "HermiteSpline" ) : CubicSplineBase( name ) {}

    //!
    //! Spline destructor.
    //!
    ~HermiteSpline() override {}

    //! Return spline type (as number)
    [[nodiscard]] SplineType1D type() const override { return SplineType1D::HERMITE; }

    // --------------------------- VIRTUALS -----------------------------------

    void build() override { m_search.must_reset(); }

    //!
    //!
    //! Setup a spline using a `GenericContainer`
    //!
    //! - gc("xdata")  vector with the `x` coordinate of the data
    //! - gc("ydata")  vector with the `y` coordinate of the data
    //! - gc("ypdata") vector with the `y` derivative of the data
    //!
    void setup( GenericContainer const & gc ) override;
  };

}  // namespace Splines

#endif  // SPLINE_HERMITE_HXX
