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

#pragma once

#ifndef SPLINE_VANLEER_HXX
#define SPLINE_VANLEER_HXX

/*\
 |  __     __          _                   ____        _ _
 |  \ \   / /_ _ _ __ | |    ___  ___ _ __/ ___| _ __ | (_)_ __   ___
 |   \ \ / / _` | '_ \| |   / _ \/ _ \ '__\___ \| '_ \| | | '_ \ / _ \
 |    \ V / (_| | | | | |__|  __/  __/ |   ___) | |_) | | | | | |  __/
 |     \_/ \__,_|_| |_|_____\___|\___|_|  |____/| .__/|_|_|_| |_|\___|
 |                                              |_|
\*/

namespace Splines
{

  //!
  //! Van Leer spline class
  //!
  class VanLeerSpline final : public CubicSplineBase
  {
  public:
    using CubicSplineBase::build;
    using CubicSplineBase::reserve;

    //!
    //! Build an empty spline of `VanLeerSpline` type
    //!
    //! \param name the name of the spline
    //!
    explicit VanLeerSpline( string_view name = "VanLeerSpline" ) : CubicSplineBase( name ) {}

    //!
    //! spline destructor
    //!
    ~VanLeerSpline() override {}

    //!
    //! Return spline type (as number)
    //!
    [[nodiscard]] SplineType1D type() const override { return SplineType1D::VANLEER; }

    // --------------------------- VIRTUALS -----------------------------------

    void build() override;
    void setup( GenericContainer const & gc ) override;
  };

}  // namespace Splines

#endif
