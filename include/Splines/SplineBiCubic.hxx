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

#ifndef SPLINE_BICUBIC_HXX
#define SPLINE_BICUBIC_HXX

namespace Splines
{

  /*\
   |   ____  _  ____      _     _      ____        _ _
   |  | __ )(_)/ ___|   _| |__ (_) ___/ ___| _ __ | (_)_ __   ___
   |  |  _ \| | |  | | | | '_ \| |/ __\___ \| '_ \| | | '_ \ / _ \
   |  | |_) | | |__| |_| | |_) | | (__ ___) | |_) | | | | | |  __/
   |  |____/|_|\____\__,_|_.__/|_|\___|____/| .__/|_|_|_| |_|\___|
   |                                        |_|
  \*/
  //!
  //! Cubic spline base class
  //!
  class BiCubicSpline final : public BiCubicSplineBase
  {
    using BiCubicSplineBase::mDX;
    using BiCubicSplineBase::mDXY;
    using BiCubicSplineBase::mDY;

    void make_spline() override;

  public:
    using BiCubicSplineBase::eval;

    //!
    //! Build an empty spline of `BiCubicSpline` type
    //!
    //! \param type spline type
    //! \param name the name of the spline
    //!
    explicit BiCubicSpline( Spline_sub_type sub_type = Spline_sub_type::PCHIP, string_view name = "BiCubicSpline" )
      : BiCubicSplineBase( sub_type, name )
    {
    }

    //!
    //! Build an empty spline of `BiCubicSpline` type
    //!
    //! \param name the name of the spline
    //!
    explicit BiCubicSpline( string_view name ) : BiCubicSplineBase( Spline_sub_type::PCHIP, name ) {}

    //!
    //! Spline destructor.
    //!
    ~BiCubicSpline() override {}

    void write_to_stream( ostream_type & s ) const override;

    [[nodiscard]] char const * type_name() const override { return "BiCubic"; }
  };

}  // namespace Splines

#endif

//
// EOF: SplineBiCubic.hxx
//
