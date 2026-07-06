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

//
// FILE: SplineBiCubicBase.hxx
//

#pragma once

#ifndef SPLINE_BICUBIC_BASE_HXX
#define SPLINE_BICUBIC_BASE_HXX

namespace Splines
{

  /*\
   |   ____  _  ____      _     _      ____        _ _            ____
   |  | __ )(_)/ ___|   _| |__ (_) ___/ ___| _ __ | (_)_ __   ___| __ )  __ _ ___  ___
   |  |  _ \| | |  | | | | '_ \| |/ __\___ \| '_ \| | | '_ \ / _ \  _ \ / _` / __|/ _ \
   |  | |_) | | |__| |_| | |_) | | (__ ___) | |_) | | | | | |  __/ |_) | (_| \__ \  __/
   |  |____/|_|\____\__,_|_.__/|_|\___|____/| .__/|_|_|_| |_|\___|____/ \__,_|___/\___|
   |                                        |_|
  \*/

  //!
  //! Bi-cubic spline base class
  //!
  class BiCubicSplineBase : public SplineSurf
  {
  protected:
    MatC mDX;
    MatC mDY;
    MatC mDXY;

    Spline_sub_type m_sub_type;

    using SplineSurf::m_nx;
    using SplineSurf::m_ny;

    using SplineSurf::mX;
    using SplineSurf::mY;
    using SplineSurf::mZ;

    void load( integer const i, integer const j, Mat4x4 & bili3 ) const;

    real_type & Dx_node_ref( integer const i, integer const j ) { return mDX.coeffRef( i, j ); }
    real_type & Dy_node_ref( integer const i, integer const j ) { return mDY.coeffRef( i, j ); }
    real_type & Dxy_node_ref( integer const i, integer const j ) { return mDXY.coeffRef( i, j ); }

  public:
    using SplineSurf::eval;

    //! spline constructor
    explicit BiCubicSplineBase(
      Spline_sub_type sub_type = Spline_sub_type::PCHIP,
      string_view     name     = "BiCubicSplineBase" )
      : SplineSurf( name ), m_sub_type( sub_type )
    {
    }

    ~BiCubicSplineBase() override {}

    //!
    //! \name Estimated derivatives at interpolation nodes
    //!
    ///@{

    //!
    //! Estimated `x` derivatives at node `(i,j)`
    //!
    [[nodiscard]] real_type Dx_node( integer const i, integer const j ) const { return mDX.coeff( i, j ); }

    //!
    //! Estimated `y` derivatives at node `(i,j)`
    //!
    [[nodiscard]] real_type Dy_node( integer const i, integer const j ) const { return mDY.coeff( i, j ); }

    //!
    //! Estimated mixed `xy` derivatives at node `(i,j)`
    //!
    [[nodiscard]] real_type Dxy_node( integer const i, integer const j ) const { return mDXY.coeff( i, j ); }

    ///@}

    //!
    //! \name Evaluate
    //!
    ///@{
    //!
    //! Evaluate spline at point \f$ (x,y) \f$
    //!
    [[nodiscard]] real_type eval( real_type const x, real_type const y ) const override;

    //!
    //! Evaluate spline with derivative at point \f$ (x,y) \f$
    //!
    //! - `d[0]` the value of the spline
    //! - `d[1]` the value of the spline `x` derivative
    //! - `d[2]` the value of the spline `y` derivative
    //!
    void D( real_type const x, real_type const y, real_type d[3] ) const override;

    //!
    //! Evaluate spline with derivative at point \f$ (x,y) \f$
    //!
    //! - `d[0]` the value of the spline
    //! - `d[1]` the value of the spline `x` derivative
    //! - `d[2]` the value of the spline `y` derivative
    //! - `d[3]` the value of the spline `x` second derivative
    //! - `d[4]` the value of the spline `xy` mixed derivative
    //! - `d[5]` the value of the spline `y` second derivative
    //!
    void DD( real_type const x, real_type const y, real_type dd[6] ) const override;

    //!
    //! Evaluate spline `x`  derivative at point \f$ (x,y) \f$
    //!
    [[nodiscard]] real_type Dx( real_type const x, real_type const y ) const override;

    //!
    //! Evaluate spline `y`  derivative at point \f$ (x,y) \f$
    //!
    [[nodiscard]] real_type Dy( real_type const x, real_type const y ) const override;

    //!
    //! Evaluate spline `x` second derivative at point \f$ (x,y) \f$
    //!
    [[nodiscard]] real_type Dxx( real_type const x, real_type const y ) const override;

    //!
    //! Evaluate spline `xy` mixed derivative at point \f$ (x,y) \f$
    //!
    [[nodiscard]] real_type Dxy( real_type const x, real_type const y ) const override;

    //!
    //! Evaluate spline `y` second derivative at point \f$ (x,y) \f$
    //!
    [[nodiscard]] real_type Dyy( real_type const x, real_type const y ) const override;

#ifdef AUTODIFF_SUPPORT
    //!
    //! \name Autodiff
    //!
    ///@{
    [[nodiscard]] autodiff::dual1st eval( autodiff::dual1st const & x, autodiff::dual1st const & y ) const override;
    [[nodiscard]] autodiff::dual2nd eval( autodiff::dual2nd const & x, autodiff::dual2nd const & y ) const override;
    ///@}
#endif

    ///@}
  };

}  // namespace Splines

#endif

//
// EOF: SplineBiCubicBase.hxx
//
