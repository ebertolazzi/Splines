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
// FILE: SplineBiCubicBase.cc
//

#define AUTODIFF_SUPPORT
#include "Splines/Splines.hh"

namespace Splines
{

  void BiCubicSplineBase::load( integer const i, integer const j, Mat4x4 & bili3 ) const
  {
    // Copia a blocchi
    // Invece di copiare scalarmente, copiamo 4 blocchi 2x2.
    // Eigen ottimizzerà queste operazioni usando istruzioni SIMD.

    // Quadrante in alto a sinistra: Z
    bili3.topLeftCorner<2, 2>() = mZ.block<2, 2>( i, j );

    // Quadrante in alto a destra: DY
    bili3.topRightCorner<2, 2>() = mDY.block<2, 2>( i, j );

    // Quadrante in basso a sinistra: DX
    bili3.bottomLeftCorner<2, 2>() = mDX.block<2, 2>( i, j );

    // Quadrante in basso a destra: DXY
    bili3.bottomRightCorner<2, 2>() = mDXY.block<2, 2>( i, j );
  }

  real_type BiCubicSplineBase::eval( real_type const x, real_type const y ) const
  {
    Mat4x4 bili3;
    Vec4   u, v;

    auto [i, j, dx, dy, DX, DY] = find_patch( x, y );

    Hermite3( dx, DX, u.data() );
    Hermite3( dy, DY, v.data() );

    load( i, j, bili3 );

    return u.dot( bili3 * v );
  }

  void BiCubicSplineBase::D( real_type const x, real_type const y, real_type d[3] ) const
  {
    Mat4x4 M;
    Vec4   u, u_D, v, v_D;

    auto [i, j, dx, dy, DX, DY] = find_patch( x, y );

    Hermite3( dx, DX, u.data() );
    Hermite3_D( dx, DX, u_D.data() );

    Hermite3( dy, DY, v.data() );
    Hermite3_D( dy, DY, v_D.data() );

    load( i, j, M );

    Vec4 Mv = M * v;

    d[0] = u.dot( Mv );
    d[1] = u_D.dot( Mv );
    d[2] = u.dot( M * v_D );
  }

  void BiCubicSplineBase::DD( real_type const x, real_type const y, real_type dd[6] ) const
  {
    Mat4x4 M;
    Vec4   u, u_D, u_DD, v, v_D, v_DD;

    auto [i, j, dx, dy, DX, DY] = find_patch( x, y );

    Hermite3( dx, DX, u.data() );
    Hermite3_D( dx, DX, u_D.data() );
    Hermite3_DD( dx, DX, u_DD.data() );

    Hermite3( dy, DY, v.data() );
    Hermite3_D( dy, DY, v_D.data() );
    Hermite3_DD( dy, DY, v_DD.data() );

    load( i, j, M );

    Vec4 Mv   = M * v;
    Vec4 Mv_D = M * v_D;

    dd[0] = u.dot( Mv );
    dd[1] = u_D.dot( Mv );
    dd[2] = u.dot( Mv_D );
    dd[3] = u_DD.dot( Mv );
    dd[4] = u_D.dot( Mv_D );
    dd[5] = u.dot( M * v_DD );
  }

  real_type BiCubicSplineBase::Dx( real_type const x, real_type const y ) const
  {
    Mat4x4 M;
    Vec4   u_D, v;

    auto [i, j, dx, dy, DX, DY] = find_patch( x, y );

    Hermite3_D( dx, DX, u_D.data() );
    Hermite3( dy, DY, v.data() );

    load( i, j, M );

    return u_D.dot( M * v );
  }

  real_type BiCubicSplineBase::Dy( real_type const x, real_type const y ) const
  {
    Mat4x4 M;
    Vec4   u, v_D;

    auto [i, j, dx, dy, DX, DY] = find_patch( x, y );

    Hermite3( dx, DX, u.data() );
    Hermite3_D( dy, DY, v_D.data() );

    load( i, j, M );

    return u.dot( M * v_D );
  }

  real_type BiCubicSplineBase::Dxx( real_type const x, real_type const y ) const
  {
    Mat4x4 M;
    Vec4   u_DD, v;

    auto [i, j, dx, dy, DX, DY] = find_patch( x, y );

    Hermite3_DD( dx, DX, u_DD.data() );
    Hermite3( dy, DY, v.data() );

    load( i, j, M );

    return u_DD.dot( M * v );
  }

  real_type BiCubicSplineBase::Dxy( real_type const x, real_type const y ) const
  {
    Mat4x4 M;
    Vec4   u_D, v_D;

    auto [i, j, dx, dy, DX, DY] = find_patch( x, y );

    Hermite3_D( dx, DX, u_D.data() );
    Hermite3_D( dy, DY, v_D.data() );

    load( i, j, M );

    return u_D.dot( M * v_D );
  }

  real_type BiCubicSplineBase::Dyy( real_type const x, real_type const y ) const
  {
    Mat4x4 M;
    Vec4   u, v_DD;

    auto [i, j, dx, dy, DX, DY] = find_patch( x, y );

    Hermite3( dx, DX, u.data() );
    Hermite3_DD( dy, DY, v_DD.data() );

    load( i, j, M );

    return u.dot( M * v_DD );
  }

  autodiff::dual1st BiCubicSplineBase::eval( autodiff::dual1st const & x, autodiff::dual1st const & y ) const
  {
    Mat4x4            M;
    autodiff::dual1st u[4], v[4], Mv[4];

    std::pair<integer, real_type> X( 0, x ), Y( 0, y );
    m_search_x.find( X );
    m_search_y.find( Y );

    integer const i = X.first;
    integer const j = Y.first;

    autodiff::dual1st XX, YY;
    XX.val  = X.second;
    XX.grad = x.grad;
    YY.val  = Y.second;
    YY.grad = y.grad;

    autodiff::dual1st const dx = XX - mX.coeff( i );
    autodiff::dual1st const dy = YY - mY.coeff( j );
    real_type const         DX = mX.coeff( i + 1 ) - mX.coeff( i );
    real_type const         DY = mY.coeff( j + 1 ) - mY.coeff( j );

    Hermite3<autodiff::dual1st>( dx, DX, u );
    Hermite3<autodiff::dual1st>( dy, DY, v );

    load( i, j, M );

    Mv[0] = M.coeff( 0, 0 ) * v[0] + M.coeff( 0, 1 ) * v[1] + M.coeff( 0, 2 ) * v[2] + M.coeff( 0, 3 ) * v[3];

    Mv[1] = M.coeff( 1, 0 ) * v[0] + M.coeff( 1, 1 ) * v[1] + M.coeff( 1, 2 ) * v[2] + M.coeff( 1, 3 ) * v[3];

    Mv[2] = M.coeff( 2, 0 ) * v[0] + M.coeff( 2, 1 ) * v[1] + M.coeff( 2, 2 ) * v[2] + M.coeff( 2, 3 ) * v[3];

    Mv[3] = M.coeff( 3, 0 ) * v[0] + M.coeff( 3, 1 ) * v[1] + M.coeff( 3, 2 ) * v[2] + M.coeff( 3, 3 ) * v[3];

    return u[0] * Mv[0] + u[1] * Mv[1] + u[2] * Mv[2] + u[3] * Mv[3];
  }

  autodiff::dual2nd BiCubicSplineBase::eval( autodiff::dual2nd const & x, autodiff::dual2nd const & y ) const
  {
    Mat4x4            M;
    autodiff::dual2nd u[4], v[4], Mv[4];

    std::pair<integer, real_type> X( 0, x ), Y( 0, y );
    m_search_x.find( X );
    m_search_y.find( Y );

    integer const i = X.first;
    integer const j = Y.first;

    autodiff::dual2nd XX, YY;

    XX.val.val   = X.second;
    XX.val.grad  = x.val.grad;
    XX.grad.val  = x.grad.val;
    XX.grad.grad = x.grad.grad;

    YY.val.val   = Y.second;
    YY.val.grad  = y.val.grad;
    YY.grad.val  = y.grad.val;
    YY.grad.grad = y.grad.grad;

    autodiff::dual2nd const dx = XX - mX.coeff( i );
    autodiff::dual2nd const dy = YY - mY.coeff( j );
    real_type const         DX = mX.coeff( i + 1 ) - mX.coeff( i );
    real_type const         DY = mY.coeff( j + 1 ) - mY.coeff( j );

    Hermite3<autodiff::dual2nd>( dx, DX, u );
    Hermite3<autodiff::dual2nd>( dy, DY, v );

    load( i, j, M );

    Mv[0] = M.coeff( 0, 0 ) * v[0] + M.coeff( 0, 1 ) * v[1] + M.coeff( 0, 2 ) * v[2] + M.coeff( 0, 3 ) * v[3];

    Mv[1] = M.coeff( 1, 0 ) * v[0] + M.coeff( 1, 1 ) * v[1] + M.coeff( 1, 2 ) * v[2] + M.coeff( 1, 3 ) * v[3];

    Mv[2] = M.coeff( 2, 0 ) * v[0] + M.coeff( 2, 1 ) * v[1] + M.coeff( 2, 2 ) * v[2] + M.coeff( 2, 3 ) * v[3];

    Mv[3] = M.coeff( 3, 0 ) * v[0] + M.coeff( 3, 1 ) * v[1] + M.coeff( 3, 2 ) * v[2] + M.coeff( 3, 3 ) * v[3];

    return u[0] * Mv[0] + u[1] * Mv[1] + u[2] * Mv[2] + u[3] * Mv[3];
  }

}  // namespace Splines

//
// EOF: SplineBiCubicBase.cc
//
