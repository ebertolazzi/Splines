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

#include "Splines.hh"

//! Various kind of splines
namespace Splines {

  void
  Spline2D::build(
    SplineType2D    tp,
    real_type const x[], integer incx,
    real_type const y[], integer incy,
    real_type const z[], integer ldZ,
    integer nx, integer ny,
    bool fortran_storage,
    bool transposed
  ) {
    switch ( tp ) {
    case BILINEAR_TYPE:
      static_cast<BilinearSpline*>(m_pSpline2D)->build(
        x, incx, y, incy, z, ldZ, nx, ny,
        fortran_storage, transposed
      );
      break;
    case BICUBIC_TYPE:
      static_cast<BiCubicSpline*>(m_pSpline2D)->build(
        x, incx, y, incy, z, ldZ, nx, ny,
        fortran_storage, transposed
      );
      break;
    case BIQUINTIC_TYPE:
      static_cast<BiQuinticSpline*>(m_pSpline2D)->build(
        x, incx, y, incy, z, ldZ, nx, ny,
        fortran_storage, transposed
      );
      break;
    case AKIMA2D_TYPE:
      static_cast<Akima2Dspline*>(m_pSpline2D)->build(
        x, incx, y, incy, z, ldZ, nx, ny,
        fortran_storage, transposed
      );
      break;
    }
  }

  void
  Spline2D::build(
    SplineType2D              tp,
    vector<real_type> const & x,
    vector<real_type> const & y,
    vector<real_type> const & z,
    bool fortran_storage,
    bool transposed
  ) {
    switch ( tp ) {
    case BILINEAR_TYPE:
      static_cast<BilinearSpline*>(m_pSpline2D)->build(
        x, y, z, fortran_storage, transposed
      );
      break;
    case BICUBIC_TYPE:
      static_cast<BiCubicSpline*>(m_pSpline2D)->build(
        x, y, z, fortran_storage, transposed
      );
      break;
    case BIQUINTIC_TYPE:
      static_cast<BiQuinticSpline*>(m_pSpline2D)->build(
        x, y, z, fortran_storage, transposed
      );
      break;
    case AKIMA2D_TYPE:
      static_cast<Akima2Dspline*>(m_pSpline2D)->build(
        x, y, z, fortran_storage, transposed
      );
      break;
    }
  }

  void
  Spline2D::build(
    SplineType2D    tp,
    real_type const z[],
    integer         ldZ,
    integer         nx,
    integer         ny,
    bool fortran_storage,
    bool transposed
  ) {
    switch ( tp ) {
    case BILINEAR_TYPE:
      static_cast<BilinearSpline*>(m_pSpline2D)->build(
        z, ldZ, nx, ny, fortran_storage, transposed
      );
      break;
    case BICUBIC_TYPE:
      static_cast<BiCubicSpline*>(m_pSpline2D)->build(
        z, ldZ, nx, ny, fortran_storage, transposed
      );
      break;
    case BIQUINTIC_TYPE:
      static_cast<BiQuinticSpline*>(m_pSpline2D)->build(
        z, ldZ, nx, ny, fortran_storage, transposed
      );
      break;
    case AKIMA2D_TYPE:
      static_cast<Akima2Dspline*>(m_pSpline2D)->build(
        z, ldZ, nx, ny, fortran_storage, transposed
      );
      break;
    }
  }

  void
  Spline2D::build(
    SplineType2D              tp,
    vector<real_type> const & z,
    integer                   nx,
    integer                   ny,
    bool fortran_storage,
    bool transposed
  ) {
    switch ( tp ) {
    case BILINEAR_TYPE:
      static_cast<BilinearSpline*>(m_pSpline2D)->build(
        z, nx, ny, fortran_storage, transposed
      );
      break;
    case BICUBIC_TYPE:
      static_cast<BiCubicSpline*>(m_pSpline2D)->build(
        z, nx, ny, fortran_storage, transposed
      );
      break;
    case BIQUINTIC_TYPE:
      static_cast<BiQuinticSpline*>(m_pSpline2D)->build(
        z, nx, ny, fortran_storage, transposed
      );
      break;
    case AKIMA2D_TYPE:
      static_cast<Akima2Dspline*>(m_pSpline2D)->build(
        z, nx, ny, fortran_storage, transposed
      );
      break;
    }
  }

}
