/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 1998                                                      |
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

/**
 * 
 */

namespace Splines {

  using namespace std ; // load standard namspace

  void
  SplineSet::build ( indexType  const nspl,
                     indexType  const npts,
                     char       const *headers[],
                     SplineType const stype[],
                     valueType  const X[],
                     valueType  const *Y[],
                     bool       const rp_policy[] ) {
    SPLINE_ASSERT( nspl > 0, "SplineSet::build expected positive nspl = " << nspl ) ;
    SPLINE_ASSERT( npts > 1, "SplineSet::build expected npts = " << npts << " greather than 1" ) ;
    _nspl = nspl ;
    _npts = npts ;
    // alloco spazio per slines
    splines.resize(sizeType(nspl)) ;
    indexType mem = npts ;
    for ( indexType i = 0 ; i < nspl ; ++i ) {
      switch (stype[i]) {
        case CONSTANT: case LINEAR:
          mem += npts ;
        break;

        case AKIMA: case BESSEL: case PCHIP: case CUBIC:
          mem += 2*npts ; // Y, Yp
        break;

        case QUINTIC:
          mem += 3*npts ; // Y, Yp, Ypp
        break;
        
        default:
          SPLINE_ASSERT( false, "SplineSet::build unknwn type = " << stype[i] ) ;
        break;
      }
    }

    baseValue.allocate( sizeType(mem + 2*nspl) ) ;
    basePointer.allocate( sizeType(3*nspl) ) ;

    _Y    = basePointer(_nspl) ;
    _Yp   = basePointer(_nspl) ;
    _Ypp  = basePointer(_nspl) ;
    _X    = baseValue(_npts) ;
    _Ymin = baseValue(_nspl) ;
    _Ymax = baseValue(_nspl) ;
    std::copy( X, X+npts, _X ) ;
    for ( indexType i = 0 ; i < nspl ; ++i ) {
      valueType *& pY   = _Y[i] ;
      valueType *& pYp  = _Yp[i] ;
      valueType *& pYpp = _Ypp[i] ;
      pY = baseValue( _npts ) ;
      std::copy( Y[i], Y[i]+npts, pY ) ;
      _Ymin[i] = *std::min_element( pY,pY+npts ) ;
      _Ymax[i] = *std::max_element( pY,pY+npts ) ;
      pYpp = pYp = nullptr ;
      switch ( stype[i] ) {
        case QUINTIC:
          pYpp = baseValue( _npts ) ;
          pYp = baseValue( _npts ) ;
        break ;

        case AKIMA: case BESSEL: case PCHIP: case CUBIC:
          pYp = baseValue( _npts ) ;
        break;

        case CONSTANT: case LINEAR:
        break;
      }
      string h = headers[i] ;
      Spline * & s = splines[sizeType(i)] ;
      switch (stype[i]) {
        case CONSTANT:
          s = new ConstantSpline(h) ;
          static_cast<ConstantSpline*>(s)->reserve_external( sizeType(npts), _X, pY ) ;
          static_cast<ConstantSpline*>(s)->build( _X, pY, sizeType(npts) ) ;
        break;

        case LINEAR:
          s = new LinearSpline(h) ;
          static_cast<LinearSpline*>(s)->reserve_external( sizeType(npts), _X, pY ) ;
          static_cast<LinearSpline*>(s)->build( _X, pY, sizeType(npts) ) ;
        break;

        case AKIMA:
          s = new AkimaSpline(h) ;
          static_cast<AkimaSpline*>(s)->reserve_external( sizeType(npts), _X, pY, pYp ) ;
          static_cast<AkimaSpline*>(s)->build( _X, pY, sizeType(npts) ) ;
        break ;

        case BESSEL:
          s = new BesselSpline(h) ;
          static_cast<BesselSpline*>(s)->reserve_external( sizeType(npts), _X, pY, pYp ) ;
          static_cast<BesselSpline*>(s)->build( _X, pY, sizeType(npts) ) ;
        break ;

        case PCHIP:
          s = new PchipSpline(h) ;
          static_cast<PchipSpline*>(s)->reserve_external( sizeType(npts), _X, pY, pYp ) ;
          static_cast<PchipSpline*>(s)->build( _X, pY, sizeType(npts) ) ;
        break ;

        case CUBIC:
          s = new CubicSpline(h) ;
          static_cast<CubicSpline*>(s)->reserve_external( sizeType(npts), _X, pY, pYp ) ;
          static_cast<CubicSpline*>(s)->build( _X, pY, sizeType(npts) ) ;
        break;

        case QUINTIC:
          s = new QuinticSpline(h) ;
          static_cast<QuinticSpline*>(s)->reserve_external( sizeType(npts), _X, pY, pYp, pYpp ) ;
          static_cast<QuinticSpline*>(s)->build( _X, pY, sizeType(npts) ) ;
        break;
      }
    }
    
  }

}
