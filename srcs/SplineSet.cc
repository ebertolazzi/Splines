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
        case CONSTANT_TYPE: case LINEAR_TYPE:
          mem += npts ;
        break;

        case AKIMA_TYPE: case BESSEL_TYPE: case PCHIP_TYPE: case CUBIC_TYPE:
          mem += 2*npts ; // Y, Yp
        break;

        case QUINTIC_TYPE:
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
        case QUINTIC_TYPE:
          pYpp = baseValue( _npts ) ;
          pYp = baseValue( _npts ) ;
        break ;

        case AKIMA_TYPE: case BESSEL_TYPE: case PCHIP_TYPE: case CUBIC_TYPE:
          pYp = baseValue( _npts ) ;
        break;

        case CONSTANT_TYPE: case LINEAR_TYPE:
        break;

        case SPLINE_SET_TYPE:
        break;
      }
      string h = headers[i] ;
      Spline * & s = splines[sizeType(i)] ;
      switch (stype[i]) {
        case CONSTANT_TYPE:
          s = new ConstantSpline(h) ;
          static_cast<ConstantSpline*>(s)->reserve_external( sizeType(npts), _X, pY ) ;
          static_cast<ConstantSpline*>(s)->build( _X, pY, sizeType(npts) ) ;
        break;

        case LINEAR_TYPE:
          s = new LinearSpline(h) ;
          static_cast<LinearSpline*>(s)->reserve_external( sizeType(npts), _X, pY ) ;
          static_cast<LinearSpline*>(s)->build( _X, pY, sizeType(npts) ) ;
        break;

        case AKIMA_TYPE:
          s = new AkimaSpline(h) ;
          static_cast<AkimaSpline*>(s)->reserve_external( sizeType(npts), _X, pY, pYp ) ;
          static_cast<AkimaSpline*>(s)->build( _X, pY, sizeType(npts) ) ;
        break ;

        case BESSEL_TYPE:
          s = new BesselSpline(h) ;
          static_cast<BesselSpline*>(s)->reserve_external( sizeType(npts), _X, pY, pYp ) ;
          static_cast<BesselSpline*>(s)->build( _X, pY, sizeType(npts) ) ;
        break ;

        case PCHIP_TYPE:
          s = new PchipSpline(h) ;
          static_cast<PchipSpline*>(s)->reserve_external( sizeType(npts), _X, pY, pYp ) ;
          static_cast<PchipSpline*>(s)->build( _X, pY, sizeType(npts) ) ;
        break ;

        case CUBIC_TYPE:
          s = new CubicSpline(h) ;
          static_cast<CubicSpline*>(s)->reserve_external( sizeType(npts), _X, pY, pYp ) ;
          static_cast<CubicSpline*>(s)->build( _X, pY, sizeType(npts) ) ;
        break;

        case QUINTIC_TYPE:
          s = new QuinticSpline(h) ;
          static_cast<QuinticSpline*>(s)->reserve_external( sizeType(npts), _X, pY, pYp, pYpp ) ;
          static_cast<QuinticSpline*>(s)->build( _X, pY, sizeType(npts) ) ;
        break;

        case SPLINE_SET_TYPE:
          SPLINE_ASSERT( false, "SPLINE_SET_TYPE not allowed as spline type\nin SplineSet::build for " << i << "-th spline" ) ;
        break ;
      }
    }
    
  }

  #ifdef SPLINES_USE_GENERIC_CONTAINER
  void
  SplineSet::build( GC::GenericContainer const & gc ) {
    /*
    // gc["headers"]
    // gc["spline_type"]
    // gc["units"]
    // gc["data"]
    // gc["independent"]
    //
    */
    SPLINE_ASSERT( gc.exists("spline_type"), "[" << _name << "] SplineSet::build, missing `spline_type` field!") ;
    SPLINE_ASSERT( gc.exists("headers") ,    "[" << _name << "] SplineSet::build, missing `headers` field!") ;
    SPLINE_ASSERT( gc.exists("data") ,       "[" << _name << "] SplineSet::build, missing `data` field!") ;
  
    GC::GenericContainer const & gc_headers     = gc("headers") ;
    GC::GenericContainer const & gc_spline_type = gc("spline_type") ;
    GC::GenericContainer const & gc_data        = gc("data") ;
    
    SPLINE_ASSERT( GC::GC_VEC_STRING == gc_headers.get_type(),
                   "Field `headers` expected to be of type `vec_string_type` found: ` " <<
                   gc_headers.get_type_name() << "`" ) ;

    SPLINE_ASSERT( GC::GC_VEC_STRING == gc_spline_type.get_type(),
                   "Field `spline_type` expected to be of type `vec_string_type` found: ` " <<
                   gc_spline_type.get_type_name() << "`" ) ;

    SPLINE_ASSERT( GC::GC_MAT_REAL == gc_data.get_type(),
                   "Field `spline_type` expected to be of type `mat_real_type` found: ` " <<
                   gc_data.get_type_name() << "`" ) ;

    indexType independent = 0 ;
    if ( gc.exists("independent") ) {
      GC::GenericContainer const & gc_independent = gc("independent") ;
      SPLINE_ASSERT( GC::GC_INTEGER == gc_independent.get_type(),
                     "Field `independent` expected to be of type `int_type` found: ` " <<
                     gc_independent.get_type_name() << "`" ) ;
      independent = gc_independent.get_int() ;
    }

    GC::mat_real_type   const & data    = gc_data.get_mat_real() ;
    GC::vec_string_type const & headers = gc_headers.get_vec_string() ;

    indexType const nspl = headers.size() ;
    #ifdef SPLINE_USE_ALLOCA
    valueType * X0 = (valueType*)alloca( nkk*sizeof(valueType) ) ;
	  valueType * Y0 = (valueType*)alloca( nkk*sizeof(valueType) ) ;
    char const ** headers_strs = (char const**)alloca( nspl*sizeof(char const *) ) ;
    valueType const ** Y       = (valueType const**)alloca( nspl*sizeof(valueType const *) ) ;
    SplineType * stype         = (SplineType*)alloca( nspl*sizeof(SplineType) ) ;
    #else
    char      const * headers_strs[nspl] ;
    valueType const * Y[nspl] ;
    SplineType        stype[nspl] ;
	  #endif

    for ( indexType i = 0 ; i < nspl ; ++i ) {
      headers_strs[i] = headers[i].c_str() ;
      Y[i] = &data(0,i) ;
      string n = gc_spline_type(i).get_string() ;
      std::transform(n.begin(), n.end(), n.begin(), ::tolower) ;
      SplineType & st = stype[i] ;
      if      ( n == "constant" ) st = CONSTANT_TYPE ;
      else if ( n == "linear"   ) st = LINEAR_TYPE ;
      else if ( n == "akima"    ) st = AKIMA_TYPE ;
      else if ( n == "bessel"   ) st = BESSEL_TYPE ;
      else if ( n == "pchip"    ) st = PCHIP_TYPE ;
      else if ( n == "cubic"    ) st = CUBIC_TYPE ;
      else if ( n == "quintic"  ) st = QUINTIC_TYPE ;
      else {
        SPLINE_ASSERT( false, "[" << _name << "] SplineSet::build\ntype = " << n << " unkonwn for " << i << "-th spline" );
      }
    }

    build( data.numCols(),
           data.numRows(),
           headers_strs,
           stype,
           Y[independent],
           Y,
           nullptr ) ;

  }
  #endif
}
