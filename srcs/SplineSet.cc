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
#include <limits>
#include <cmath>

/**
 * 
 */

namespace Splines {

  using std::abs ;
  using std::sqrt ;

  void
  SplineSet::info( std::basic_ostream<char> & s ) const {
    s << "SplineSet[" << name() << "] N.points = "
      << _npts << " N.splines = " << _nspl << '\n' ;
    for ( sizeType i = 0 ; i < _nspl ; ++i ) {
      s << "Spline N." << i ;
      switch ( is_monotone[i] ) {
        case -2: s << " with NON monotone data\n" ; break ;
        case -1: s << " is NOT monotone\n"        ; break ;
        case  0: s << " is monotone\n"            ; break ;
        case  1: s << " is strictly monotone\n"   ; break ;
        default: SPLINE_ASSERT( false, "SplineSet::info classification: " << is_monotone[i] << " not in range {-2,-1,0,1}" ) ;
      }
      splines[i]->info(s) ;
    }
  }

  sizeType
  SplineSet::getPosition( char const * hdr ) const {
    map<string,sizeType>::const_iterator it = header_to_position.find(hdr) ;
    SPLINE_ASSERT( it != header_to_position.end(), "Spline [" << hdr << "] not found!" ) ;
    return it->second ;
  }

  void
  SplineSet::build( indexType  const nspl,
                    indexType  const npts,
                    char       const *headers[],
                    SplineType const stype[],
                    valueType  const X[],
                    valueType  const *Y[],
                    valueType  const *Yp[] ) {
    SPLINE_ASSERT( nspl > 0, "SplineSet::build expected positive nspl = " << nspl ) ;
    SPLINE_ASSERT( npts > 1, "SplineSet::build expected npts = " << npts << " greather than 1" ) ;
    _nspl = sizeType(nspl) ;
    _npts = sizeType(npts) ;
    // allocate memory
    splines.resize(_nspl) ;
    is_monotone.resize(_nspl) ;
    indexType mem = npts ;
    for ( sizeType spl = 0 ; spl < nspl ; ++spl ) {
      switch (stype[spl]) {
        case CONSTANT_TYPE:
        case LINEAR_TYPE:
          mem += npts ;
        break;

        case CUBIC_BASE_TYPE:
        case CUBIC_TYPE:
        case AKIMA_TYPE:
        case BESSEL_TYPE:
        case PCHIP_TYPE:
          mem += 2*npts ; // Y, Yp
        break;

        case QUINTIC_TYPE:
          mem += 3*npts ; // Y, Yp, Ypp
        break;
        
        default:
          SPLINE_ASSERT( false,
                         "SplineSet::build\nAt spline N. " << spl <<
                         " named " << headers[spl] <<
                         " unknwn type = " << stype[spl] ) ;
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
    for ( sizeType spl = 0 ; spl < nspl ; ++spl ) {
      valueType *& pY   = _Y[spl] ;
      valueType *& pYp  = _Yp[spl] ;
      valueType *& pYpp = _Ypp[spl] ;
      pY = baseValue( _npts ) ;
      std::copy( Y[spl], Y[spl]+npts, pY ) ;
      if ( stype[spl] == CONSTANT_TYPE ) {
        _Ymin[spl] = *std::min_element( pY,pY+npts-1 ) ;
        _Ymax[spl] = *std::max_element( pY,pY+npts-1 ) ;        
      } else {
        _Ymin[spl] = *std::min_element( pY,pY+npts ) ;
        _Ymax[spl] = *std::max_element( pY,pY+npts ) ;
      }
      pYpp = pYp = nullptr ;
      switch ( stype[spl] ) {
        case QUINTIC_TYPE:
          pYpp = baseValue( _npts ) ;
          pYp  = baseValue( _npts ) ;
        break ;

        case CUBIC_BASE_TYPE:
          SPLINE_ASSERT( Yp != nullptr,
                         "SplineSet::build\nAt spline N. " << spl << " named " << headers[spl] <<
                         "\n7th argument of SplineSet::build must be nonnull pointer for `cubic_base` spline" ) ;
          SPLINE_ASSERT( Yp[spl] != nullptr,
                         "SplineSet::build\nAt spline N. " << spl << " named " << headers[spl] <<
                         "\n7th argument Yp[" << spl << "] argument of SplineSet::build must be nonnull pointer for `cubic_base` spline" ) ;
          pYp = baseValue( _npts ) ;
          std::copy( Yp[spl], Yp[spl]+npts, pYp ) ; // copy values
        break;

        case CUBIC_TYPE:
        case AKIMA_TYPE:
        case BESSEL_TYPE:
        case PCHIP_TYPE:
          pYp = baseValue( _npts ) ;
        break;

        case CONSTANT_TYPE:
        case LINEAR_TYPE:
        break;

        case SPLINE_SET_TYPE:
        break;
      }
      string h = headers[spl] ;
      Spline * & s = splines[spl] ;
      
      is_monotone[spl] = -1 ;
      switch (stype[spl]) {
        case CONSTANT_TYPE:
          s = new ConstantSpline(h) ;
          static_cast<ConstantSpline*>(s)->reserve_external( _npts, _X, pY ) ;
          static_cast<ConstantSpline*>(s)->build( _X, pY, _npts ) ;
        break;

        case LINEAR_TYPE:
          s = new LinearSpline(h) ;
          static_cast<LinearSpline*>(s)->reserve_external( _npts, _X, pY ) ;
          static_cast<LinearSpline*>(s)->build( _X, pY, _npts ) ;
          // check monotonicity of data
          { indexType flag = 1 ;
            for ( sizeType j = 1 ; j < _npts ; ++j ) {
              if ( pY[j-1] > pY[j] ) { flag = -1 ; break ; } // non monotone data
              if ( pY[j-1] == pY[j] && _X[j-1] < _X[j] ) flag = 0 ; // non strict monotone
            }
            is_monotone[spl] = flag ;
          }
        break;

        case CUBIC_BASE_TYPE:
          s = new CubicSplineBase(h) ;
          static_cast<CubicSplineBase*>(s)->reserve_external( _npts, _X, pY, pYp ) ;
          static_cast<CubicSplineBase*>(s)->build( _X, pY, pYp, _npts ) ;
          is_monotone[spl] = checkCubicSplineMonotonicity( _X, pY, pYp, npts ) ;
        break;

        case CUBIC_TYPE:
          s = new CubicSpline(h) ;
          static_cast<CubicSpline*>(s)->reserve_external( _npts, _X, pY, pYp ) ;
          static_cast<CubicSpline*>(s)->build( _X, pY, _npts ) ;
          is_monotone[spl] = checkCubicSplineMonotonicity( _X, pY, pYp, npts ) ;
        break;

        case AKIMA_TYPE:
          s = new AkimaSpline(h) ;
          static_cast<AkimaSpline*>(s)->reserve_external( _npts, _X, pY, pYp ) ;
          static_cast<AkimaSpline*>(s)->build( _X, pY, _npts ) ;
          is_monotone[spl] = checkCubicSplineMonotonicity( _X, pY, pYp, npts ) ;
        break ;

        case BESSEL_TYPE:
          s = new BesselSpline(h) ;
          static_cast<BesselSpline*>(s)->reserve_external( _npts, _X, pY, pYp ) ;
          static_cast<BesselSpline*>(s)->build( _X, pY, _npts ) ;
          is_monotone[spl] = checkCubicSplineMonotonicity( _X, pY, pYp, npts ) ;
        break ;

        case PCHIP_TYPE:
          s = new PchipSpline(h) ;
          static_cast<PchipSpline*>(s)->reserve_external( _npts, _X, pY, pYp ) ;
          static_cast<PchipSpline*>(s)->build( _X, pY, _npts ) ;
          is_monotone[spl] = checkCubicSplineMonotonicity( _X, pY, pYp, npts ) ;
        break ;

        case QUINTIC_TYPE:
          s = new QuinticSpline(h) ;
          static_cast<QuinticSpline*>(s)->reserve_external( _npts, _X, pY, pYp, pYpp ) ;
          static_cast<QuinticSpline*>(s)->build( _X, pY, _npts ) ;
        break;

        case SPLINE_SET_TYPE:
          SPLINE_ASSERT( false,
                         "SplineSet::build\nAt spline N. " << spl << " named " << headers[spl] <<
                         "\nSPLINE_SET_TYPE not allowed as spline type\nin SplineSet::build for " << spl << "-th spline" ) ;
        break ;
        
        default:
          SPLINE_ASSERT( false,
                         "SplineSet::build\nAt spline N. " << spl << " named " << headers[spl] <<
                         "\ntype " << stype[spl] << " not recognized as spline type\nin SplineSet::build for " << spl << "-th spline" ) ;
        break;
      }
      header_to_position[s->name()] = spl ;
    }
    
  }

  #ifdef SPLINES_USE_GENERIC_CONTAINER

  using GenericContainerNamepace::GC_VEC_REAL ;
  using GenericContainerNamepace::GC_VEC_STRING ;
  using GenericContainerNamepace::GC_MAT_REAL ;
  using GenericContainerNamepace::GC_VECTOR ;
  using GenericContainerNamepace::GC_INTEGER ;
  using GenericContainerNamepace::mat_real_type ;
  using GenericContainerNamepace::vec_real_type ;
  using GenericContainerNamepace::vec_string_type ;
  using GenericContainerNamepace::vector_type ;

  void
  SplineSet::setup( GenericContainer const & gc ) {
    /*
    // gc["headers"]
    // gc["spline_type"]
    // gc["data"]
    // gc["independent"]
    //
    */
    SPLINE_ASSERT( gc.exists("spline_type"), "[SplineSet[" << _name << "]::setup] missing `spline_type` field!") ;
    SPLINE_ASSERT( gc.exists("headers") ,    "[SplineSet[" << _name << "]::setup] missing `headers` field!") ;
    SPLINE_ASSERT( gc.exists("xdata") ,      "[SplineSet[" << _name << "]::setup] missing `xdata` field!") ;
    SPLINE_ASSERT( gc.exists("ydata") ,      "[SplineSet[" << _name << "]::setup] missing `ydata` field!") ;
  
    GenericContainer const & gc_headers     = gc("headers") ;
    GenericContainer const & gc_spline_type = gc("spline_type") ;
    GenericContainer const & gc_ydata       = gc("ydata") ;
    GenericContainer const & gc_xdata       = gc("xdata") ;

    SPLINE_ASSERT( GC_VEC_STRING == gc_headers.get_type(),
                   "[SplineSet[" << _name << "]::setup] field `headers` expected to be of type `vec_string_type` found: ` " <<
                   gc_headers.get_type_name() << "`" ) ;

    SPLINE_ASSERT( GC_VEC_STRING == gc_spline_type.get_type(),
                   "[SplineSet[" << _name << "]::setup] field `spline_type` expected to be of type `vec_string_type` found: ` " <<
                   gc_spline_type.get_type_name() << "`" ) ;

    SPLINE_ASSERT( GC_VEC_REAL == gc_xdata.get_type(),
                   "[SplineSet[" << _name << "]::setup] field `xdata` expected to be of type `vec_real` found: ` " <<
                   gc_xdata.get_type_name() << "`" ) ;

    vec_real_type   const & X       = gc_xdata.get_vec_real() ;
    vec_string_type const & headers = gc_headers.get_vec_string() ;

    sizeType const nspl = sizeType(headers.size()) ;
    #ifdef SPLINE_USE_ALLOCA
    char const ** headers_strs = (char const**)alloca( nspl*sizeof(char const *) ) ;
    valueType const ** Y       = (valueType const**)alloca( nspl*sizeof(valueType const *) ) ;
    SplineType * stype         = (SplineType*)alloca( nspl*sizeof(SplineType) ) ;
    #else
    char      const * headers_strs[nspl] ;
    valueType const * Y[nspl] ;
    SplineType        stype[nspl] ;
	  #endif

    sizeType nrow ;
    if ( GC_MAT_REAL == gc_ydata.get_type() ) {
      mat_real_type const & data = gc_ydata.get_mat_real() ;
      SPLINE_ASSERT( nspl == data.numCols(),
                     "[SplineSet[" << _name << "]::setup] number of headers [" << nspl <<
                     "] differs the numeber of columns [" << data.numCols() << "] in data" ) ;
      for ( sizeType i = 0 ; i < nspl ; ++i ) Y[i] = &data(0,i) ;
      nrow = data.numRows() ;
    } else if ( GC_VECTOR == gc_ydata.get_type() ) {
      nrow = 0 ;
      vector_type const & data = gc_ydata.get_vector() ;
      SPLINE_ASSERT( data.size() == nspl,
                    "[SplineSet[" << _name << "]::setup] field `data` expected of size " << nspl <<
                    " found of size " << data.size()  ) ;
      for ( sizeType i = 0 ; i < nspl ; ++i ) {
        GenericContainer const & datai = data[i] ;
        SPLINE_ASSERT( GC_VEC_REAL == datai.get_type(),
                       "[SplineSet::setup] data[" << i << "] expected to be of type `vec_real_type` found: ` " <<
                        datai.get_type_name() << "`" ) ;
        Y[i] = &datai.get_vec_real().front() ;
        sizeType len = sizeType(datai.get_vec_real().size()) ;
        if ( i == 0 ) {
          nrow = len ;
        } else {
          if ( stype[i] == CONSTANT_TYPE ) ++len ; // constant spline uses n-1 points
          if ( nrow > len ) nrow = len ;
        }
      }
    } else {
      SPLINE_ASSERT( false,
                     "[SplineSet[" << _name << "]::setup] field `data` expected to be of type `mat_real_type` or `vector_type` found: ` " <<
                      gc_ydata.get_type_name() << "`" ) ;
    }

    for ( sizeType i = 0 ; i < nspl ; ++i ) {
      headers_strs[i] = headers[i].c_str() ;
      string n = gc_spline_type.get_string(i) ;
      std::transform(n.begin(), n.end(), n.begin(), ::tolower) ;
      SplineType & st = stype[i] ;
      if      ( n == spline_type[CONSTANT_TYPE]   ) st = CONSTANT_TYPE ;
      else if ( n == spline_type[LINEAR_TYPE]     ) st = LINEAR_TYPE ;
      //else if ( n == spline_type[CUBIC_BASE_TYPE] ) st = CUBIC_BASE_TYPE ; NOT YET IMPLEMENTED
      else if ( n == spline_type[CUBIC_TYPE]      ) st = CUBIC_TYPE ;
      else if ( n == spline_type[AKIMA_TYPE]      ) st = AKIMA_TYPE ;
      else if ( n == spline_type[BESSEL_TYPE]     ) st = BESSEL_TYPE ;
      else if ( n == spline_type[PCHIP_TYPE]      ) st = PCHIP_TYPE ;
      else if ( n == spline_type[QUINTIC_TYPE]    ) st = QUINTIC_TYPE ;
      else {
        SPLINE_ASSERT( false, "[" << _name << "] SplineSet::build\ntype = " << n << " unkonwn for " << i << "-th spline" );
      }
    }

    build( nspl,
           nrow,
           headers_strs,
           stype,
           &X.front(),
           Y,
           nullptr ) ;

  }
  #endif

  // vectorial values
  void
  SplineSet::getHeaders( vector<string> & h ) const {
    h.resize(_nspl) ;
    for ( sizeType i = 0 ; i < _nspl ; ++i )
      h[i] = splines[i]->name() ;
  }

  void
  SplineSet::eval( valueType x, vector<valueType> & vals ) const {
    vals.resize(_nspl) ;
    for ( sizeType i = 0 ; i < _nspl ; ++i )
      vals[i] = (*splines[i])(x) ;
  }

  void
  SplineSet::eval( valueType x, valueType vals[], indexType incy ) const {
    sizeType ii = 0 ;
    for ( sizeType i = 0 ; i < _nspl ; ++i, ii += sizeType(incy) )
      vals[ii] = (*splines[i])(x) ;
  }

  void
  SplineSet::eval_D( valueType x, vector<valueType> & vals ) const {
    vals.resize(_nspl) ;
    for ( sizeType i = 0 ; i < _nspl ; ++i )
      vals[i] = splines[i]->D(x) ;
  }

  void
  SplineSet::eval_D( valueType x, valueType vals[], indexType incy ) const {
    sizeType ii = 0 ;
    for ( sizeType i = 0 ; i < _nspl ; ++i, ii += sizeType(incy) )
      vals[ii] = splines[i]->D(x) ;
  }

  void
  SplineSet::eval_DD( valueType x, vector<valueType> & vals ) const {
    vals.resize(_nspl) ;
    for ( sizeType i = 0 ; i < _nspl ; ++i )
      vals[i] = splines[i]->DD(x) ;
  }

  void
  SplineSet::eval_DD( valueType x, valueType vals[], indexType incy ) const {
    sizeType ii = 0 ;
    for ( sizeType i = 0 ; i < _nspl ; ++i, ii += sizeType(incy) )
      vals[ii] = splines[i]->DD(x) ;
  }

  void
  SplineSet::eval_DDD( valueType x, vector<valueType> & vals ) const {
    vals.resize(_nspl) ;
    for ( sizeType i = 0 ; i < _nspl ; ++i )
      vals[i] = splines[i]->DDD(x) ;
  }

  void
  SplineSet::eval_DDD( valueType x, valueType vals[], indexType incy ) const {
    sizeType ii = 0 ;
    for ( sizeType i = 0 ; i < _nspl ; ++i, ii += sizeType(incy) )
      vals[ii] = splines[i]->DDD(x) ;
  }

  // vectorial values
  
  Spline const *
  SplineSet::intersect( sizeType spl, valueType zeta, valueType & x ) const {
    SPLINE_ASSERT( spl >= 0 && spl < _nspl,
                   "Spline n." << spl << " is not in SplineSet") ;
    SPLINE_ASSERT( is_monotone[sizeType(spl)]>0,
                   "Spline n." << spl << " is not monotone and can't be used as independent") ;
    Spline const * S = splines[sizeType(spl)] ;
    // cerco intervallo intersezione
    valueType const * X = _Y[sizeType(spl)] ;
    SPLINE_ASSERT( zeta >= X[0] && zeta <= X[_npts-1],
                   "SplineSet, evaluation at zeta = " << zeta <<
                   " is out of range: [" << X[0] << ", " << X[_npts-1] << "]" ) ;

    sizeType interval = sizeType(lower_bound( X, X+_npts, zeta ) - X) ;
    if ( interval > 0 ) --interval ;
    if ( X[interval] == X[interval+1] ) ++interval ; // degenerate interval for duplicated nodes
    if ( interval >= _npts-1 ) interval = _npts-2 ;

    // compute intersection
    valueType a  = _X[interval] ;
    valueType b  = _X[interval+1] ;
    valueType ya = X[interval] ;
    valueType yb = X[interval+1] ;
    valueType DX = b-a ;
    valueType DY = yb-ya ;
    SPLINE_ASSERT( zeta >= ya && zeta <= yb,
                   "SplineSet, Bad interval [ " << ya << "," << yb << "] for zeta = " << zeta ) ;
    SPLINE_ASSERT( a < b,
                   "SplineSet, Bad x interval [ " << a << "," << b << "]" ) ;
    if ( S->type() == LINEAR_TYPE ) {
      x = a + (b-a)*(zeta-ya)/(yb-ya) ;
    } else {
      valueType const * dX = _Yp[sizeType(spl)] ;
      valueType        dya = dX[interval] ;
      valueType        dyb = dX[interval+1] ;
      valueType coeffs[4] = { ya-zeta, dya, (3*DY/DX-2*dya-dyb)/DX, (dyb+dya-2*DY/DX)/(DX*DX) } ;
      valueType real[3], imag[3] ;
      pair<int,int> icase = cubicRoots( coeffs, real, imag ) ;
      SPLINE_ASSERT( icase.first > 0,
                     "SplineSet, No intersection found with independent spline at zeta = " << zeta ) ;
      // cerca radice buona
      bool ok = false ;
      for ( indexType i = 0 ; i < icase.first && !ok ; ++i ) {
        ok = real[i] >= 0 && real[i] <= DX ;
        if ( ok ) x = a + real[i] ;
      }
      SPLINE_ASSERT( ok, "SplineSet, failed to find intersection with independent spline at zeta = " << zeta ) ;
    }
    return S ;
  }

  void
  SplineSet::eval2( sizeType            spl,
                    valueType           zeta,
                    vector<valueType> & vals ) const {
    valueType x ;
    intersect( spl, zeta, x ) ;
    vals.resize(_nspl) ;
    for ( sizeType i = 0 ; i < _nspl ; ++i )
      vals[i] = (*splines[i])(x) ;
  }

  void
  SplineSet::eval2( sizeType  spl,
                    valueType zeta,
                    valueType vals[],
                    indexType incy ) const {
    valueType x ;
    intersect( spl, zeta, x ) ;
    sizeType ii = 0 ;
    for ( sizeType i = 0 ; i < _nspl ; ++i, ii += sizeType(incy) )
      vals[ii] = (*splines[i])(x) ;
  }

  void
  SplineSet::eval2_D( sizeType            spl,
                      valueType           zeta,
                      vector<valueType> & vals ) const {
    valueType x ;
    Spline const * S = intersect( spl, zeta, x ) ;
    valueType ds = S->D(x) ;
    vals.resize(_nspl) ;
    for ( sizeType i = 0 ; i < _nspl ; ++i )
      vals[i] = splines[i]->D(x)/ds ;
  }

  void
  SplineSet::eval2_D( sizeType  spl,
                      valueType zeta,
                      valueType vals[],
                      indexType incy ) const {
    valueType x ;
    Spline const * S = intersect( spl, zeta, x ) ;
    valueType ds = S->D(x) ;
    sizeType ii = 0 ;
    for ( sizeType i = 0 ; i < _nspl ; ++i, ii += sizeType(incy) )
      vals[ii] = splines[i]->D(x)/ds ;
  }

  void
  SplineSet::eval2_DD( sizeType            spl,
                       valueType           zeta,
                       vector<valueType> & vals ) const {
    valueType x ;
    Spline const * S = intersect( spl, zeta, x ) ;
    valueType dt  = 1/S->D(x) ;
    valueType dt2 = dt*dt ;
    valueType ddt = -S->DD(x)*(dt*dt2) ;
    vals.resize(_nspl) ;
    for ( sizeType i = 0 ; i < _nspl ; ++i )
      vals[i] = splines[i]->DD(x)*dt2+splines[i]->D(x)*ddt ;
  }

  void
  SplineSet::eval2_DD( sizeType  spl,
                       valueType zeta,
                       valueType vals[],
                       indexType incy ) const {
    valueType x ;
    Spline const * S = intersect( spl, zeta, x ) ;
    valueType dt  = 1/S->D(x) ;
    valueType dt2 = dt*dt ;
    valueType ddt = -S->DD(x)*(dt*dt2) ;
    sizeType ii = 0 ;
    for ( sizeType i = 0 ; i < _nspl ; ++i, ii += sizeType(incy) )
      vals[ii] = splines[i]->DD(x)*dt2+splines[i]->D(x)*ddt ;
  }

  void
  SplineSet::eval2_DDD( sizeType            spl,
                        valueType           zeta,
                        vector<valueType> & vals ) const {
    valueType x ;
    Spline const * S = intersect( spl, zeta, x ) ;
    valueType dt  = 1/S->D(x) ;
    valueType dt3 = dt*dt*dt ;
    valueType ddt = -S->DD(x)*dt3 ;
    valueType dddt = 3*(ddt*ddt)/dt-S->DDD(x)*(dt*dt3) ;
    vals.resize(_nspl) ;
    for ( sizeType i = 0 ; i < _nspl ; ++i )
      vals[i] = splines[i]->DDD(x)*dt3 +
                3*splines[i]->DD(x)*dt*ddt +
                splines[i]->D(x)*dddt ;
  }

  void
  SplineSet::eval2_DDD( sizeType  spl,
                        valueType zeta,
                        valueType vals[],
                        indexType incy ) const {
    valueType x ;
    Spline const * S = intersect( spl, zeta, x ) ;
    valueType dt  = 1/S->D(x) ;
    valueType dt3 = dt*dt*dt ;
    valueType ddt = -S->DD(x)*dt3 ;
    valueType dddt = 3*(ddt*ddt)/dt-S->DDD(x)*(dt*dt3) ;
    sizeType ii = 0 ;
    for ( sizeType i = 0 ; i < _nspl ; ++i, ii += sizeType(incy) )
      vals[ii] = splines[i]->DDD(x)*dt3 +
                 3*splines[i]->DD(x)*dt*ddt +
                 splines[i]->D(x)*dddt ;
  }

  valueType
  SplineSet::eval2( valueType    zeta,
                    char const * indep,
                    char const * name ) const {
    vector<valueType> vals ;
    eval2( getPosition(indep), zeta, vals ) ;
    return vals[unsigned(getPosition(name))] ;
  }

  valueType
  SplineSet::eval2_D( valueType    zeta,
                      char const * indep,
                      char const * name ) const {
    vector<valueType> vals ;
    eval2_D( getPosition(indep), zeta, vals ) ;
    return vals[unsigned(getPosition(name))] ;
  }

  valueType
  SplineSet::eval2_DD( valueType    zeta,
                       char const * indep,
                       char const * name ) const {
    vector<valueType> vals ;
    eval2_DD( getPosition(indep), zeta, vals ) ;
    return vals[unsigned(getPosition(name))] ;
  }

  valueType
  SplineSet::eval2_DDD( valueType    zeta,
                        char const * indep,
                        char const * name ) const {
    vector<valueType> vals ;
    eval2_DDD( getPosition(indep), zeta, vals ) ;
    return vals[unsigned(getPosition(name))] ;
  }

  valueType
  SplineSet::eval2( valueType zeta, sizeType indep, sizeType spl ) const {
    vector<valueType> vals ;
    eval2( indep, zeta, vals ) ;
    return vals[spl] ;
  }

  valueType
  SplineSet::eval2_D( valueType zeta, sizeType indep, sizeType spl ) const {
    vector<valueType> vals ;
    eval2_D( indep, zeta, vals ) ;
    return vals[spl] ;
  }

  valueType
  SplineSet::eval2_DD( valueType zeta, sizeType indep, sizeType spl ) const {
    vector<valueType> vals ;
    eval2_DD( indep, zeta, vals ) ;
    return vals[spl] ;
  }

  valueType
  SplineSet::eval2_DDD( valueType zeta, sizeType indep, sizeType spl ) const {
    vector<valueType> vals ;
    eval2_DDD( indep, zeta, vals ) ;
    return vals[spl] ;
  }
}
