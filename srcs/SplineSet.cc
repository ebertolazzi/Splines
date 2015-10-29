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
  using std::cbrt ;

  static valueType const machineEps = std::numeric_limits<valueType>::epsilon() ;
  static valueType const m_2pi      = 6.28318530717958647692528676656  ; // 2*pi

  //! quadratic polinomial roots
  /*!
    Compute the roots of the polynomial
    
    \f[ a_0 + a_1 z + a_2 z^2 \f]
    
    and store the results is `real` and `imag`.
    It is assumed that \f$ a_2 \f$ is nonzero.
  */
  /*
    Converted to be compatible with ELF90 by Alan Miller
    amiller @ bigpond.net.au
    WWW-page: http://users.bigpond.net.au/amiller
    Latest revision - 27 February 1997
  */

  static
  int
  quadraticRoots( valueType const a[3],
                  valueType       real[2], 
                  valueType       imag[2] ) {

    // A x^2 + B x + C
    valueType const & C = a[0] ;
    valueType const & B = a[1] ;
    valueType const & A = a[2] ;

    real[0] = real[1] = imag[0] = imag[1] = 0 ;
    
    int res = 0 ;
    if ( a[0] == 0 ) {
      real[1] = -B/A ;
      res = 1 ; // una singola radice reale
    } else {
      valueType twoA = 2*A ;
      valueType d    = B*B - 4*A*C ;
      valueType absd = abs(d) ;
      if ( absd <= 2*machineEps*B*B ) { 
        real[0] = real[1] = -B/twoA ; // EQUAL REAL ROOTS
        res = 3 ; // 2 radici reali coincidenti
      } else {
        valueType r = sqrt(absd) ;
        if ( d < 0 ) { // COMPLEX ROOTS
          real[0] = real[1] = -B/twoA ;
          imag[0] = abs(r/twoA) ;
          imag[1] = -imag[0] ;
          res = 4 ; // 2 radici complesse coniugate
        } else {
          // DISTINCT REAL ROOTS
          if ( B == 0  ) {
            real[0] = abs(r/twoA) ;
            real[1] = -real[0] ;
          } else {
            valueType w = -B ;
            if ( w > 0 ) w += r ; else w -= r ;
            w *= 0.5 ;
            real[0] = C/w ;
            real[1] = w/A ;
          }
          res = 2 ; // 2 radici reali distinte
        }
      }
    }
    return res ;
  }
  
  //! cubic polinomial roots
  /*!
    Compute the roots of the polynomial
    
    \f[ a_0 + a_1 z + a_2 z^2 + a_3 z^3 \f]
    
    and store the results is `real` and `imag`.
    It is assumed that \f$ a_3 \f$ is nonzero.
  */

  static
  int
  cubicRoots( valueType const a[4],
              valueType       real[3], 
              valueType       imag[3] ) {

    // initialize roots
    real[0] = real[1] = real[2] = 
    imag[0] = imag[1] = imag[2] = 0 ;

    // trivial case
    if ( a[0] == 0 ) return quadraticRoots( a+1, real+1, imag+1 ) ; // quadratica degenerata

    // x^3 + A x^2 + B x + C
    valueType const C = a[0]/a[3] ;
    valueType const B = a[1]/a[3] ;
    valueType const A = a[2]/a[3] ;
    
    // p(y-A/3) = y^3 + p*y + q
    valueType const A3 = A/3 ;
    valueType const p  = B-A*A3 ;
    valueType const q  = C+A3*(2*(A3*A3)-B) ;
    
    // scaling equation p(S*z)/S^3 = z^3 + 3*(p/S^2/3)*z + 2*(q/S^3/2)
    valueType const S = max( sqrt(abs(p)), cbrt(abs(q)) ) ;

    // check for a triple root
    if ( S <= machineEps ) {
      real[0] = real[1] = real[2] = -A3 ;
      return 5 ; // 3 radici reali coincidenti
    }
    
    valueType const P     = (p/3)/S/S ;
    valueType const sqrtP = sqrt(abs(p/3))/S ;
    valueType const Q     = (q/2)/S/S/S ;

    valueType const d     = P*P*P + Q*Q ;
    valueType const sqrtd = sqrt(abs(d)) ;

    int res = 0 ;
    if ( sqrtd < abs(q)*machineEps ) {
      // P^3 = - Q^2
      // (x+2*a)(x-a)^2 = x^3 - 3*x*a^2 + 2*a^3
      // cioè -a^2 = P, a^3 = Q ==> a = sqrt(-P)
      valueType tmp = Q > 0 ? sqrtP : -sqrtP ;
      real[1] = real[2] = tmp ;
      real[0] = -2*tmp ;
      res = 6 ; // 3 radici reali, 2 coincidenti
    } else if ( d > 0 ) {
      // w1 = (- Q + sqrt( P^3 + Q^2 ))^(1/3)
      // w2 = (- Q - sqrt( P^3 + Q^2 ))^(1/3)
      valueType w1, w2 ;
      if ( Q > 0 ) {
        w2 = - pow( sqrtd + Q, 1.0 / 3.0 ) ;
        w1 = - P / w2 ;
      } else {
        w1 =   pow( sqrtd - Q, 1.0 / 3.0 ) ;
        w2 = - P / w1 ;
      }
      real[0] = w1 + w2 ;
      real[1] = real[2] = -0.5*real[0] ;
      imag[1] = (w1-w2)*sqrt(3.0/4.0) ;
      imag[2] = -imag[1] ;
      res = 8 ; // 1 reale 2 complesse coniugate
    } else { // 3 radici reali
      // w1 = (- Q + I*sqrt(|P^3 + Q^2|) )^(1/3)
      // w2 = (- Q - I*sqrt(|P^3 + Q^2|) )^(1/3)
      valueType angle  = atan2( sqrtd, -Q ) ;
      if ( angle < 0 ) angle += m_2pi ;
      angle /= 3 ;
      valueType re = sqrtP * cos(angle) ;
      valueType im = sqrtP * sin(angle) ;
      //if ( Q > 0 ) re = -re ;
      real[0]  = 2*re ;
      real[1]  = real[2] = -re ;
      real[1] += sqrt(3.0) * im ;
      real[2] -= sqrt(3.0) * im ;
      res = 7 ; // 3 radici reali distinte
    }

    // scalo radici
    real[0] *= S ; real[1] *= S ; real[2] *= S ;
    imag[0] *= S ; imag[1] *= S ; imag[2] *= S ;
    
    // traslo radici
    real[0] -= A3 ;
    real[1] -= A3 ;
    real[2] -= A3 ;

    return res ;
  }

  void
  SplineSet::info( std::basic_ostream<char> & s ) const {
    s << "SplineSet[" << name() << "] N.points = "
      << _npts << " N.splines = " << _nspl << '\n' ;
    for ( sizeType i = 0 ; i < _nspl ; ++i )
      s << "Spline N." << i << " " << header(i)
        << (is_monotone[i]>0?" is monotone" : " is NOT monotone" )
        << '\n' ;
  }

  void
  SplineSet::build ( indexType  const nspl,
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
    for ( indexType i = 0 ; i < nspl ; ++i ) {
      switch (stype[i]) {
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
          pYp  = baseValue( _npts ) ;
        break ;

        case CUBIC_BASE_TYPE:
          SPLINE_ASSERT( Yp != nullptr,
                         "7th argument of SplineSet::build must be nonnull pointer for `cubic_base` spline" ) ;
          SPLINE_ASSERT( Yp[i] != nullptr,
                         "7th argument Yp[" << i << "] argument of SplineSet::build must be nonnull pointer for `cubic_base` spline" ) ;
          pYp = baseValue( _npts ) ;
          std::copy( Yp[i], Yp[i]+npts, pYp ) ; // copy values
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
      string h = headers[i] ;
      Spline * & s = splines[sizeType(i)] ;
      
      is_monotone[sizeType(i)] = -1 ;
      switch (stype[sizeType(i)]) {
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
              if ( Y[j-1] > Y[j] ) { flag = -1 ; break ; } // non monotone data
              if ( Y[j-1] == Y[j] && X[j-1] < X[j] ) flag = 0 ; // non strict monotone
            }
            is_monotone[sizeType(i)] = flag ;
          }
        break;

        case CUBIC_BASE_TYPE:
          s = new CubicSplineBase(h) ;
          static_cast<CubicSplineBase*>(s)->reserve_external( _npts, _X, pY, pYp ) ;
          static_cast<CubicSplineBase*>(s)->build( _X, pY, pYp, _npts ) ;
          is_monotone[sizeType(i)] = checkCubicSplineMonotonicity( _X, pY, pYp, npts ) ;
        break;

        case CUBIC_TYPE:
          s = new CubicSpline(h) ;
          static_cast<CubicSpline*>(s)->reserve_external( _npts, _X, pY, pYp ) ;
          static_cast<CubicSpline*>(s)->build( _X, pY, _npts ) ;
          is_monotone[sizeType(i)] = checkCubicSplineMonotonicity( _X, pY, pYp, npts ) ;
        break;

        case AKIMA_TYPE:
          s = new AkimaSpline(h) ;
          static_cast<AkimaSpline*>(s)->reserve_external( _npts, _X, pY, pYp ) ;
          static_cast<AkimaSpline*>(s)->build( _X, pY, _npts ) ;
          is_monotone[sizeType(i)] = checkCubicSplineMonotonicity( _X, pY, pYp, npts ) ;
        break ;

        case BESSEL_TYPE:
          s = new BesselSpline(h) ;
          static_cast<BesselSpline*>(s)->reserve_external( _npts, _X, pY, pYp ) ;
          static_cast<BesselSpline*>(s)->build( _X, pY, _npts ) ;
          is_monotone[sizeType(i)] = checkCubicSplineMonotonicity( _X, pY, pYp, npts ) ;
        break ;

        case PCHIP_TYPE:
          s = new PchipSpline(h) ;
          static_cast<PchipSpline*>(s)->reserve_external( _npts, _X, pY, pYp ) ;
          static_cast<PchipSpline*>(s)->build( _X, pY, _npts ) ;
          is_monotone[sizeType(i)] = checkCubicSplineMonotonicity( _X, pY, pYp, npts ) ;
        break ;

        case QUINTIC_TYPE:
          s = new QuinticSpline(h) ;
          static_cast<QuinticSpline*>(s)->reserve_external( _npts, _X, pY, pYp, pYpp ) ;
          static_cast<QuinticSpline*>(s)->build( _X, pY, _npts ) ;
        break;

        case SPLINE_SET_TYPE:
          SPLINE_ASSERT( false, "SPLINE_SET_TYPE not allowed as spline type\nin SplineSet::build for " << i << "-th spline" ) ;
        break ;
        
        default:
          SPLINE_ASSERT( false, "type " << stype[i] << " not recognized as spline type\nin SplineSet::build for " << i << "-th spline" ) ;
        break;
      }
    }
    
  }

  #ifdef SPLINES_USE_GENERIC_CONTAINER
  void
  SplineSet::build( GC::GenericContainer const & gc ) {
    /*
    // gc["headers"]
    // gc["spline_type"]
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

    sizeType const nspl = sizeType(headers.size()) ;
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

    for ( sizeType i = 0 ; i < nspl ; ++i ) {
      headers_strs[i] = headers[i].c_str() ;
      Y[i] = &data(0,i) ;
      string n = gc_spline_type(i).get_string() ;
      std::transform(n.begin(), n.end(), n.begin(), ::tolower) ;
      SplineType & st = stype[i] ;
      if      ( n == "constant"   ) st = CONSTANT_TYPE ;
      else if ( n == "linear"     ) st = LINEAR_TYPE ;
      //else if ( n == "cubic_base" ) st = CUBIC_BASE_TYPE ; NOT yet supported
      else if ( n == "cubic"      ) st = CUBIC_TYPE ;
      else if ( n == "akima"      ) st = AKIMA_TYPE ;
      else if ( n == "bessel"     ) st = BESSEL_TYPE ;
      else if ( n == "pchip"      ) st = PCHIP_TYPE ;
      else if ( n == "quintic"    ) st = QUINTIC_TYPE ;
      else {
        SPLINE_ASSERT( false, "[" << _name << "] SplineSet::build\ntype = " << n << " unkonwn for " << i << "-th spline" );
      }
    }

    build( indexType(data.numCols()),
           indexType(data.numRows()),
           headers_strs,
           stype,
           Y[independent],
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
  SplineSet::intersect( indexType spl, valueType zeta, valueType & x ) const {
    SPLINE_ASSERT( spl >= 0 && spl < _nspl,
                   "Spline n." << spl << " is not in SplineSet") ;
    SPLINE_ASSERT( is_monotone[sizeType(spl)]>0,
                   "Spline n." << spl << " is not monotone") ;
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
    if ( S->type() == LINEAR_TYPE ) {
      x = a + (b-a)*(zeta-ya)/(yb-ya) ;
    } else {
      valueType const * dX = _Yp[sizeType(spl)] ;
      valueType dya = dX[interval] ;
      valueType dyb = dX[interval+1] ;
      valueType const coeffs[4] = { ya-zeta, dya, (3*DY/DX-2*dya-dyb)/DX, (dyb+dya-2*DY/DX)/(DX*DX) } ;
      valueType real[3], imag[3] ;
      int icase = cubicRoots( coeffs, real, imag ) ;
      cout << "coeffs " << coeffs[0] << " " << coeffs[1] << " " << coeffs[2] << " " << coeffs[3] << '\n' ;
      cout << "DX " << DX << " zeta = " << zeta << '\n' ;
      int nr = 0 ;
      switch ( icase ) {
      case 4:// 2 radici complesse coniugate
        break ;
      case 1:// una singola radice reale
      case 3:// 2 radici reali coincidenti
      case 5:// 3 radici reali coincidenti
      case 8:// 1 reale 2 complesse coniugate
        nr = 1 ;
        break ;
      case 2:// 2 radici reali distinte
      case 6:// 3 radici reali, 2 coincidenti
        nr = 2 ;
        break ;
      case 7:// 3 radici reali distinte
        nr = 3 ;
        break ;
      }
      // cerca radice buona
      bool ok = false ;
      for ( indexType i = 0 ; i < nr && !ok ; ++i ) {
        ok = real[i] >= 0 && real[i] <= DX ;
        if ( ok ) x = a + real[i] ;
      }
      SPLINE_ASSERT( ok, "SplineSet, failed to find intersection with independent spline at zeta = " << zeta << " icase = " << icase ) ;
    }
    return S ;
  }

  void
  SplineSet::eval2( indexType           spl,
                    valueType           zeta,
                    vector<valueType> & vals ) const {
    valueType x ;
    intersect( spl, zeta, x ) ;
    vals.resize(_nspl) ;
    for ( sizeType i = 0 ; i < _nspl ; ++i )
      vals[i] = (*splines[i])(x) ;
  }

  void
  SplineSet::eval2( indexType spl,
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
  SplineSet::eval2_D( indexType           spl,
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
  SplineSet::eval2_D( indexType spl,
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
  SplineSet::eval2_DD( indexType           spl,
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
  SplineSet::eval2_DD( indexType spl,
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
  SplineSet::eval2_DDD( indexType           spl,
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
  SplineSet::eval2_DDD( indexType spl,
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
}
