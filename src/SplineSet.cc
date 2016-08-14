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
#include <limits>
#include <cmath>

/**
 * 
 */

namespace Splines {

  using std::abs ;
  using std::sqrt ;

  //! spline constructor
  SplineSet::SplineSet( string const & name )
  : _name(name)
  , baseValue(name+"_values")
  , basePointer(name+"_pointers")
  , _npts(0)
  , _nspl(0)
  , _X(nullptr)
  , _Y(nullptr)
  , _Yp(nullptr)
  , _Ypp(nullptr)
  , _Ymin(nullptr)
  , _Ymax(nullptr)
  , lastInterval(0)
  {}

  //! spline destructor
  SplineSet::~SplineSet() {
    baseValue.free() ;
    basePointer.free() ;
  }

  void
  SplineSet::info( std::basic_ostream<char> & s ) const {
    s << "SplineSet[" << name() << "] n.points = "
      << _npts << " n.splines = " << _nspl << '\n' ;
    for ( sizeType i = 0 ; i < _nspl ; ++i ) {
      s << "\nSpline n." << i ;
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

  void
  SplineSet::dump_table( std::basic_ostream<char> & stream, sizeType num_points ) const {
    vector<valueType> vals ;
    stream << 's' ;
    for ( sizeType i = 0 ; i < numSplines() ; ++i ) stream << '\t' << header(i) ;
    stream << '\n' ;
  
    for ( sizeType j = 0 ; j < num_points ; ++j ) {
      valueType s = xMin() + ((xMax()-xMin())*j)/(num_points-1) ;
      this->eval( s, vals ) ;
      stream << s ;
      for ( sizeType i = 0 ; i < numSplines() ; ++i ) stream << '\t' << vals[i] ;
      stream << '\n' ;
    }
  }

  sizeType
  SplineSet::getPosition( char const * hdr ) const {
    map<string,sizeType>::const_iterator it = header_to_position.find(hdr) ;
    SPLINE_ASSERT( it != header_to_position.end(), "SplineSet::getPosition(\"" << hdr << "\") not found!" ) ;
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
      case QUINTIC_TYPE:
        mem += npts ; // Y, Yp, Ypp
      case CUBIC_BASE_TYPE:
      case CUBIC_TYPE:
      case AKIMA_TYPE:
      case BESSEL_TYPE:
      case PCHIP_TYPE:
        mem += npts ; // Y, Yp
      case CONSTANT_TYPE:
      case LINEAR_TYPE:
        mem += npts ;
      break;
      case SPLINE_SET_TYPE:
      case BSPLINE_TYPE:
        SPLINE_ASSERT( false,
                       "SplineSet::build\nAt spline n. " << spl <<
                       " named " << headers[spl] <<
                       " cannot be done for type = " << stype[spl] ) ;
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
      case CUBIC_BASE_TYPE:
        SPLINE_ASSERT( Yp != nullptr,
                       "SplineSet::build\nAt spline n. " << spl << " named " << headers[spl] <<
                       "\n7th argument of SplineSet::build must be nonnull pointer for `cubic_base` spline" ) ;
        SPLINE_ASSERT( Yp[spl] != nullptr,
                       "SplineSet::build\nAt spline n. " << spl << " named " << headers[spl] <<
                       "\n7th argument Yp[" << spl << "] argument of SplineSet::build must be nonnull pointer for `cubic_base` spline" ) ;
        pYp = baseValue( _npts ) ;
        std::copy( Yp[spl], Yp[spl]+npts, pYp ) ; // copy values
        break;

      case QUINTIC_TYPE:
        pYpp = baseValue( _npts ) ;
      case CUBIC_TYPE:
      case AKIMA_TYPE:
      case BESSEL_TYPE:
      case PCHIP_TYPE:
        pYp = baseValue( _npts ) ;
      case CONSTANT_TYPE:
      case LINEAR_TYPE:
      case SPLINE_SET_TYPE:
        break;
      case BSPLINE_TYPE: // nothing to do (not implemented)
        break ;
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
        is_monotone[spl] = checkCubicSplineMonotonicity( _X, pY, pYp, _npts ) ;
        break;

      case CUBIC_TYPE:
        s = new CubicSpline(h) ;
        static_cast<CubicSpline*>(s)->reserve_external( _npts, _X, pY, pYp ) ;
        static_cast<CubicSpline*>(s)->build( _X, pY, _npts ) ;
        is_monotone[spl] = checkCubicSplineMonotonicity( _X, pY, pYp, _npts ) ;
        break;

      case AKIMA_TYPE:
        s = new AkimaSpline(h) ;
        static_cast<AkimaSpline*>(s)->reserve_external( _npts, _X, pY, pYp ) ;
        static_cast<AkimaSpline*>(s)->build( _X, pY, _npts ) ;
        is_monotone[spl] = checkCubicSplineMonotonicity( _X, pY, pYp, _npts ) ;
        break ;

      case BESSEL_TYPE:
        s = new BesselSpline(h) ;
        static_cast<BesselSpline*>(s)->reserve_external( _npts, _X, pY, pYp ) ;
        static_cast<BesselSpline*>(s)->build( _X, pY, _npts ) ;
        is_monotone[spl] = checkCubicSplineMonotonicity( _X, pY, pYp, _npts ) ;
        break ;

      case PCHIP_TYPE:
        s = new PchipSpline(h) ;
        static_cast<PchipSpline*>(s)->reserve_external( _npts, _X, pY, pYp ) ;
        static_cast<PchipSpline*>(s)->build( _X, pY, _npts ) ;
        is_monotone[spl] = checkCubicSplineMonotonicity( _X, pY, pYp, _npts ) ;
        break ;

      case QUINTIC_TYPE:
        s = new QuinticSpline(h) ;
        static_cast<QuinticSpline*>(s)->reserve_external( _npts, _X, pY, pYp, pYpp ) ;
        static_cast<QuinticSpline*>(s)->build( _X, pY, _npts ) ;
        break;

      case SPLINE_SET_TYPE:
        SPLINE_ASSERT( false,
                       "SplineSet::build\nAt spline n. " << spl << " named " << headers[spl] <<
                       "\nSPLINE_SET_TYPE not allowed as spline type\nin SplineSet::build for " << spl << "-th spline" ) ;
      case BSPLINE_TYPE:
        SPLINE_ASSERT( false,
                       "SplineSet::build\nAt spline n. " << spl << " named " << headers[spl] <<
                       "\nBSPLINE_TYPE not allowed as spline type\nin SplineSet::build for " << spl << "-th spline" ) ;
      }
      header_to_position[s->name()] = spl ;
    }

  }

  #ifdef SPLINES_USE_GENERIC_CONTAINER

  using GenericContainerNamespace::GC_INTEGER ;
  using GenericContainerNamespace::GC_VEC_BOOL ;
  using GenericContainerNamespace::GC_VEC_INTEGER ;
  using GenericContainerNamespace::GC_VEC_REAL ;
  using GenericContainerNamespace::GC_VEC_STRING ;
  using GenericContainerNamespace::GC_VECTOR ;
  using GenericContainerNamespace::GC_MAP ;
  using GenericContainerNamespace::GC_MAT_REAL ;
  using GenericContainerNamespace::mat_real_type ;
  using GenericContainerNamespace::vec_int_type ;
  using GenericContainerNamespace::vec_real_type ;
  using GenericContainerNamespace::vec_string_type ;
  using GenericContainerNamespace::vector_type ;
  using GenericContainerNamespace::map_type ;

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
    GenericContainer const & gc_spline_type = gc("spline_type") ;
    SPLINE_ASSERT( GC_VEC_STRING == gc_spline_type.get_type(),
                   "[SplineSet[" << _name << "]::setup] field `spline_type` expected to be of type `vec_string_type` found: ` " <<
                   gc_spline_type.get_type_name() << "`" ) ;
    vec_string_type const & spline_type_vec = gc_spline_type.get_vec_string() ;
    _nspl = sizeType(spline_type_vec.size()) ;

    SPLINE_ASSERT( gc.exists("xdata") , "[SplineSet[" << _name << "]::setup] missing `xdata` field!") ;
    GenericContainer const & gc_xdata = gc("xdata") ;
    SPLINE_ASSERT( GC_VEC_REAL == gc_xdata.get_type() || GC_VEC_INTEGER == gc_xdata.get_type(),
                   "[SplineSet[" << _name << "]::setup] field `xdata` expected to be of type `vec_real` found: ` " <<
                   gc_xdata.get_type_name() << "`" ) ;

    _npts = sizeType( GC_VEC_REAL == gc_xdata.get_type() ?
                      gc_xdata.get_vec_real().size() : gc_xdata.get_vec_int().size() ) ;

    vector<SplineType> stype(_nspl) ;
    vec_string_type    headers(_nspl) ;

    // allocate memory
    splines.resize(_nspl) ;
    is_monotone.resize(_nspl) ;
    sizeType mem = _npts ;
    for ( sizeType spl = 0 ; spl < _nspl ; ++spl ) {

      bool found = false ;
      sizeType j = 0 ;
      string n = spline_type_vec[spl] ;
      std::transform(n.begin(), n.end(), n.begin(), ::tolower) ;
      while ( !found && Splines::spline_type[j] != nullptr ) {
        found = Splines::spline_type[j] == n ;
        if ( found ) stype[spl] = SplineType(j) ;
        ++j ;
      }
      SPLINE_ASSERT( found, "[SplineSet[" << _name << "]::setup] spline type: ``" << spline_type_vec[spl] << "'' unknown" ) ;

      switch (stype[spl]) {
      case QUINTIC_TYPE:
        mem += _npts ; // Y, Yp, Ypp
      case CUBIC_BASE_TYPE:
      case CUBIC_TYPE:
      case AKIMA_TYPE:
      case BESSEL_TYPE:
      case PCHIP_TYPE:
        mem += _npts ; // Y, Yp
      case CONSTANT_TYPE:
      case LINEAR_TYPE:
        mem += _npts ;
        break;
      case SPLINE_SET_TYPE:
      case BSPLINE_TYPE:
        SPLINE_ASSERT( false,
                       "SplineSet::build\nAt spline n. " << spl <<
                       " named " << headers[spl] <<
                       " unsupported type = " << stype[spl] ) ;
      }
    }

    baseValue.allocate( sizeType(mem + 2*_nspl) ) ;
    basePointer.allocate( sizeType(3*_nspl) ) ;

    _Y    = basePointer(_nspl) ;
    _Yp   = basePointer(_nspl) ;
    _Ypp  = basePointer(_nspl) ;
    _X    = baseValue(_npts) ;
    _Ymin = baseValue(_nspl) ;
    _Ymax = baseValue(_nspl) ;

    if ( GC_VEC_REAL == gc_xdata.get_type() ) {
      vec_real_type const & Xdata = gc_xdata.get_vec_real() ;
      std::copy( Xdata.begin(), Xdata.end(), _X ) ;
    } else {
      vec_int_type const & Xdata = gc_xdata.get_vec_int() ;
      std::copy( Xdata.begin(), Xdata.end(), _X ) ;
    }

    for ( sizeType spl = 0 ; spl < _nspl ; ++spl ) {
      _Yp[spl]  = nullptr ;
      _Ypp[spl] = nullptr ;
      switch (stype[spl]) {
      case QUINTIC_TYPE:
        _Ypp[spl] = baseValue( _npts ) ;
      case CUBIC_BASE_TYPE:
      case CUBIC_TYPE:
      case AKIMA_TYPE:
      case BESSEL_TYPE:
      case PCHIP_TYPE:
        _Yp[spl] = baseValue( _npts ) ;
      case CONSTANT_TYPE:
      case LINEAR_TYPE:
        _Y[spl] = baseValue( _npts ) ;
      case SPLINE_SET_TYPE:
        break;
      case BSPLINE_TYPE: // nothing to do (not implemented)
        break;
      }
    }

    SPLINE_ASSERT( gc.exists("ydata") , "[SplineSet[" << _name << "]::setup] missing `ydata` field!") ;
    GenericContainer const & gc_ydata = gc("ydata") ;

    if ( GC_MAT_REAL == gc_ydata.get_type() || GC_VECTOR == gc_ydata.get_type() ) {
      SPLINE_ASSERT( gc.exists("headers"), "[SplineSet[" << _name << "]::setup] missing `headers` field!") ;
      GenericContainer const & gc_headers = gc("headers") ;
      SPLINE_ASSERT( GC_VEC_STRING == gc_headers.get_type(),
                     "[SplineSet[" << _name << "]::setup] field `headers` expected to be of type `vec_string_type` found: ` " <<
                     gc_headers.get_type_name() << "`" ) ;
      vec_string_type const & headers_in = gc_headers.get_vec_string() ;
      SPLINE_ASSERT( headers_in.size() == _nspl,
                     "[SplineSet[" << _name << "]::setup] field `headers` expected to be of size " << _nspl <<
                      " found of size " << headers_in.size() ) ;
      std::copy( headers_in.begin(), headers_in.end(), headers.begin() ) ;
    }

    if ( GC_MAT_REAL == gc_ydata.get_type() ) {
      mat_real_type const & data = gc_ydata.get_mat_real() ;
      SPLINE_ASSERT( _nspl == data.numCols(),
                     "[SplineSet[" << _name << "]::setup] number of headers [" << _nspl <<
                     "] differs from the number of columns [" << data.numCols() << "] in data" ) ;
      SPLINE_ASSERT( _npts == data.numRows(),
                     "[SplineSet[" << _name << "]::setup] number of points [" << _npts <<
                     "] differs from the numeber of rows [" << data.numCols() << "] in data" ) ;
      for ( sizeType i = 0 ; i < _nspl ; ++i ) std::copy(&data(0,i),&data(_npts,i),_Y[i]) ;
    } else if ( GC_VECTOR == gc_ydata.get_type() ) {
      vector_type const & data = gc_ydata.get_vector() ;
      SPLINE_ASSERT( _nspl == data.size(),
                     "[SplineSet[" << _name << "]::setup] field `data` expected of size " << _nspl <<
                     " found of size " << data.size()  ) ;
      for ( sizeType i = 0 ; i < _nspl ; ++i ) {
        GenericContainer const & datai = data[i] ;
        SPLINE_ASSERT( GC_VEC_REAL == datai.get_type() || GC_VEC_INTEGER == datai.get_type(),
                       "[SplineSet::setup] data[" << i << "] expected to be of type `vec_real_type` found: ` " <<
                        datai.get_type_name() << "`" ) ;
        sizeType nrow = _npts ;
        if ( stype[i] == CONSTANT_TYPE ) --nrow ; // constant spline uses n-1 points
        if ( GC_VEC_REAL == datai.get_type() ) {
          vec_real_type const & yi = datai.get_vec_real() ;
          SPLINE_ASSERT( yi.size() >= nrow,
                         "[SplineSet::setup] data[" << i << "] expected to be of size >= " << nrow << 
                         " found of size " << yi.size() ) ;
          std::copy( yi.begin(), yi.begin() + nrow, _Y[i] ) ;
        } else {
          vec_int_type const & yi = datai.get_vec_int() ;
          SPLINE_ASSERT( yi.size() >= nrow,
                         "[SplineSet::setup] data[" << i << "] expected to be of size >= " << nrow << 
                         " found of size " << yi.size() ) ;
          std::copy( yi.begin(), yi.begin() + nrow, _Y[i] ) ;
        }
      }
    } else if ( GC_MAP == gc_ydata.get_type() ) {
      map_type const & data = gc_ydata.get_map() ;
      SPLINE_ASSERT( data.size() == _nspl,
                    "[SplineSet[" << _name << "]::setup] field `data` expected of size " << _nspl <<
                    " found of size " << data.size() ) ;
      map_type::const_iterator im = data.begin() ;
      for ( sizeType spl = 0 ; im != data.end() ; ++im, ++spl ) {
        headers[spl] = im->first ;
        GenericContainer const & datai = im->second ;
        SPLINE_ASSERT( GC_VEC_REAL    == datai.get_type() ||
                       GC_VEC_INTEGER == datai.get_type() ||
                       GC_VEC_BOOL    == datai.get_type(),
                       "[SplineSet::setup] ydata[" << spl << "] expected to be of type `vec_real_type` found: `" <<
                        datai.get_type_name() << "`" ) ;
        sizeType nrow = _npts ;
        if ( stype[spl] == CONSTANT_TYPE ) --nrow ; // constant spline uses n-1 points
        if ( GC_VEC_REAL == datai.get_type() ) {
          vec_real_type const & yi = datai.get_vec_real() ;
          std::copy( yi.begin(), yi.begin() + nrow, _Y[spl] ) ;
        } else {
          vec_int_type const & yi = datai.get_vec_int() ;
          std::copy( yi.begin(), yi.begin() + nrow, _Y[spl] ) ;
        }
      }
    } else {
      SPLINE_ASSERT( false,
                     "[SplineSet[" << _name << "]::setup] field `data` expected to be of type `mat_real_type`, `vector_type` or `map_type' found: ` " <<
                      gc_ydata.get_type_name() << "`" ) ;
    }

    for ( sizeType spl = 0 ; spl < _nspl ; ++spl ) {

      valueType *& pY   = _Y[spl] ;
      valueType *& pYp  = _Yp[spl] ;
      valueType *& pYpp = _Ypp[spl] ;

      sizeType nrow = _npts ;
      if ( stype[spl] == CONSTANT_TYPE ) --nrow ; // constant spline uses n-1 points
      _Ymin[spl] = *std::min_element( pY,pY+nrow ) ;
      _Ymax[spl] = *std::max_element( pY,pY+nrow ) ;        

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
        is_monotone[spl] = checkCubicSplineMonotonicity( _X, pY, pYp, _npts ) ;
        break;
      case CUBIC_TYPE:
        s = new CubicSpline(h) ;
        static_cast<CubicSpline*>(s)->reserve_external( _npts, _X, pY, pYp ) ;
        static_cast<CubicSpline*>(s)->build( _X, pY, _npts ) ;
        is_monotone[spl] = checkCubicSplineMonotonicity( _X, pY, pYp, _npts ) ;
        break;
      case AKIMA_TYPE:
        s = new AkimaSpline(h) ;
        static_cast<AkimaSpline*>(s)->reserve_external( _npts, _X, pY, pYp ) ;
        static_cast<AkimaSpline*>(s)->build( _X, pY, _npts ) ;
        is_monotone[spl] = checkCubicSplineMonotonicity( _X, pY, pYp, _npts ) ;
        break ;
      case BESSEL_TYPE:
        s = new BesselSpline(h) ;
        static_cast<BesselSpline*>(s)->reserve_external( _npts, _X, pY, pYp ) ;
        static_cast<BesselSpline*>(s)->build( _X, pY, _npts ) ;
        is_monotone[spl] = checkCubicSplineMonotonicity( _X, pY, pYp, _npts ) ;
        break ;
      case PCHIP_TYPE:
        s = new PchipSpline(h) ;
        static_cast<PchipSpline*>(s)->reserve_external( _npts, _X, pY, pYp ) ;
        static_cast<PchipSpline*>(s)->build( _X, pY, _npts ) ;
        is_monotone[spl] = checkCubicSplineMonotonicity( _X, pY, pYp, _npts ) ;
        break ;
      case QUINTIC_TYPE:
        s = new QuinticSpline(h) ;
        static_cast<QuinticSpline*>(s)->reserve_external( _npts, _X, pY, pYp, pYpp ) ;
        static_cast<QuinticSpline*>(s)->build( _X, pY, _npts ) ;
        break;
      case SPLINE_SET_TYPE:
        SPLINE_ASSERT( false,
                       "SplineSet::build\nAt spline n. " << spl << " named " << headers[spl] <<
                       "\nSPLINE_SET_TYPE not allowed as spline type\nin SplineSet::build for " << spl << "-th spline" ) ;
      case BSPLINE_TYPE:
        SPLINE_ASSERT( false,
                       "SplineSet::build\nAt spline n. " << spl << " named " << headers[spl] <<
                       "\nBSPLINE_TYPE not allowed as spline type\nin SplineSet::build for " << spl << "-th spline" ) ;
      }
      header_to_position[s->name()] = spl ;
    }
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
