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
#include "SplinesUtils.hh"

#include <limits>
#include <cmath>

#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wimplicit-fallthrough"
#endif

/**
 * 
 */

namespace Splines {

  using std::abs;
  using std::sqrt;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer
  SplineSet::Treap::search( std::string const & id ) const {
    // binary search
    size_t U = data.size();
    size_t L = 0;
    while ( U-L > 1 ) {
      size_t pos = (L+U)>>1;
      std::string const & id_pos = data[pos].first;
      if ( id_pos < id ) L = pos; else U = pos;
    }
    if ( data[L].first == id ) return integer(L);
    if ( data[U].first == id ) return integer(U);
    return -1; // non trovato
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineSet::Treap::insert( std::string const & id, integer position ) {
    size_t pos = data.size();
    data.push_back(DATA_TYPE(id,position));
    while ( pos > 0 ) {
      size_t pos1 = pos-1;
      data[pos].first  = data[pos1].first;
      data[pos].second = data[pos1].second;
      if ( data[pos1].first < id ) break;
      pos = pos1;
    }
    data[pos] = DATA_TYPE(id,position);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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
  {
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  //! spline destructor
  SplineSet::~SplineSet() {
    baseValue.free();
    basePointer.free();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineSet::info( ostream_type & s ) const {
    s << "SplineSet[" << name() << "] n.points = "
      << _npts << " n.splines = " << _nspl << '\n';
    for ( size_t i = 0; i < size_t(_nspl); ++i ) {
      s << "\nSpline n." << i;
      switch ( is_monotone[i] ) {
        case -2: s << " with NON monotone data\n"; break;
        case -1: s << " is NOT monotone\n";        break;
        case  0: s << " is monotone\n";            break;
        case  1: s << " is strictly monotone\n";   break;
        default: SPLINE_DO_ERROR(
          "SplineSet::info classification: " << is_monotone[i] <<
          " not in range {-2,-1,0,1}"
        )
      }
      splines[i]->info(s);
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineSet::dump_table( ostream_type & stream, integer num_points ) const {
    vector<real_type> vals;
    stream << 's';
    for ( integer i = 0; i < numSplines(); ++i )
      stream << '\t' << header(i);
    stream << '\n';

    for ( integer j = 0; j < num_points; ++j ) {
      real_type s = xMin() + ((xMax()-xMin())*j)/(num_points-1);
      this->eval( s, vals );
      stream << s;
      for ( integer i = 0; i < numSplines(); ++i )
        stream << '\t' << vals[size_t(i)];
      stream << '\n';
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer
  SplineSet::getPosition( char const * hdr ) const {
    integer pos = header_to_position.search(hdr);
    SPLINE_ASSERT(
      pos >= 0,
      "SplineSet::getPosition(\"" << hdr << "\") not found!\n" <<
      "available keys: " << name_list()
    )
    return pos;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineSet::build(
    integer            nspl,
    integer            npts,
    char         const *headers[],
    SplineType1D const stype[],
    real_type    const X[],
    real_type    const *Y[],
    real_type    const *Yp[]
  ) {
    SPLINE_ASSERT(
      nspl > 0,
      "SplineSet::build expected positive nspl = " << nspl
    )
    SPLINE_ASSERT(
      npts > 1,
      "SplineSet::build expected npts = " << npts << " greather than 1"
    )
    this->_nspl = nspl;
    this->_npts = npts;
    // allocate memory
    splines.resize(size_t(this->_nspl));
    is_monotone.resize(size_t(this->_nspl));
    integer mem = npts;
    for ( integer spl = 0; spl < nspl; ++spl ) {
      switch (stype[size_t(spl)]) {
      case QUINTIC_TYPE:
        mem += npts; // Y, Yp, Ypp
      case CUBIC_TYPE:
      case AKIMA_TYPE:
      case BESSEL_TYPE:
      case PCHIP_TYPE:
      case HERMITE_TYPE:
        mem += npts; // Y, Yp
      case CONSTANT_TYPE:
      case LINEAR_TYPE:
        mem += npts;
      break;
      case SPLINE_SET_TYPE:
      case SPLINE_VEC_TYPE:
      //default:
        SPLINE_ASSERT(
          false,
          "SplineSet::build\nAt spline n. " << spl <<
          " named " << headers[spl] <<
          " cannot be done for type = " << stype[spl]
        )
      }
    }

    this->baseValue   . allocate( size_t(mem + 2*nspl) );
    this->basePointer . allocate( size_t(3*nspl) );

    this->_Y    = this->basePointer(size_t(this->_nspl));
    this->_Yp   = this->basePointer(size_t(this->_nspl));
    this->_Ypp  = this->basePointer(size_t(this->_nspl));
    this->_X    = this->baseValue(size_t(this->_npts));
    this->_Ymin = this->baseValue(size_t(this->_nspl));
    this->_Ymax = this->baseValue(size_t(this->_nspl));

    std::copy( X, X+npts, this->_X );
    for ( size_t spl = 0; spl < size_t(nspl); ++spl ) {
      real_type *& pY   = this->_Y[spl];
      real_type *& pYp  = this->_Yp[spl];
      real_type *& pYpp = this->_Ypp[spl];
      pY = baseValue(size_t(this->_npts));
      std::copy( Y[spl], Y[spl]+npts, pY );
      if ( stype[spl] == CONSTANT_TYPE ) {
        this->_Ymin[spl] = *std::min_element( pY, pY+npts-1 );
        this->_Ymax[spl] = *std::max_element( pY, pY+npts-1 );
      } else {
        this->_Ymin[spl] = *std::min_element( pY, pY+npts );
        this->_Ymax[spl] = *std::max_element( pY, pY+npts );
      }
      pYpp = pYp = nullptr;
      switch ( stype[size_t(spl)] ) {
      case QUINTIC_TYPE:
        pYpp = baseValue(size_t(this->_npts));
      case CUBIC_TYPE:
      case AKIMA_TYPE:
      case BESSEL_TYPE:
      case PCHIP_TYPE:
      case HERMITE_TYPE:
        pYp = baseValue(size_t(this->_npts));
        if ( stype[spl] == HERMITE_TYPE ) {
          SPLINE_ASSERT(
            Yp != nullptr && Yp[spl] != nullptr,
            "SplineSet::build\nAt spline n. " << spl <<
            " named " << headers[spl] << "\nexpect to find derivative values"
          )
          std::copy( Yp[spl], Yp[spl]+npts, pYp );
        }
      case CONSTANT_TYPE:
      case LINEAR_TYPE:
      case SPLINE_SET_TYPE:
      case SPLINE_VEC_TYPE:
      //default:
        break;
      }
      string h = headers[spl];
      Spline * & s = splines[spl];

      is_monotone[spl] = -1;
      switch (stype[size_t(spl)]) {
      case CONSTANT_TYPE:
        s = new ConstantSpline(h);
        static_cast<ConstantSpline*>(s)->reserve_external( this->_npts, this->_X, pY );
        static_cast<ConstantSpline*>(s)->npts = _npts;
        static_cast<ConstantSpline*>(s)->build();
        break;

      case LINEAR_TYPE:
        s = new LinearSpline(h);
        static_cast<LinearSpline*>(s)->reserve_external( this->_npts, this->_X, pY );
        static_cast<LinearSpline*>(s)->npts = _npts;
        static_cast<LinearSpline*>(s)->build();
        // check monotonicity of data
        { integer flag = 1;
          for ( integer j = 1; j < this->_npts; ++j ) {
            if ( pY[j-1] > pY[j] ) { flag = -1; break; } // non monotone data
            if ( isZero(pY[j-1]-pY[j]) && this->_X[j-1] < this->_X[j] ) flag = 0; // non strict monotone
          }
          is_monotone[spl] = flag;
        }
        break;

      case CUBIC_TYPE:
        s = new CubicSpline(h);
        static_cast<CubicSpline*>(s)->reserve_external( this->_npts, _X, pY, pYp );
        static_cast<CubicSpline*>(s)->npts = this->_npts;
        static_cast<CubicSpline*>(s)->build();
        is_monotone[spl] = checkCubicSplineMonotonicity( _X, pY, pYp, _npts );
        break;

      case AKIMA_TYPE:
        s = new AkimaSpline(h);
        static_cast<AkimaSpline*>(s)->reserve_external( this->_npts, _X, pY, pYp );
        static_cast<AkimaSpline*>(s)->npts = this->_npts;
        static_cast<AkimaSpline*>(s)->build();
        is_monotone[spl] = checkCubicSplineMonotonicity( this->_X, pY, pYp, this->_npts );
        break;

      case BESSEL_TYPE:
        s = new BesselSpline(h);
        static_cast<BesselSpline*>(s)->reserve_external( this->_npts, this->_X, pY, pYp );
        static_cast<BesselSpline*>(s)->npts = this->_npts;
        static_cast<BesselSpline*>(s)->build();
        is_monotone[spl] = checkCubicSplineMonotonicity( this->_X, pY, pYp, this->_npts );
        break;

      case PCHIP_TYPE:
        s = new PchipSpline(h);
        static_cast<PchipSpline*>(s)->reserve_external( this->_npts, this->_X, pY, pYp );
        static_cast<PchipSpline*>(s)->npts = this->_npts;
        static_cast<PchipSpline*>(s)->build();
        is_monotone[spl] = checkCubicSplineMonotonicity( this->_X, pY, pYp, this->_npts );
        break;

      case HERMITE_TYPE:
        s = new HermiteSpline(h);
        static_cast<CubicSpline*>(s)->reserve_external( this->_npts, this->_X, pY, pYp );
        static_cast<CubicSpline*>(s)->npts = this->_npts;
        static_cast<CubicSpline*>(s)->build();
        is_monotone[spl] = checkCubicSplineMonotonicity( this->_X, pY, pYp, this->_npts );
        break;

      case QUINTIC_TYPE:
        s = new QuinticSpline(h);
        static_cast<QuinticSpline*>(s)->reserve_external( this->_npts, this->_X, pY, pYp, pYpp );
        static_cast<QuinticSpline*>(s)->npts = this->_npts;
        static_cast<QuinticSpline*>(s)->build();
        break;

      case SPLINE_SET_TYPE:
      case SPLINE_VEC_TYPE:
      //default:
        SPLINE_ASSERT(
          false,
          "SplineSet::build\nAt spline n. " << spl << " named " << headers[spl] <<
          "\n" << stype[size_t(spl)] <<
          " not allowed as spline type\nin SplineSet::build for " << spl <<
          "-th spline"
        )
      }
      this->header_to_position.insert( s->name(), integer(spl) );
    }

    this->baseValue   . must_be_empty( "SplineSet::build, baseValue" );
    this->basePointer . must_be_empty( "SplineSet::build, basePointer" );

  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineSet::getHeaders( vector<string> & h ) const {
    h.resize(size_t(this->_nspl));
    for ( size_t i = 0; i < size_t(this->_nspl); ++i )
      h[i] = this->splines[i]->name();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineSet::eval( real_type x, vector<real_type> & vals ) const {
    vals.resize(size_t(this->_nspl));
    for ( size_t i = 0; i < size_t(this->_nspl); ++i )
      vals[i] = (*this->splines[i])(x);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineSet::eval( real_type x, real_type vals[], integer incy ) const {
    size_t ii = 0;
    for ( size_t i = 0; i < size_t(this->_nspl); ++i, ii += size_t(incy) )
      vals[ii] = (*this->splines[i])(x);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineSet::eval_D( real_type x, vector<real_type> & vals ) const {
    vals.resize(size_t(this->_nspl));
    for ( size_t i = 0; i < size_t(this->_nspl); ++i )
      vals[i] = this->splines[i]->D(x);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineSet::eval_D( real_type x, real_type vals[], integer incy ) const {
    size_t ii = 0;
    for ( size_t i = 0; i < size_t(this->_nspl); ++i, ii += size_t(incy) )
      vals[ii] = this->splines[i]->D(x);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineSet::eval_DD( real_type x, vector<real_type> & vals ) const {
    vals.resize(size_t(this->_nspl));
    for ( size_t i = 0; i < size_t(this->_nspl); ++i )
      vals[i] = this->splines[i]->DD(x);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineSet::eval_DD( real_type x, real_type vals[], integer incy ) const {
    size_t ii = 0;
    for ( size_t i = 0; i < size_t(this->_nspl); ++i, ii += size_t(incy) )
      vals[ii] = this->splines[i]->DD(x);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineSet::eval_DDD( real_type x, vector<real_type> & vals ) const {
    vals.resize(size_t(_nspl));
    for ( size_t i = 0; i < size_t(this->_nspl); ++i )
      vals[i] = this->splines[i]->DDD(x);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineSet::eval_DDD( real_type x, real_type vals[], integer incy ) const {
    size_t ii = 0;
    for ( size_t i = 0; i < size_t(this->_nspl); ++i, ii += size_t(incy) )
      vals[ii] = this->splines[i]->DDD(x);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // vectorial values

  Spline const *
  SplineSet::intersect(
    integer     spl,
    real_type   zeta,
    real_type & x
  ) const {
    SPLINE_ASSERT(
      spl >= 0 && spl < _nspl,
      "Spline n." << spl << " is not in SplineSet"
    )
    SPLINE_ASSERT(
      this->is_monotone[size_t(spl)]>0,
      "Spline n." << spl << " is not monotone and can't be used as independent"
    )
    Spline const * S = this->splines[size_t(spl)];
    // cerco intervallo intersezione
    real_type const * X = this->_Y[size_t(spl)];
    SPLINE_ASSERT(
      zeta >= X[0] && zeta <= X[size_t(this->_npts-1)],
      "SplineSet, evaluation at zeta = " << zeta <<
      " is out of range: [" << X[0] << ", " << X[size_t(_npts-1)] << "]"
    )

    integer interval = integer(lower_bound( X, X+this->_npts, zeta ) - X);
    if ( interval > 0 ) --interval;
    if ( isZero(X[size_t(interval)]-X[size_t(interval+1)]) ) ++interval; // degenerate interval for duplicated nodes
    if ( interval >= this->_npts-1 ) interval = this->_npts-2;

    // compute intersection
    real_type a  = this->_X[size_t(interval)];
    real_type b  = this->_X[size_t(interval+1)];
    real_type ya = X[size_t(interval)];
    real_type yb = X[size_t(interval+1)];
    real_type DX = b-a;
    real_type DY = yb-ya;
    SPLINE_ASSERT(
      zeta >= ya && zeta <= yb,
      "SplineSet, Bad interval [ " << ya << "," << yb << "] for zeta = " << zeta
    )
    SPLINE_ASSERT(
      a < b,
      "SplineSet, Bad x interval [ " << a << "," << b << "]"
    )
    if ( S->type() == LINEAR_TYPE ) {
      x = a + (b-a)*(zeta-ya)/(yb-ya);
    } else {
      real_type const * dX = this->_Yp[size_t(spl)];
      real_type        dya = dX[interval];
      real_type        dyb = dX[interval+1];
      real_type coeffs[4] = { ya-zeta, dya, (3*DY/DX-2*dya-dyb)/DX, (dyb+dya-2*DY/DX)/(DX*DX) };
      real_type real[3], imag[3];
      pair<int,int> icase = cubicRoots( coeffs, real, imag );
      SPLINE_ASSERT(
        icase.first > 0,
        "SplineSet, No intersection found with independent spline at zeta = " << zeta
      )
      // cerca radice buona
      bool ok = false;
      for ( integer i = 0; i < icase.first && !ok; ++i ) {
        ok = real[i] >= 0 && real[i] <= DX;
        if ( ok ) x = a + real[i];
      }
      SPLINE_ASSERT(
        ok,
        "SplineSet, failed to find intersection with independent spline at zeta = " << zeta
      )
    }
    return S;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineSet::eval2(
    integer             spl,
    real_type           zeta,
    vector<real_type> & vals
  ) const {
    real_type x;
    intersect( spl, zeta, x );
    vals.resize(size_t(this->_nspl));
    for ( size_t i = 0; i < size_t(this->_nspl); ++i )
      vals[i] = (*this->splines[i])(x);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineSet::eval2(
    integer   spl,
    real_type zeta,
    real_type vals[],
    integer   incy
  ) const {
    real_type x;
    intersect( spl, zeta, x );
    size_t ii = 0;
    for ( size_t i = 0; i < size_t(this->_nspl); ++i, ii += size_t(incy) )
      vals[ii] = (*this->splines[i])(x);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineSet::eval2_D(
    integer             spl,
    real_type           zeta,
    vector<real_type> & vals
  ) const {
    real_type x;
    Spline const * S = intersect( spl, zeta, x );
    real_type ds = S->D(x);
    vals.resize(size_t(this->_nspl));
    for ( size_t i = 0; i < size_t(this->_nspl); ++i )
      vals[i] = this->splines[i]->D(x)/ds;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineSet::eval2_D(
    integer   spl,
    real_type zeta,
    real_type vals[],
    integer   incy
  ) const {
    real_type x;
    Spline const * S = intersect( spl, zeta, x );
    real_type ds = S->D(x);
    size_t ii = 0;
    for ( size_t i = 0; i < size_t(this->_nspl); ++i, ii += size_t(incy) )
      vals[ii] = this->splines[i]->D(x)/ds;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineSet::eval2_DD(
    integer             spl,
    real_type           zeta,
    vector<real_type> & vals
  ) const {
    real_type x;
    Spline const * S = intersect( spl, zeta, x );
    real_type dt  = 1/S->D(x);
    real_type dt2 = dt*dt;
    real_type ddt = -S->DD(x)*(dt*dt2);
    vals.resize(size_t(this->_nspl));
    for ( size_t i = 0; i < size_t(this->_nspl); ++i )
      vals[i] = this->splines[i]->DD(x)*dt2 + this->splines[i]->D(x)*ddt;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineSet::eval2_DD(
    integer   spl,
    real_type zeta,
    real_type vals[],
    integer   incy
  ) const {
    real_type x;
    Spline const * S = intersect( spl, zeta, x );
    real_type dt  = 1/S->D(x);
    real_type dt2 = dt*dt;
    real_type ddt = -S->DD(x)*(dt*dt2);
    size_t ii = 0;
    for ( size_t i = 0; i < size_t(this->_nspl); ++i, ii += size_t(incy) )
      vals[ii] = this->splines[i]->DD(x)*dt2 + this->splines[i]->D(x)*ddt;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineSet::eval2_DDD(
    integer             spl,
    real_type           zeta,
    vector<real_type> & vals
  ) const {
    real_type x;
    Spline const * S = intersect( spl, zeta, x );
    real_type dt  = 1/S->D(x);
    real_type dt3 = dt*dt*dt;
    real_type ddt = -S->DD(x)*dt3;
    real_type dddt = 3*(ddt*ddt)/dt-S->DDD(x)*(dt*dt3);
    vals.resize( size_t(this->_nspl) );
    for ( size_t i = 0; i < size_t(this->_nspl); ++i )
      vals[i] = this->splines[i]->DDD(x)*dt3 +
                3*this->splines[i]->DD(x)*dt*ddt +
                this->splines[i]->D(x)*dddt;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineSet::eval2_DDD(
    integer   spl,
    real_type zeta,
    real_type vals[],
    integer   incy
  ) const {
    real_type x;
    Spline const * S = intersect( spl, zeta, x );
    real_type dt  = 1/S->D(x);
    real_type dt3 = dt*dt*dt;
    real_type ddt = -S->DD(x)*dt3;
    real_type dddt = 3*(ddt*ddt)/dt-S->DDD(x)*(dt*dt3);
    size_t ii = 0;
    for ( size_t i = 0; i < size_t(this->_nspl); ++i, ii += size_t(incy) )
      vals[ii] = this->splines[i]->DDD(x)*dt3 +
                 3*this->splines[i]->DD(x)*dt*ddt +
                 this->splines[i]->D(x)*dddt;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  SplineSet::eval2(
    real_type    zeta,
    char const * indep,
    char const * name
  ) const {
    vector<real_type> vals;
    this->eval2( this->getPosition(indep), zeta, vals );
    return vals[size_t(this->getPosition(name))];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  SplineSet::eval2_D(
    real_type    zeta,
    char const * indep,
    char const * name
  ) const {
    vector<real_type> vals;
    this->eval2_D( this->getPosition(indep), zeta, vals );
    return vals[size_t(this->getPosition(name))];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  SplineSet::eval2_DD(
    real_type    zeta,
    char const * indep,
    char const * name
  ) const {
    vector<real_type> vals;
    this->eval2_DD( this->getPosition(indep), zeta, vals );
    return vals[size_t(this->getPosition(name))];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  SplineSet::eval2_DDD(
    real_type    zeta,
    char const * indep,
    char const * name
  ) const {
    vector<real_type> vals;
    this->eval2_DDD( this->getPosition(indep), zeta, vals );
    return vals[size_t(this->getPosition(name))];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  SplineSet::eval2( real_type zeta, integer indep, integer spl ) const {
    vector<real_type> vals;
    this->eval2( indep, zeta, vals );
    return vals[size_t(spl)];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  SplineSet::eval2_D( real_type zeta, integer indep, integer spl ) const {
    vector<real_type> vals;
    this->eval2_D( indep, zeta, vals );
    return vals[size_t(spl)];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  SplineSet::eval2_DD( real_type zeta, integer indep, integer spl ) const {
    vector<real_type> vals;
    this->eval2_DD( indep, zeta, vals );
    return vals[size_t(spl)];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  SplineSet::eval2_DDD( real_type zeta, integer indep, integer spl ) const {
    vector<real_type> vals;
    this->eval2_DDD( indep, zeta, vals );
    return vals[size_t(spl)];
  }
}
