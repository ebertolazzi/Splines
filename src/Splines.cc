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

/*!
 |
 | \date     October 28, 2015
 | \version  5.2
 | \note     first release Jan 12, 1998
 |
 | \author   Enrico Bertolazzi
 |
 | \par      Affiliation:
 |           Department of Industrial Engineering<br>
 |           University of Trento<br>
 |           Via Sommarive 9, I -- 38123 Trento, Italy <br>
 |           enrico.bertolazzi@unitn.it
 |
\*/

#include "Splines.hh"
#include <cmath>
#include <limits> // std::numeric_limits

#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#endif

#ifdef SPLINES_OS_OSX
  #define UNW_LOCAL_ONLY
  #include <cxxabi.h>
  #include <libunwind.h>
#endif

//! Various kind of splines
namespace Splines {

  using std::abs;

  /*
   backtrace() from:
   https://eli.thegreenplace.net/2015/programmatic-access-to-the-call-stack-in-c/

   to get the line from address
   addr2line 0x400968 -e libunwind_backtrace
  */

  #ifndef SPLINES_OS_OSX
  void backtrace( ostream_type & ) {}
  #else
  void
  backtrace( ostream_type & ost ) {
    unw_cursor_t cursor;
    unw_context_t context;

    // Initialize cursor to current frame for local unwinding.
    unw_getcontext(&context);
    unw_init_local(&cursor, &context);

    // Unwind frames one by one, going up the frame stack.
    while ( unw_step(&cursor) > 0 ) {
      unw_word_t offset, pc;
      unw_get_reg(&cursor, UNW_REG_IP, &pc);
      if ( pc == 0 ) break;
      ost << "0x" << std::hex << pc << ":" << std::dec;
      char sym[256];
      if ( unw_get_proc_name(&cursor, sym, sizeof(sym), &offset) == 0 ) {
        char* nameptr = sym;
        int status;
        char* demangled = abi::__cxa_demangle(sym, nullptr, nullptr, &status);
        if ( status == 0 ) nameptr = demangled;
        ost << " (" << nameptr << "+0x" << std::hex << offset << ")\n" << std::dec;
        std::free(demangled);
      } else {
        ost << " -- error: unable to obtain symbol name for this frame\n";
      }
    }
  }
  #endif

  // cbrt is not available on WINDOWS? or C++ < C++11?
  #ifdef _MSC_VER
    using std::sqrt;
    using std::pow;
    static
    inline
    real_type
    cbrt( real_type x )
    { return pow( x, 1.0/3.0 ); }
  #else
    using std::sqrt;
    using std::pow;
    #if __cplusplus <= 199711L
      static
      inline
      real_type
      cbrt( real_type x )
      { return pow( x, 1.0/3.0 ); }
    #else
      using std::pow;
    #endif
  #endif

  /*\
   |  ____                                _        _          _   _
   | |  _ \ __ _ _ __ __ _ _ __ ___   ___| |_ _ __(_)______ _| |_(_) ___  _ __
   | | |_) / _` | '__/ _` | '_ ` _ \ / _ \ __| '__| |_  / _` | __| |/ _ \| '_ \
   | |  __/ (_| | | | (_| | | | | | |  __/ |_| |  | |/ / (_| | |_| | (_) | | | |
   | |_|   \__,_|_|  \__,_|_| |_| |_|\___|\__|_|  |_/___\__,_|\__|_|\___/|_| |_|
   |
  \*/

  void
  uniform(
    integer         /* dim */,
    integer         npts,
    real_type const []/* pnts    */,
    integer         /* ld_pnts */,
    real_type       t[]
  ) {
    t[0]      = 0;
    t[npts-1] = 1;
    for ( integer k = 1; k < npts-1; ++k )
      t[k] = static_cast<real_type>(k)/static_cast<real_type>(npts);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  chordal(
    integer         dim,
    integer         npts,
    real_type const pnts[],
    integer         ld_pnts,
    real_type       t[]
  ) {
    t[0] = 0;
    real_type const * p0 = pnts;
    for ( integer k = 1; k < npts; ++k ) {
      real_type const * p1 = p0 + ld_pnts;
      real_type dst = 0;
      for ( integer j = 0; j < dim; ++j ) {
        real_type c = p1[j] - p0[j];
        dst += c*c;
      }
      t[k] = t[k-1] + sqrt(dst);
    }
    for ( integer k = 1; k < npts-1; ++k ) t[k] /= t[npts-1];
    t[npts-1] = 1;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  centripetal(
    integer         dim,
    integer         npts,
    real_type const pnts[],
    integer         ld_pnts,
    real_type const alpha,
    real_type       t[]
  ) {
    t[0] = 0;
    real_type const * p0 = pnts;
    for ( integer k = 1; k < npts; ++k ) {
      real_type const * p1 = p0 + ld_pnts;
      real_type dst = 0;
      for ( integer j = 0; j < dim; ++j ) {
        real_type c = p1[j] - p0[j];
        dst += c*c;
      }
      t[k] = t[k-1] + pow(dst,alpha/2);
    }
    for ( integer k = 1; k < npts-1; ++k ) t[k] /= t[npts-1];
    t[npts-1] = 1;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  universal(
    integer         dim,
    integer         npts,
    real_type const pnts[],
    integer         ld_pnts,
    real_type       t[]
  ); // to be done

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  FoleyNielsen(
    integer         dim,
    integer         npts,
    real_type const pnts[],
    integer         ld_pnts,
    real_type       t[]
  ); // to be done

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  FangHung(
    integer         dim,
    integer         npts,
    real_type const pnts[],
    integer         ld_pnts,
    real_type       t[]
  ); // to be done

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  char const * spline_type[] = {
    "constant",    // 0
    "linear",      // 1
    "cubic",       // 2
    "akima",       // 3
    "bessel",      // 4
    "pchip",       // 5
    "quintic",     // 6
    "hermite",     // 7
    "spline set",  // 8
    "spline vec",  // 9
    nullptr
  };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  SplineType
  string_to_splineType( std::string const & nin ) {
    std::string n = nin;
    std::transform(n.begin(), n.end(), n.begin(), ::tolower);
    for ( size_t j = 0; spline_type[j] != nullptr; ++j ) {
      if ( Splines::spline_type[j] == n ) return SplineType(j);
    }
    std::ostringstream ost;
    ost << "string_to_splineType(" << n << ") unknown type\n";
    throw std::runtime_error(ost.str());
  }

  static real_type const machineEps = std::numeric_limits<real_type>::epsilon();
  static real_type const m_2pi      = 6.28318530717958647692528676656; // 2*pi

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

  // num real roots, num complex root
  pair<int,int>
  quadraticRoots(
    real_type const a[3],
    real_type       real[2],
    real_type       imag[2]
  ) {

    // A x^2 + B x + C
    real_type const & C = a[0];
    real_type const & B = a[1];
    real_type const & A = a[2];

    real[0] = real[1] = imag[0] = imag[1] = 0;

    pair<int,int> res(0,0);
    if ( isZero(a[0]) ) {
      real[0] = -B/A;
      res.first = 1; // una singola radice reale
    } else {
      real_type twoA = 2*A;
      real_type d    = B*B - 4*A*C;
      real_type absd = abs(d);
      if ( absd <= 2*machineEps*B*B ) { 
        real[0] = -B/twoA; // EQUAL REAL ROOTS
        res.first = 1; // 2 radici reali coincidenti
      } else {
        real_type r = sqrt(absd);
        if ( d < 0 ) { // COMPLEX ROOTS
          real[0] = real[1] = -B/twoA;
          imag[0] = abs(r/twoA);
          imag[1] = -imag[0];
          res.second = 2; // 2 radici complesse coniugate
        } else {
          // DISTINCT REAL ROOTS
          if ( isZero(B) ) {
            real[0] = abs(r/twoA);
            real[1] = -real[0];
          } else {
            real_type w = -B;
            if ( w > 0 ) w += r; else w -= r;
            w *= 0.5;
            real[0] = C/w;
            real[1] = w/A;
          }
          res.first = 2; // 2 radici reali distinte
        }
      }
    }
    return res;
  }
  
  //! cubic polinomial roots
  /*!
    Compute the roots of the polynomial
    
    \f[ a_0 + a_1 z + a_2 z^2 + a_3 z^3 \f]
    
    and store the results is `real` and `imag`.
    It is assumed that \f$ a_3 \f$ is nonzero.
  */

  pair<int,int>
  cubicRoots(
    real_type const a[4],
    real_type       real[3],
    real_type       imag[3]
  ) {

    // initialize roots
    real[0] = real[1] = real[2] = 
    imag[0] = imag[1] = imag[2] = 0;

    // trivial case
    if ( isZero(a[0]) ) {
      pair<int,int> res = quadraticRoots( a+1, real+1, imag+1 ); // quadratica degenerata
      ++res.first;
      return res;
    }

    // trivial case
    if ( isZero(a[3]) ) return quadraticRoots( a, real, imag ); // cubica degenerata

    // x^3 + A x^2 + B x + C
    real_type const C = a[0]/a[3];
    real_type const B = a[1]/a[3];
    real_type const A = a[2]/a[3];
    
    // p(y-A/3) = y^3 + p*y + q
    real_type const A3 = A/3;
    real_type const p  = B-A*A3;
    real_type const q  = C+A3*(2*(A3*A3)-B);
    
    // scaling equation p(S*z)/S^3 = z^3 + 3*(p/S^2/3)*z + 2*(q/S^3/2)
    real_type const S = max( sqrt(abs(p)), cbrt(abs(q)) );

    // check for a triple root
    if ( S <= machineEps ) {
      real[0] = -A3;
      return pair<int,int>(1,0); // 3 radici reali coincidenti
    }

    real_type const P     = (p/3)/S/S;
    real_type const sqrtP = sqrt(abs(p/3))/S;
    real_type const Q     = (q/2)/S/S/S;

    real_type const d     = P*P*P + Q*Q;
    real_type const sqrtd = sqrt(abs(d));

    pair<int,int> res(0,0);
    if ( sqrtd < abs(q)*machineEps ) {
      // P^3 = - Q^2
      // (x+2*a)(x-a)^2 = x^3 - 3*x*a^2 + 2*a^3
      // cioè -a^2 = P, a^3 = Q ==> a = sqrt(-P)
      real_type tmp = Q > 0 ? sqrtP : -sqrtP;
      real[0] = tmp;
      real[1] = -2*tmp;
      res.first = 2; // 3 radici reali, 2 coincidenti
    } else if ( d > 0 ) {
      // w1 = (- Q + sqrt( P^3 + Q^2 ))^(1/3)
      // w2 = (- Q - sqrt( P^3 + Q^2 ))^(1/3)
      real_type w1, w2;
      if ( Q > 0 ) {
        w2 = - pow( sqrtd + Q, 1.0 / 3.0 );
        w1 = - P / w2;
      } else {
        w1 =   pow( sqrtd - Q, 1.0 / 3.0 );
        w2 = - P / w1;
      }
      real[0] = w1 + w2;
      real[1] =
      real[2] = -0.5*real[0];
      imag[1] = (w1-w2)*sqrt(3.0/4.0);
      imag[2] = -imag[1];
      res.first  = 1;
      res.second = 2; // 1 reale 2 complesse coniugate
    } else { // 3 radici reali
      // w1 = (- Q + I*sqrt(|P^3 + Q^2|) )^(1/3)
      // w2 = (- Q - I*sqrt(|P^3 + Q^2|) )^(1/3)
      real_type angle  = atan2( sqrtd, -Q );
      if ( angle < 0 ) angle += m_2pi;
      angle /= 3;
      real_type re = sqrtP * cos(angle);
      real_type im = sqrtP * sin(angle);
      //if ( Q > 0 ) re = -re;
      real[0]  = 2*re;
      real[1]  = real[2] = -re;
      real[1] += sqrt(3.0) * im;
      real[2] -= sqrt(3.0) * im;
      res.first = 3; // 3 radici reali distinte
    }

    for ( integer i = 0; i < res.first+res.second; ++i ) {
      // scalo radici
      real[i] *= S;
      imag[i] *= S;
      // traslo radici
      real[i] -= A3;
    }

    return res;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  /*!
   | Check if cubic spline with this data is monotone,
   |  -2 non monotone data, -1 no, 0 yes, 1 strictly monotone
  \*/
  integer
  checkCubicSplineMonotonicity(
    real_type const X[],
    real_type const Y[],
    real_type const Yp[],
    integer         npts
  ) {
    // check monotonicity of data: (assuming X monotone)
    integer flag = 1;
    for ( size_t i = 1; i < size_t(npts); ++i ) {
      if ( Y[i-1] > Y[i] ) return -2; // non monotone data
      if ( isZero(Y[i-1]-Y[i]) && X[i-1] < X[i] ) flag = 0; // non strict monotone
    }
    // pag 146 Methods of Shape-Preserving Spline Approximation, K
    for ( size_t i = 1; i < size_t(npts); ++i ) {
      if ( X[i] <= X[i-1] ) continue; // skip duplicate points
      real_type dd = (Y[i]-Y[i-1])/(X[i]-X[i-1]);
      real_type m0 = Yp[i-1]/dd;
      real_type m1 = Yp[i]/dd;
      if ( m0 < 0 || m1 < 0 ) return -1; // non monotone
      if ( m0 <= 3 && m1 <= 3 ) {
        if ( flag > 0 && i > 1              && (isZero(m0) || isZero(m0-3) ) ) flag = 0;
        if ( flag > 0 && i < size_t(npts-1) && (isZero(m1) || isZero(m1-3) ) ) flag = 0;
      } else {
        real_type tmp1 = 2*m0+m1-3;
        real_type tmp2 = 2*(m0+m1-2);
        real_type tmp3 = m0*tmp2-(tmp1*tmp1);
        if ( tmp2 >= 0 ) {
          if ( tmp3 < 0 ) return -1; // non monotone spline
        } else {
          if ( tmp3 > 0 ) return -1;
        }
        if ( isZero(tmp3) ) flag = 0;
      }
    }
    return flag; // passed all check
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  updateInterval(
    integer       & lastInterval,
    real_type       x,
    real_type const X[],
    integer         npts
  ) {

    if ( npts <= 2 ) { lastInterval = 0; return; } // nothing to search

    // find the interval of the support of the B-spline
    real_type const * XL = X+lastInterval;
    if ( XL[1] < x ) { // x on the right
      if ( x >= X[npts-2] ) { // x in [X[npt-2],X[npts-1]]
        lastInterval = npts-2; // last interval
      } else if ( x < XL[2] ) { // x in (XL[1],XL[2])
        ++lastInterval;
      } else { // x >= XL[2] search the right interval
        real_type const * XE = X+npts;
        lastInterval += integer(std::lower_bound( XL, XE, x )-XL);
        real_type const * XX = X+lastInterval;
        if ( x < XX[0] || isZero(XX[0]-XX[1]) ) --lastInterval;
      }
    } else if ( x < XL[0] ) { // on the left
      if ( x <= X[1] ) { // x in [X[0],X[1]]
        lastInterval = 0; // first interval
      } else if ( XL[-1] <= x ) { // x in [XL[-1],XL[0])
        --lastInterval;
      } else {
        lastInterval = integer(std::lower_bound( X, XL, x )-X);
        real_type const * XX = X+lastInterval;
        if ( x < XX[0] || isZero(XX[0]-XX[1]) ) --lastInterval;
      }
    } else {
      // x in the interval [ XL[0], XL[1] ] nothing to do
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Spline::build(
    real_type const x[], integer incx,
    real_type const y[], integer incy,
    integer n
  ) {
    reserve( n );
    for ( integer i = 0; i < n; ++i ) X[i] = x[i*incx];
    for ( integer i = 0; i < n; ++i ) Y[i] = y[i*incy];
    npts = n;
    build();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Spline::info( ostream_type & s ) const {
    s << "Spline `" << _name
      << "` of type: " << type_name()
      << " of order: " << order();
    if ( npts > 0 )
      s << "\nxMin = " << xMin() << " xMax = " << xMax()
        << "\nyMin = " << yMin() << " yMax = " << yMax();
    s << '\n';
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Spline::pushBack( real_type x, real_type y ) {
    if ( npts > 0 ) {
      SPLINE_ASSERT(
        x >= X[size_t(npts-1)], // ammetto punti doppi
        "Spline::pushBack, non monotone insert at insert N. " << npts <<
        "\nX[ " << npts-1 << "] = " << X[size_t(npts-1)] <<
        "\nX[ " << npts   << "] = " << x
      )
    }
    if ( npts_reserved == 0 ) {
      reserve( 2 );
    } else if ( npts >= npts_reserved ) {
      // riallocazione & copia
      integer saved_npts = npts; // salvo npts perche reserve lo azzera
      vector<real_type> Xsaved, Ysaved;
      Xsaved.resize( size_t(npts) );
      Ysaved.resize( size_t(npts) );

      std::copy( X, X+npts, Xsaved.begin() );
      std::copy( Y, Y+npts, Ysaved.begin() );
      reserve( (npts+1) * 2 );
      npts = saved_npts;
      std::copy( Xsaved.begin(), Xsaved.end(), X );
      std::copy( Ysaved.begin(), Ysaved.end(), Y );
    }
    X[npts] = x;
    Y[npts] = y;
    ++npts;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Spline::setOrigin( real_type x0 ) {
    real_type Tx = x0 - X[0];
    real_type *ix = X;
    while ( ix < X+npts ) *ix++ += Tx;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Spline::setRange( real_type xmin, real_type xmax ) {
    SPLINE_ASSERT(
      xmax > xmin,
      "Spline::setRange( " << xmin <<
      " , " << xmax << " ) bad range "
    );
    real_type S  = (xmax - xmin) / ( X[npts-1] - X[0] );
    real_type Tx = xmin - S * X[0];
    for( real_type *ix = X; ix < X+npts; ++ix ) *ix = *ix * S + Tx;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Spline::dump(
    ostream_type & s,
    integer        nintervals,
    char const     header[]
  ) const {
    s << header << '\n';
    real_type dx = (xMax()-xMin())/nintervals;
    for ( integer i = 0; i <= nintervals; ++i ) {
      real_type x = xMin() + i*dx;
      s << x << '\t' << (*this)(x) << '\n';
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  curvature( real_type s, Spline const & X, Spline const & Y ) {
    real_type x_1 = X.D(s);
    real_type x_2 = X.DD(s);
    real_type y_1 = Y.D(s);
    real_type y_2 = Y.DD(s);
    real_type t4  = x_1 * x_1;
    real_type t5  = y_1 * y_1;
    real_type t6  = t4 + t5;
    real_type t7  = sqrt(t6);
    return (x_1 * y_2 - y_1 * x_2) / ( t6 * t7 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  curvature_D( real_type s, Spline const & X, Spline const & Y ) {
    real_type x_1 = X.D(s);
    real_type x_2 = X.DD(s);
    real_type x_3 = X.DDD(s);
    real_type y_1 = Y.D(s);
    real_type y_2 = Y.DD(s);
    real_type y_3 = Y.DDD(s);
    real_type t1  = x_1 * x_1;
    real_type t9  = x_2 * x_2;
    real_type t13 = y_1 * y_1;
    real_type t17 = y_2 * y_2;
    real_type t26 = t1 + t13;
    real_type t27 = t26 * t26;
    real_type t28 = sqrt(t26);
    real_type aa  = y_3 * x_1 - y_1 * x_3;
    real_type bb  = 3 * y_2 * x_2;
    return ( t1 * ( aa - bb ) + t13 * ( aa + bb )
             + 3 * x_1 * y_1 * ( t9 - t17 ) ) / ( t28 * t27 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  curvature_DD( real_type s, Spline const & X, Spline const & Y ) {
    real_type x_1 = X.D(s);
    real_type x_2 = X.DD(s);
    real_type x_3 = X.DDD(s);
    real_type x_4 = X.DDDD(s);
    real_type y_1 = Y.D(s);
    real_type y_2 = Y.DD(s);
    real_type y_3 = Y.DDD(s);
    real_type y_4 = Y.DDDD(s);

    real_type t1  = y_1 * y_1;
    real_type t2  = t1 * t1;
    real_type t12 = x_2 * x_2;
    real_type t13 = t12 * x_2;
    real_type t15 = x_1 * x_3;
    real_type t16 = 9 * t15;
    real_type t17 = y_2 * y_2;
    real_type t21 = x_1 * x_1;
    real_type t22 = x_4 * t21;
    real_type t26 = 9 * x_1 * y_2 * y_3;
    real_type t30 = t21 * x_1;
    real_type t40 = t17 * y_2;
    real_type t64 = t21 + t1;
    real_type t65 = t64 * t64;
    real_type t67 = sqrt(t64);
    return (
      t2 * (x_1 * y_4 + 4 * x_2 * y_3 + 5 * x_3 * y_2 - y_1 * x_4)
      + t1 * y_1 * (3 * t13 + x_2 * (t16 - 12 * t17) - 2 * t22 - t26)
      + t1 * ( x_1 * ( 12 * t40 - 33 * y_2 * t12 )
               + t21 * ( y_2 * x_3 - y_3 * x_2)
               + 2 * y_4 * t30 )
      - y_1 * t21 * ( 12 * t13 - x_2 * (t16+33*t17) + t22 + t26 )
      + t30 * ( y_2 * (12 * t12 - 4 * t15) + y_4 * t21
                -5 * x_1 * x_2 * y_3 - 3 * t40)
    ) / (t67*t65*t64);
  }

  /*
  //    ____  ____   ____                               _
  //   / ___|/ ___| / ___| _   _ _ __  _ __   ___  _ __| |_
  //  | |  _| |     \___ \| | | | '_ \| '_ \ / _ \| '__| __|
  //  | |_| | |___   ___) | |_| | |_) | |_) | (_) | |  | |_
  //   \____|\____| |____/ \__,_| .__/| .__/ \___/|_|   \__|
  //                            |_|   |_|
  */

  #ifndef SPLINES_DO_NOT_USE_GENERIC_CONTAINER

  using GenericContainerNamespace::GC_VEC_REAL;
  using GenericContainerNamespace::vec_real_type;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Spline::setup( GenericContainer const & gc ) {
    /*
    // gc["x"]
    // gc["y"]
    //
    */
    SPLINE_ASSERT(
      gc.exists("x"),
      "Spline[" << _name << "]::setup missing `x` field!"
    );
    SPLINE_ASSERT(
      gc.exists("y"),
      "Spline[" << _name << "]::setup missing `y` field!"
    );

    GenericContainer const & gc_x = gc("x");
    GenericContainer const & gc_y = gc("y");

    vec_real_type x, y;
    {
      std::ostringstream ost;
      ost << "Spline[" << _name << "]::setup, field `x'";
      gc_x.copyto_vec_real ( x, ost.str().c_str() );
    }
    {
      std::ostringstream ost;
      ost << "Spline[" << _name << "]::setup, field `y'";
      gc_y.copyto_vec_real ( y, ost.str().c_str() );
    }
    build( x, y );
  }
  #endif

}
