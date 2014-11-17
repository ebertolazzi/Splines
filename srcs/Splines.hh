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

#ifndef SPLINES_HH
#define SPLINES_HH

#include <iostream>
#include <vector>
#include <algorithm>

//
// file: Splines
//

/*!
 \mainpage  Splines
 \author    Enrico Bertolazzi (enrico.bertolazzi@unitn.it), homepage: http://www.ing.unitn.it/~bertolaz
 \version   1.0.0
 \date      2013
 \copyright GNU Public License.

 \details

 This library is available at

 - https://github.com/ebertolazzi/Splines
 - https://bitbucket.org/ebertolazzi/splines

 Splines
 =======

 `Splines` is a set of C++ classes which implements varios spline interpolation.
 The classes are the following:

 - ConstantsSpline, for piecewise constants functions
 - LinearSpline, for piecewise linear interpolation
 - CubicSpline, for classical cubic spline interpolation
 - AkimaSpline, for Akima "non oscillatory" spline interpolation
 - BesselSpline, for Bessel "non oscillatory" spline interpolation
 - PchipSpline,
 - QuinticSpline, Simple quintic spline based oin PCHIP with 4th
 derivative set to 0 at nodal points

 **References**

 - F.N. Fritsch and R.E. Carlson,
 Monotone Piecewise Cubic Interpolation,
 SIAM Journal of Numerical Analysis, Vol. 17, No. 2, pp. 238-246,
 April 1980.

 **Usage**

 The usage is simple:

 ~~~~~~~~~~~~~
 #include "Splines.hh"
 using namespace SplinesLoad ;

 ....

 CubicSpline spline ;
 double x[] = {1,2,3,4} ;
 double y[] = {3,1,1,3} ;
 spline.build(x,y,4) ; // build a cubic spline with 4 points

 cout << spline(1.1) << '\n';     // spline at x = 1.1
 cout << spline.D(1.1) << '\n';   // spline first derivative at x = 1.1
 cout << spline.DD(1.1) << '\n';  // spline second derivative at x = 1.1
 cout << spline.DDD(1.1) << '\n'; // spline third derivative at x = 1.1
 ~~~~~~~~~~~~~

 splines can be built incrementally

 ~~~~~~~~~~~~~
 #include "Splines.hh"
 using namespace SplinesLoad ;

 ....

 CubicSpline spline ;

 spline.pushBack( 1, 3 ) ;
 spline.pushBack( 2, 1 ) ;
 spline.pushBack( 3, 1 ) ;
 spline.pushBack( 4, 3 ) ;
 spline.build() ;

 cout << spline(1.1) << '\n';     // spline at x = 1.1
 cout << spline.D(1.1) << '\n';   // spline first derivative at x = 1.1
 cout << spline.DD(1.1) << '\n';  // spline second derivative at x = 1.1
 cout << spline.DDD(1.1) << '\n'; // spline third derivative at x = 1.1
 ~~~~~~~~~~~~~

 or by using standard vector

 ~~~~~~~~~~~~~
 #include "Splines.hh"
 #include <vector>
 using namespace SplinesLoad ;
 using namespace std ;

 ....

 CubicSpline spline ;
 std::vector x, y ;
 x.push_back(1) ; y.push_back(3) ;
 x.push_back(2) ; y.push_back(1) ;
 x.push_back(3) ; y.push_back(1) ;
 x.push_back(4) ; y.push_back(3) ;
 spline.build(x,y) ;

 cout << spline(1.1) << '\n';     // spline at x = 1.1
 cout << spline.D(1.1) << '\n';   // spline first derivative at x = 1.1
 cout << spline.DD(1.1) << '\n';  // spline second derivative at x = 1.1
 cout << spline.DDD(1.1) << '\n'; // spline third derivative at x = 1.1
 ~~~~~~~~~~~~~

 **Compile and tests**

 Edit makefile file to match compiler of your OS and do:

 make

 To run the test

 make run

 To generate documentation (using DOXYGEN: http://www.stack.nl/~dimitri/doxygen/index.html)

 make doc

 **DOXYGEN documentation**

 Available at: http://www.ing.unitn.it/~bertolaz/4-software/Splines/index.html

*/

// if C++ < C++11 define nullptr
#if __cplusplus <= 199711L
  #ifndef nullptr
    #include <cstddef>
    #define nullptr NULL
  #endif
#endif

#ifndef SPLINE_ASSERT
  #include <stdexcept>
  #include <sstream>
  #define SPLINE_ASSERT(COND,MSG)           \
    if ( !(COND) ) {                        \
      std::ostringstream ost ;              \
      ost << "On line: " << __LINE__        \
          << " file: " << __FILE__          \
          << '\n' << MSG << '\n' ;          \
      throw std::runtime_error(ost.str()) ; \
    }
#endif

#ifndef SPLINE_WARNING
  #include <stdexcept>
  #include <sstream>
  #define SPLINE_WARNING(COND,MSG)         \
    if ( !(COND) ) {                       \
      std::cout << "On line: " << __LINE__ \
                << " file: " << __FILE__   \
                << MSG << '\n' ;           \
    }
#endif

#ifdef DEBUG
  #define SPLINE_CHECK_NAN( PTR, MSG, DIM ) Splines::checkNaN( PTR, MSG, DIM )
#else
  #define SPLINE_CHECK_NAN( PTR, MSG, DIM )
#endif

//! Various kind of splines
namespace Splines {

  using namespace ::std ; // load standard namspace

  typedef double            valueType ;
  typedef unsigned          sizeType ;
  typedef int               indexType ;
  typedef vector<valueType> VectorOfValues ;

  /*       _               _    _   _       _   _
  //   ___| |__   ___  ___| | _| \ | | __ _| \ | |
  //  / __| '_ \ / _ \/ __| |/ /  \| |/ _` |  \| |
  // | (__| | | |  __/ (__|   <| |\  | (_| | |\  |
  //  \___|_| |_|\___|\___|_|\_\_| \_|\__,_|_| \_|
  */
  void
  checkNaN( valueType const pv[],
            char      const v_name[],
            sizeType  const DIM ) ;

  /*
  //   _   _                     _ _       
  //  | | | | ___ _ __ _ __ ___ (_) |_ ___ 
  //  | |_| |/ _ \ '__| '_ ` _ \| | __/ _ \
  //  |  _  |  __/ |  | | | | | | | ||  __/
  //  |_| |_|\___|_|  |_| |_| |_|_|\__\___|
  */
  void Hermite3    ( valueType const x, valueType const H, valueType base[4]     ) ;
  void Hermite3_D  ( valueType const x, valueType const H, valueType base_D[4]   ) ;
  void Hermite3_DD ( valueType const x, valueType const H, valueType base_DD[4]  ) ;
  void Hermite3_DDD( valueType const x, valueType const H, valueType base_DDD[4] ) ;

  void Hermite5      ( valueType const x, valueType const H, valueType base[6]       ) ;
  void Hermite5_D    ( valueType const x, valueType const H, valueType base_D[6]     ) ;
  void Hermite5_DD   ( valueType const x, valueType const H, valueType base_DD[6]    ) ;
  void Hermite5_DDD  ( valueType const x, valueType const H, valueType base_DDD[6]   ) ;
  void Hermite5_DDDD ( valueType const x, valueType const H, valueType base_DDDD[6]  ) ;
  void Hermite5_DDDDD( valueType const x, valueType const H, valueType base_DDDDD[6] ) ;

  /*
  //   ____  _ _ _
  //  | __ )(_) (_)_ __   ___  __ _ _ __
  //  |  _ \| | | | '_ \ / _ \/ _` | '__|
  //  | |_) | | | | | | |  __/ (_| | |
  //  |____/|_|_|_|_| |_|\___|\__,_|_|
  */
  valueType
  bilinear( valueType const p[4],
            valueType const M[4][4],
            valueType const q[4] ) ;

  /*
  //   ____        _ _            
  //  / ___| _ __ | (_)_ __   ___ 
  //  \___ \| '_ \| | | '_ \ / _ \
  //   ___) | |_) | | | | | |  __/
  //  |____/| .__/|_|_|_| |_|\___|
  //        |_|                   
  */
  //! Spline Management Class
  /*!
   * 
   * \date     October 11, 2011
   * \version  1.0
   * \note     first release October 11, 2011
   *
   * \author   Enrico Bertolazzi
   *
   * \par      Affiliation:
   *           Department of Mechanics and Structures Engineering <br>
   *           University of Trento <br>
   *           via Mesiano 77, I -- 38050 Trento, Italy <br>
   *           enrico.bertolazzi@ing.unitn.it
   *
   */
  class Spline {
  protected:

    sizeType       npts ;
    VectorOfValues X, Y ;

    sizeType search( valueType x ) const ;
    mutable sizeType lastInterval ;

    void
    allocate( valueType const x[], sizeType incx,
              valueType const y[], sizeType incy,
              sizeType n ) ;

  public:

    //! spline constructor
    Spline()
    : npts(0)
    , X()
    , Y()
    , lastInterval(0)
    {}

    //! spline destructor
    virtual 
    ~Spline()
    {}

    //! return the number of support points of the spline.
    sizeType numPoints(void) const { return npts ; }

    //! return the i-th node of the spline (x component).
    valueType xNode( sizeType i ) const { return X[i] ; }

    //! return the i-th node of the spline (y component).
    valueType yNode( sizeType i ) const { return Y[i] ; }

    //! Add a support point (x,y) to the spline.
    void pushBack( valueType x, valueType y ) ;

    //! Drop a support point to the spline.
    void dropBack() ;

    virtual
    void
    build (void) = 0 ;

    virtual
    void
    build ( valueType const x[], valueType const y[], sizeType n ) = 0 ;

    virtual
    void
    build ( valueType const x[], sizeType incx,
            valueType const y[], sizeType incy,
            sizeType n ) = 0 ;

    virtual
    void
    build ( VectorOfValues const & x, VectorOfValues const & y ) = 0 ;

    //! Cancel the support points, empty the spline.
    virtual
    void
    clear (void) {
      X.clear() ;
      Y.clear() ;
      lastInterval = 0 ;
      npts         = 0 ;
    }

    //! return x-minumum spline value
    valueType xMin() const { return X.front() ; }

    //! return x-maximum spline value
    valueType xMax() const { return X.back()  ; }

    ///////////////////////////////////////////////////////////////////////////
    //! change X-origin of the spline
    void setOrigin( valueType x0 ) ;

    //! change X-range of the spline
    void setRange( valueType xmin, valueType xmax ) ;

    ///////////////////////////////////////////////////////////////////////////
    //! Evalute spline value
    virtual valueType operator () ( valueType x ) const = 0 ;

    //! First derivative
    virtual valueType D( valueType x ) const = 0 ;

    //! Second derivative
    virtual valueType DD( valueType x ) const = 0 ;

    //! Third derivative
    virtual valueType DDD( valueType x ) const = 0 ;

    //! Print spline coefficients
    virtual void writeToStream( std::basic_ostream<char> & s ) const = 0 ;

    //! Return spline typename
    virtual char const * type_name() const = 0 ;

  } ;

  /*
  //    ____      _     _        ____        _ _              ____                 
  //   / ___|   _| |__ (_) ___  / ___| _ __ | (_)_ __   ___  | __ )  __ _ ___  ___ 
  //  | |  | | | | '_ \| |/ __| \___ \| '_ \| | | '_ \ / _ \ |  _ \ / _` / __|/ _ \
  //  | |__| |_| | |_) | | (__   ___) | |_) | | | | | |  __/ | |_) | (_| \__ \  __/
  //   \____\__,_|_.__/|_|\___| |____/| .__/|_|_|_| |_|\___| |____/ \__,_|___/\___|
  //                                  |_|
  */
  //! cubic spline base class
  class CubicSplineBase : public Spline {
  protected:
    VectorOfValues    Yp ;
    mutable valueType base[4] ;
    mutable valueType base_D[4] ;
    mutable valueType base_DD[4] ;
    mutable valueType base_DDD[4] ;

  public:
  
    //! spline constructor
    CubicSplineBase()
    : Spline()
    , Yp()
    {}
    
    virtual
    ~CubicSplineBase()
    {}

    void copySpline( CubicSplineBase const & S ) ;
    void
    allocate( valueType const x[], sizeType incx,
              valueType const y[], sizeType incy,
              sizeType n ) ;

    //! return the i-th node of the spline (y' component).
    valueType ypNode( sizeType i ) const { return Yp[i] ; }

    //! change X-range of the spline
    void setRange( valueType xmin, valueType xmax ) ;

    ///////////////////////////////////////////////////////////////////////////
    //! Evalute spline value
    virtual valueType operator () ( valueType x ) const ;

    //! First derivative
    virtual valueType D( valueType x ) const ;

    //! Second derivative
    virtual valueType DD( valueType x ) const ;

    //! Third derivative
    virtual valueType DDD( valueType x ) const ;

    //! Print spline coefficients
    virtual void writeToStream( std::basic_ostream<char> & s ) const ;

  } ;

  /*
  //      _    _    _                   ____        _ _            
  //     / \  | | _(_)_ __ ___   __ _  / ___| _ __ | (_)_ __   ___ 
  //    / _ \ | |/ / | '_ ` _ \ / _` | \___ \| '_ \| | | '_ \ / _ \
  //   / ___ \|   <| | | | | | | (_| |  ___) | |_) | | | | | |  __/
  //  /_/   \_\_|\_\_|_| |_| |_|\__,_| |____/| .__/|_|_|_| |_|\___|
  //                                         |_|                   
  */
  //! Akima spline class
  /*!
   *  Reference
   *  =========
   *  Hiroshi Akima, Journal of the ACM, Vol. 17, No. 4, October 1970, pages 589-602.
   */
  class AkimaSpline : public CubicSplineBase {
  private:
    valueType akima_one( valueType epsi,
                         valueType di_m2,
                         valueType di_m1,
                         valueType di,
                         valueType di_p1 ) const ;
  public:

    //! spline constructor
    AkimaSpline()
    : CubicSplineBase()
    {}

    //! spline destructor
    virtual
    ~AkimaSpline()
    {}

    //! Build an Akima spline from previously inserted points
    virtual
    void
    build (void) ;

    //! Build an Akima spline.
    /*!
     * \param x    vector of x-coordinates
     * \param incx access elements as x[0], x[incx], x[2*incx],...
     * \param y    vector of y-coordinates
     * \param incy access elements as y[0], y[incy], x[2*incy],...
     * \param n    total number of points
     */
    virtual
    void
    build ( valueType const x[], sizeType incx,
            valueType const y[], sizeType incy,
            sizeType n )
    { allocate( x, incx, y, incy, n ) ; build() ; }

    //! Build an Akima spline.
    /*!
     * \param x vector of x-coordinates
     * \param y vector of y-coordinates
     * \param n total number of points
     */
    virtual
    void
    build ( valueType const x[], valueType const y[], sizeType n )
    { allocate( x, 1, y, 1, n ) ; build() ; }

    //! Build an Akima spline.
    /*!
     * \param x vector of x-coordinates
     * \param y vector of y-coordinates
     */
    virtual
    void
    build ( VectorOfValues const & x, VectorOfValues const & y )
    { build( &x.front(), &y.front(), sizeType(x.size()) ) ; }

    //! Return spline typename
    virtual char const * type_name() const { return "akima" ; }

  } ;

  /*
  //   ____                     _ ____        _ _            
  //  | __ )  ___  ___ ___  ___| / ___| _ __ | (_)_ __   ___ 
  //  |  _ \ / _ \/ __/ __|/ _ \ \___ \| '_ \| | | '_ \ / _ \
  //  | |_) |  __/\__ \__ \  __/ |___) | |_) | | | | | |  __/
  //  |____/ \___||___/___/\___|_|____/| .__/|_|_|_| |_|\___|
  //                                   |_|                   
  */
  //! Bessel spline class
  class BesselSpline : public CubicSplineBase {
  public:

    //! spline constructor
    BesselSpline()
    : CubicSplineBase()
    {}

    //! spline destructor
    virtual
    ~BesselSpline()
    {}

    //! Build a Bessel spline from previously inserted points
    virtual
    void
    build (void) ;

    //! Build a Bessel spline.
    /*!
     * \param x    vector of x-coordinates
     * \param incx access elements as x[0], x[incx], x[2*incx],...
     * \param y    vector of y-coordinates
     * \param incy access elements as y[0], y[incy], x[2*incy],...
     * \param n    total number of points
     */
    virtual
    void
    build ( valueType const x[], sizeType incx,
            valueType const y[], sizeType incy,
            sizeType n )
    { allocate( x, incx, y, incy, n ) ; build() ; }

    //! Build a Bessel spline.
    /*!
     * \param x vector of x-coordinates
     * \param y vector of y-coordinates
     * \param n total number of points
     */
    virtual
    void
    build ( valueType const x[], valueType const y[], sizeType n )
    { allocate( x, 1, y, 1, n ) ; build() ; }

    //! Build a Bessel spline.
    /*!
     * \param x vector of x-coordinates
     * \param y vector of y-coordinates
     */
    virtual
    void
    build ( VectorOfValues const & x, VectorOfValues const & y )
    { build( &x.front(), &y.front(), sizeType(x.size()) ) ; }

    //! Return spline typename
    virtual char const * type_name() const { return "bessel" ; }

  } ;

  /*
  //    ____      _     _      ____        _ _            
  //   / ___|   _| |__ (_) ___/ ___| _ __ | (_)_ __   ___ 
  //  | |  | | | | '_ \| |/ __\___ \| '_ \| | | '_ \ / _ \
  //  | |__| |_| | |_) | | (__ ___) | |_) | | | | | |  __/
  //   \____\__,_|_.__/|_|\___|____/| .__/|_|_|_| |_|\___|
  //                                |_|                   
  */
  //! Cubic Spline Management Class
  /*!
   * 
   * \date     Feruary 25, 2008
   * \version  5.0
   * \note     first release Jan 12, 1998 
   *
   * \author   Enrico Bertolazzi
   *
   * \par      Affiliation:
   *           Department of Industrial Engineering <br>
   *           University of Trento <br>
   *           via Mesiano 77, I -- 38050 Trento, Italy <br>
   *           enrico.bertolazzi@ing.unitn.it
   *
   */
  class CubicSpline : public CubicSplineBase {
  private:
    valueType ddy0 ;
    valueType ddyn ;
  public:

    //! spline constructor
    CubicSpline()
    : CubicSplineBase()
    , ddy0(0)
    , ddyn(0)
    {}

    //! spline destructor
    virtual
    ~CubicSpline()
    {}

    /*!
     * \param ddy0  first boundary condition.
     *              The second derivative at initial point.
     * \param ddyn  second boundary condition.
     *              The second derivative at final point.
     */
    void
    setbc( valueType ddy0, valueType ddyn ) {
      this -> ddy0 = ddy0 ;
      this -> ddyn = ddyn ;    
    }

    //! Build a spline from previously inserted points
    virtual
    void
    build (void) ;

    //! Build a cubic spline.
    /*!
     * \param x    vector of x-coordinates
     * \param incx access elements as x[0], x[incx], x[2*incx],...
     * \param y    vector of y-coordinates
     * \param incy access elements as y[0], y[incy], x[2*incy],...
     * \param n    total number of points
     */
    virtual
    void
    build ( valueType const x[], sizeType incx,
            valueType const y[], sizeType incy,
            sizeType n )
    { allocate( x, incx, y, incy, n ) ; build() ; }

    //! Build a cubic spline.
    /*!
     * \param x vector of x-coordinates
     * \param y vector of y-coordinates
     * \param n total number of points
     */
    virtual
    void
    build ( valueType const x[], valueType const y[], sizeType n )
    { allocate( x, 1, y, 1, n ) ; build() ; }

    //! Build a cubic spline.
    /*!
     * \param x vector of x-coordinates
     * \param y vector of y-coordinates
     */
    virtual
    void
    build ( VectorOfValues const & x, VectorOfValues const & y )
    { build( &x.front(), &y.front(), sizeType(x.size()) ) ; }

    //! Return spline typename
    virtual char const * type_name() const { return "cubic" ; }

  } ;

  /*
  //   ____      _     _      ____        _ _            
  //  |  _ \ ___| |__ (_)_ __/ ___| _ __ | (_)_ __   ___ 
  //  | |_) / __| '_ \| | '_ \___ \| '_ \| | | '_ \ / _ \
  //  |  __/ (__| | | | | |_) |__) | |_) | | | | | |  __/
  //  |_|   \___|_| |_|_| .__/____/| .__/|_|_|_| |_|\___|
  //                    |_|        |_|                   
  */
  //! Pchip (Piecewise Cubic Hermite Interpolating Polynomial) spline class
  class PchipSpline : public CubicSplineBase {
  public:

    //! spline constructor
    PchipSpline()
    : CubicSplineBase()
    {}

    //! spline destructor
    virtual
    ~PchipSpline()
    {}

    //! Build a Monotone spline from previously inserted points
    virtual
    void
    build (void) ;

    //! Build a Monotone spline.
    /*!
     * \param x    vector of x-coordinates
     * \param incx access elements as x[0], x[incx], x[2*incx],...
     * \param y    vector of y-coordinates
     * \param incy access elements as y[0], y[incy], x[2*incy],...
     * \param n    total number of points
     */
    virtual
    void
    build ( valueType const x[], sizeType incx,
            valueType const y[], sizeType incy,
            sizeType n )
    { allocate( x, incx, y, incy, n ) ; build() ; }

    //! Build a Monotone spline.
    /*!
     * \param x vector of x-coordinates
     * \param y vector of y-coordinates
     * \param n total number of points
     */
    virtual
    void
    build ( valueType const x[], valueType const y[], sizeType n )
    { allocate( x, 1, y, 1, n ) ; build() ; }

    //! Build a Monotone spline.
    /*!
     * \param x vector of x-coordinates
     * \param y vector of y-coordinates
     */
    virtual
    void
    build ( VectorOfValues const & x, VectorOfValues const & y )
    { build( &x.front(), &y.front(), sizeType(x.size()) ) ; }

    //! Return spline typename
    virtual char const * type_name() const { return "pchip" ; }

  } ;

  /*
  //   _     _                       ____        _ _            
  //  | |   (_)_ __   ___  __ _ _ __/ ___| _ __ | (_)_ __   ___ 
  //  | |   | | '_ \ / _ \/ _` | '__\___ \| '_ \| | | '_ \ / _ \
  //  | |___| | | | |  __/ (_| | |   ___) | |_) | | | | | |  __/
  //  |_____|_|_| |_|\___|\__,_|_|  |____/| .__/|_|_|_| |_|\___|
  //                                      |_|                   
  */
  //! Linear spline class
  class LinearSpline : public Spline {
  public:

    LinearSpline()
    : Spline()
    {}

    virtual
    ~LinearSpline()
    {}

    //! given x and y vectors build a linear spline
    /*!
     * \param x    vector of x-coordinates
     * \param incx access elements as x[0], x[incx], x[2*incx],...
     * \param y    vector of y-coordinates
     * \param incy access elements as y[0], y[incy], x[2*incy],...
     * \param n    total number of points
     */
    virtual
    void
    build ( valueType const x[], sizeType incx,
            valueType const y[], sizeType incy,
            sizeType n ) {
      SPLINE_ASSERT( n > 1,
                     "LinearSpline::build, n =" << n <<
                     " not enought point to define a spline\n") ;
      npts = n ;
      X.resize(n) ;
      Y.resize(n) ;
      for ( sizeType i = 0 ; i < n ; ++i )
        { X[i] = x[i*incx] ; Y[i] = y[i*incy] ; }
    }

    //! given x and y vectors build a linear spline
    virtual
    void
    build( valueType const x[], valueType const y[], sizeType n ) {
      SPLINE_ASSERT( n > 1,
                     "LinearSpline::build, n =" << n <<
                     " not enought point to define a spline\n") ;
      npts = n ;
      X.resize(n) ;
      Y.resize(n) ;
      std::copy( x, x+n, X.begin() ) ;
      std::copy( y, y+n, Y.begin() ) ;
    }

    //! given x and y vectors build a linear spline
    virtual
    void
    build ( VectorOfValues const & x, VectorOfValues const & y )
    { build( &x.front(), &y.front(), sizeType(x.size()) ) ; }

    //! added for compatibility with cubic splines
    virtual
    void
    build()
    {}

    //! Evalute spline value at `x`
    virtual
    valueType
    operator () ( valueType x ) const {
      if ( x < X.front() ) return Y.front() ;
      if ( x > X.back()  ) return Y.back() ;
      sizeType i = search(x) ;
      valueType s = (x-X[i])/(X[i+1] - X[i]) ;
      return (1-s)*Y[i] + s * Y[i+1] ;
    }

    //! First derivative
    virtual
    valueType
    D( valueType x ) const {
      if ( x < X.front() ) return 0 ;
      if ( x > X.back()  ) return 0 ;
      sizeType i = search(x) ;
      return ( Y[i+1] - Y[i] ) / ( X[i+1] - X[i] ) ;
    }
    
    //! Second derivative
    virtual valueType DD( valueType ) const { return 0 ; }

    //! Third derivative
    virtual valueType DDD( valueType ) const { return 0 ; }

    //! Print spline coefficients
    virtual void writeToStream( std::basic_ostream<char> & s ) const ;

    //! Return spline typename
    virtual char const * type_name() const { return "linear" ; }

  } ;

  /*  
  //    ____                _              _       ____        _ _            
  //   / ___|___  _ __  ___| |_ __ _ _ __ | |_ ___/ ___| _ __ | (_)_ __   ___ 
  //  | |   / _ \| '_ \/ __| __/ _` | '_ \| __/ __\___ \| '_ \| | | '_ \ / _ \
  //  | |__| (_) | | | \__ \ || (_| | | | | |_\__ \___) | |_) | | | | | |  __/
  //   \____\___/|_| |_|___/\__\__,_|_| |_|\__|___/____/| .__/|_|_|_| |_|\___|
  //                                                    |_|                   
  */
  //! Picewise constants spline class
  class ConstantsSpline : public Spline {
  private:
  public:

    ConstantsSpline()
    : Spline()
    {}

    ~ConstantsSpline()
    {}

    //! Cancel the support points, empty the spline. `x0` is the initial abscissa of empty spline
    virtual
    void
    clear( valueType x0 ) {
      X.clear() ; X.push_back( x0 ) ;
      Y.clear() ;
      npts = 1 ;
    }

    /*!
     * \param x    vector of x-coordinates
     * \param incx access elements as x[0], x[incx], x[2*incx],...
     * \param y    vector of y-coordinates
     * \param incy access elements as y[0], y[incy], x[2*incy],...
     * \param n    total number of points
     */
    virtual
    void
    build ( valueType const x[], sizeType incx,
            valueType const y[], sizeType incy,
            sizeType n ) {
      SPLINE_ASSERT( n > 1,
                     "ConstantsSpline::build, n =" << n <<
                     " not enought point to define a spline\n") ;
      npts = n ;
      X.resize(n) ;
      Y.resize(n-1) ;
      for ( sizeType i = 0 ; i < n   ; ++i ) X[i] = x[i*incx] ;
      for ( sizeType i = 0 ; i < n-1 ; ++i ) Y[i] = y[i*incy] ;
    }

    //! given x and y vectors build a piecewise constants spline
    virtual
    void
    build( valueType const x[], valueType const y[], sizeType n ) {
      SPLINE_ASSERT( n > 1,
                     "ConstantsSpline::build, n =" << n <<
                     " not enought point to define a spline\n") ;
      npts = n ;
      X.resize(n) ;
      Y.resize(n-1) ;
      std::copy( x, x+n,   X.begin() ) ;
      std::copy( y, y+n-1, Y.begin() ) ;
    }

    //! Build a piecewise constants spline
    /*!
     * \param x vector of x-coordinates
     * \param y vector of y-coordinates
     */
    virtual
    void
    build ( VectorOfValues const & x, VectorOfValues const & y )
    { build( &x.front(), &y.front(), sizeType(x.size()) ) ; }

    //! added for compatibility with cubic splines
    virtual
    void
    build()
    {}

    //! Evalute spline value at `x`
    virtual
    valueType
    operator () ( valueType x ) const {
      if ( x < X.front() ) return Y.front() ;
      if ( x > X.back()  ) return Y.back() ;
      return Y[search(x)] ;
    }

    //! First derivative
    virtual valueType D( valueType ) const { return 0 ; }
    
    //! Second derivative
    virtual valueType DD( valueType ) const { return 0 ; }

    //! Third derivative
    virtual valueType DDD( valueType ) const { return 0 ; }

    //! Print spline coefficients
    virtual void writeToStream( std::basic_ostream<char> & ) const ;

    //! Return spline typename
    virtual char const * type_name() const { return "constant" ; }

  } ;

  /*
  //    ___        _       _   _      ____        _ _            ____                 
  //   / _ \ _   _(_)_ __ | |_(_) ___/ ___| _ __ | (_)_ __   ___| __ )  __ _ ___  ___ 
  //  | | | | | | | | '_ \| __| |/ __\___ \| '_ \| | | '_ \ / _ \  _ \ / _` / __|/ _ \
  //  | |_| | |_| | | | | | |_| | (__ ___) | |_) | | | | | |  __/ |_) | (_| \__ \  __/
  //   \__\_\\__,_|_|_| |_|\__|_|\___|____/| .__/|_|_|_| |_|\___|____/ \__,_|___/\___|
  //                                       |_|                                        
  //  
  */
  //! cubic quintic base class
  class QuinticSplineBase : public Spline {
  protected:
    VectorOfValues    Yp, Ypp ;
    mutable valueType base[6] ;
    mutable valueType base_D[6] ;
    mutable valueType base_DD[6] ;
    mutable valueType base_DDD[6] ;
    mutable valueType base_DDDD[6] ;
    mutable valueType base_DDDDD[6] ;

  public:
  
    //! spline constructor
    QuinticSplineBase()
    : Spline()
    , Yp()
    {}
    
    virtual
    ~QuinticSplineBase()
    {}

    void copySpline( QuinticSplineBase const & S ) ;
    void allocate( valueType const x[], sizeType incx,
                   valueType const y[], sizeType incy,
                   sizeType n ) ;

    //! return the i-th node of the spline (y' component).
    valueType ypNode( sizeType i ) const { return Yp[i] ; }

    //! change X-range of the spline
    void setRange( valueType xmin, valueType xmax ) ;

    ///////////////////////////////////////////////////////////////////////////
    //! Evalute spline value
    virtual valueType operator () ( valueType x ) const ;

    //! First derivative
    virtual valueType D( valueType x ) const ;

    //! Second derivative
    virtual valueType DD( valueType x ) const ;

    //! Third derivative
    virtual valueType DDD( valueType x ) const ;

    //! Fourth derivative
    virtual valueType DDDD( valueType x ) const ;

    //! Fifth derivative
    virtual valueType DDDDD( valueType x ) const ;

    //! Print spline coefficients
    virtual void writeToStream( std::basic_ostream<char> & s ) const ;

    //! Return spline typename
    virtual char const * type_name() const { return "quintic" ; }
  } ;
 
  /*
  //    ___        _       _   _      ____        _ _            
  //   / _ \ _   _(_)_ __ | |_(_) ___/ ___| _ __ | (_)_ __   ___ 
  //  | | | | | | | | '_ \| __| |/ __\___ \| '_ \| | | '_ \ / _ \
  //  | |_| | |_| | | | | | |_| | (__ ___) | |_) | | | | | |  __/
  //   \__\_\\__,_|_|_| |_|\__|_|\___|____/| .__/|_|_|_| |_|\___|
  //                                       |_|                   
  //  
  */
  //! Quintic spline class
  class QuinticSpline : public QuinticSplineBase {
  public:

    //! spline constructor
    QuinticSpline()
    : QuinticSplineBase()
    {}

    //! spline destructor
    virtual
    ~QuinticSpline()
    {}

    //! Build a Monotone quintic spline from previously inserted points
    void build (void) ;

    //! Build a Monotone spline.
    /*!
     * \param x    vector of x-coordinates
     * \param incx access elements as x[0], x[incx], x[2*incx],...
     * \param y    vector of y-coordinates
     * \param incy access elements as y[0], y[incy], x[2*incy],...
     * \param n    total number of points
     */
    virtual
    void
    build ( valueType const x[], sizeType incx,
            valueType const y[], sizeType incy,
            sizeType n )
    { allocate( x, incx, y, incy, n ) ; build() ; }

    //! Build a Monotone spline.
    /*!
     * \param x vector of x-coordinates
     * \param y vector of y-coordinates
     * \param n total number of points
     */
    void
    build ( valueType const x[], valueType const y[], sizeType n )
    { allocate( x, 1, y, 1, n ) ; build() ; }

    //! Build a Monotone spline.
    /*!
     * \param x vector of x-coordinates
     * \param y vector of y-coordinates
     */
    void
    build ( VectorOfValues const & x, VectorOfValues const & y )
    { build( &x.front(), &y.front(), sizeType(x.size()) ) ; }
  } ;

  /*
  //   ____        _ _            ____              __
  //  / ___| _ __ | (_)_ __   ___/ ___| _   _ _ __ / _|
  //  \___ \| '_ \| | | '_ \ / _ \___ \| | | | '__| |_
  //   ___) | |_) | | | | | |  __/___) | |_| | |  |  _|
  //  |____/| .__/|_|_|_| |_|\___|____/ \__,_|_|  |_|
  //        |_|
  */
  //! Spline Management Class
  class SplineSurf {
  protected:

    VectorOfValues X, Y, Z ;

    mutable sizeType lastInterval_x ;
    sizeType search_x( valueType x ) const ;

    mutable sizeType lastInterval_y ;
    sizeType search_y( valueType y ) const ;

    sizeType ipos( sizeType i, sizeType j ) const { return i+sizeType(X.size())*j ; }

  public:

    //! spline constructor
    SplineSurf()
    : X()
    , Y()
    , Z()
    , lastInterval_x(0)
    , lastInterval_y(0)
    {}

    //! spline destructor
    virtual 
    ~SplineSurf()
    {}

    //! Cancel the support points, empty the spline.
    virtual
    void
    clear (void) {
      X.clear() ;
      Y.clear() ;
      Z.clear() ;
      lastInterval_x = 0 ;
      lastInterval_y = 0 ;
    }

    //! return the number of support points of the spline along x direction
    sizeType numPointX(void) const { return sizeType(X.size()) ; }

    //! return the number of support points of the spline along y direction
    sizeType numPointY(void) const { return sizeType(Y.size()) ; }

    //! return the i-th node of the spline (x component).
    valueType xNode( sizeType i ) const { return X[i] ; }

    //! return the i-th node of the spline (y component).
    valueType yNode( sizeType i ) const { return Y[i] ; }

    //! return the i-th node of the spline (y component).
    valueType zNode( sizeType i, sizeType j ) const { return Z[ipos(i,j)] ; }

    //! return x-minumum spline value
    valueType xMin() const { return X.front() ; }

    //! return x-maximum spline value
    valueType xMax() const { return X.back()  ; }

    //! return x-minumum spline value
    valueType yMin() const { return Y.front() ; }

    //! return x-maximum spline value
    valueType yMax() const { return Y.back()  ; }

    ///////////////////////////////////////////////////////////////////////////

    virtual
    void
    build ( valueType const x[], sizeType incx,
            valueType const y[], sizeType incy,
            valueType const z[], sizeType incz,
            sizeType nx, sizeType ny ) = 0 ;

    virtual
    void
    build ( VectorOfValues const & x,
            VectorOfValues const & y,
            VectorOfValues const & z ) = 0 ;

    virtual
    void
    build ( valueType const z[], sizeType incz, sizeType nx, sizeType ny ) = 0 ;

    virtual
    void
    build ( VectorOfValues const & z, sizeType nx, sizeType ny ) = 0 ;

    //! Evalute spline value
    virtual valueType operator () ( valueType x, valueType y ) const = 0 ;

    //! First derivative
    virtual void D( valueType x, valueType y, valueType d[3] ) const = 0 ;
    virtual valueType Dx( valueType x, valueType y ) const = 0 ;
    virtual valueType Dy( valueType x, valueType y ) const = 0 ;

    //! Second derivative
    virtual void DD( valueType x, valueType y, valueType dd[6] ) const = 0 ;
    virtual valueType Dxx( valueType x, valueType y ) const = 0 ;
    virtual valueType Dxy( valueType x, valueType y ) const = 0 ;
    virtual valueType Dyy( valueType x, valueType y ) const = 0 ;

    //! Print spline coefficients
    virtual void writeToStream( ostream & s ) const = 0 ;

    //! Return spline typename
    virtual char const * type_name() const = 0 ;

  } ;

  /*
  //   ____  _ _ _                       ____        _ _
  //  | __ )(_) (_)_ __   ___  __ _ _ __/ ___| _ __ | (_)_ __   ___
  //  |  _ \| | | | '_ \ / _ \/ _` | '__\___ \| '_ \| | | '_ \ / _ \
  //  | |_) | | | | | | |  __/ (_| | |   ___) | |_) | | | | | |  __/
  //  |____/|_|_|_|_| |_|\___|\__,_|_|  |____/| .__/|_|_|_| |_|\___|
  //                                          |_|
  */
  //! bilinear spline base class
  class BilinearSpline : public SplineSurf {
  public:
  
    //! spline constructor
    BilinearSpline()
    : SplineSurf()
    {}
    
    virtual
    ~BilinearSpline()
    {}

    ///////////////////////////////////////////////////////////////////////////

    virtual
    void
    build ( valueType const x[], sizeType incx,
            valueType const y[], sizeType incy,
            valueType const z[], sizeType incz,
            sizeType nx, sizeType ny ) ;

    virtual
    void
    build ( VectorOfValues const & x,
            VectorOfValues const & y,
            VectorOfValues const & z ) ;

    virtual
    void
    build ( valueType const z[],  sizeType incz, sizeType nx, sizeType ny ) ;

    virtual
    void
    build ( VectorOfValues const & z, sizeType nx, sizeType ny ) ;

    //! Evalute spline value
    virtual valueType operator () ( valueType x, valueType y ) const ;

    //! First derivative
    virtual void D( valueType x, valueType y, valueType d[3] ) const ;
    virtual valueType Dx( valueType x, valueType y ) const ;
    virtual valueType Dy( valueType x, valueType y ) const ;

    //! Second derivative
    virtual void DD( valueType x, valueType y, valueType dd[6] ) const { D(x,y,dd) ; dd[3] = dd[4] = dd[5] = 0 ; }
    virtual valueType Dxx( valueType , valueType ) const { return 0 ; }
    virtual valueType Dxy( valueType , valueType ) const { return 0 ; }
    virtual valueType Dyy( valueType , valueType ) const { return 0 ; }

    //! Print spline coefficients
    virtual void writeToStream( ostream & s ) const ;

    //! Return spline typename
    virtual char const * type_name() const ;

  } ;

  /*
  //   ____  _  ____      _     _      ____        _ _            ____
  //  | __ )(_)/ ___|   _| |__ (_) ___/ ___| _ __ | (_)_ __   ___| __ )  __ _ ___  ___
  //  |  _ \| | |  | | | | '_ \| |/ __\___ \| '_ \| | | '_ \ / _ \  _ \ / _` / __|/ _ \
  //  | |_) | | |__| |_| | |_) | | (__ ___) | |_) | | | | | |  __/ |_) | (_| \__ \  __/
  //  |____/|_|\____\__,_|_.__/|_|\___|____/| .__/|_|_|_| |_|\___|____/ \__,_|___/\___|
  //                                        |_|
  */
  //! cubic spline base class
  class BiCubicSplineBase : public SplineSurf {
  protected:

    VectorOfValues    DX, DY, DXY ;
    mutable valueType u[4] ;
    mutable valueType u_D[4] ;
    mutable valueType u_DD[4] ;
    mutable valueType v[4] ;
    mutable valueType v_D[4] ;
    mutable valueType v_DD[4] ;
    mutable valueType bili[4][4] ;

    void load( sizeType i, sizeType j ) const ;
    virtual void makeSpline() = 0 ;

  public:
  
    //! spline constructor
    BiCubicSplineBase()
    : SplineSurf()
    , DX()
    , DY()
    {}
    
    virtual
    ~BiCubicSplineBase()
    {}

    valueType DxNode ( sizeType i, sizeType j ) const { return DX[ipos(i,j)] ; }
    valueType DyNode ( sizeType i, sizeType j ) const { return DY[ipos(i,j)] ; }
    valueType DxyNode( sizeType i, sizeType j ) const { return DXY[ipos(i,j)] ; }

    //! Evalute spline value
    virtual valueType operator () ( valueType x, valueType y ) const ;

    //! First derivative
    virtual void D( valueType x, valueType y, valueType d[3] ) const ;
    virtual valueType Dx( valueType x, valueType y ) const ;
    virtual valueType Dy( valueType x, valueType y ) const ;

    //! Second derivative
    virtual void DD( valueType x, valueType y, valueType dd[6] ) const ;
    virtual valueType Dxx( valueType x, valueType y ) const ;
    virtual valueType Dxy( valueType x, valueType y ) const ;
    virtual valueType Dyy( valueType x, valueType y ) const ;

    ///////////////////////////////////////////////////////////////////////////

    virtual
    void
    build ( valueType const x[], sizeType incx,
            valueType const y[], sizeType incy,
            valueType const z[], sizeType incz,
            sizeType nx, sizeType ny ) ;

    virtual
    void
    build ( VectorOfValues const & x,
            VectorOfValues const & y,
            VectorOfValues const & z ) ;

    virtual
    void
    build ( valueType const z[],  sizeType incz, sizeType nx, sizeType ny ) ;

    virtual
    void
    build ( VectorOfValues const & z, sizeType nx, sizeType ny ) ;

  } ;

  /*
  //   ____  _  ____      _     _      ____        _ _            
  //  | __ )(_)/ ___|   _| |__ (_) ___/ ___| _ __ | (_)_ __   ___
  //  |  _ \| | |  | | | | '_ \| |/ __\___ \| '_ \| | | '_ \ / _ \
  //  | |_) | | |__| |_| | |_) | | (__ ___) | |_) | | | | | |  __/
  //  |____/|_|\____\__,_|_.__/|_|\___|____/| .__/|_|_|_| |_|\___|
  //                                        |_|
  */
  //! cubic spline base class
  class BiCubicSpline : public BiCubicSplineBase {

    virtual void makeSpline() ;

  public:
  
    //! spline constructor
    BiCubicSpline()
    : BiCubicSplineBase()
    {}
    
    virtual
    ~BiCubicSpline()
    {}

    //! Print spline coefficients
    virtual void writeToStream( ostream & s ) const ;

    //! Return spline typename
    virtual char const * type_name() const ;

  } ;

  /*
  //      _    _    _                 ____  ____            _ _
  //     / \  | | _(_)_ __ ___   __ _|___ \|  _ \ ___ _ __ | (_)_ __   ___
  //    / _ \ | |/ / | '_ ` _ \ / _` | __) | | | / __| '_ \| | | '_ \ / _ \
  //   / ___ \|   <| | | | | | | (_| |/ __/| |_| \__ \ |_) | | | | | |  __/
  //  /_/   \_\_|\_\_|_| |_| |_|\__,_|_____|____/|___/ .__/|_|_|_| |_|\___|
  //                                                 |_|
  */
  //! cubic spline base class
  class Akima2Dspline : public BiCubicSplineBase {

    virtual void makeSpline() ;

  public:
  
    //! spline constructor
    Akima2Dspline()
    : BiCubicSplineBase()
    {}
    
    virtual
    ~Akima2Dspline()
    {}

    //! Print spline coefficients
    virtual void writeToStream( ostream & s ) const ;

    //! Return spline typename
    virtual char const * type_name() const ;

  } ;

}

namespace SplinesLoad {
  using Splines::Spline ;
  using Splines::AkimaSpline ;
  using Splines::BesselSpline ;
  using Splines::PchipSpline ;
  using Splines::CubicSpline ;
  using Splines::LinearSpline ;
  using Splines::ConstantsSpline ;
  using Splines::QuinticSpline ;

  using Splines::BilinearSpline ;
  using Splines::BiCubicSpline ;
  using Splines::Akima2Dspline ;
}

#endif
