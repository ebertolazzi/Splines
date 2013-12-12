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

 spline . pushBack( 1, 3 ) ;
 spline . pushBack( 2, 1 ) ;
 spline . pushBack( 3, 1 ) ;
 spline . pushBack( 4, 3 ) ;
 spline . build() ;

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
 spline . build(x,y) ;

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
          << MSG << '\n' ;                  \
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
    indexType      lastInterval ;

    sizeType
    search( valueType x ) const {
      VectorOfValues::difference_type i ;
      i = lower_bound( X.begin(), X.end(), x ) - X.begin() ;
      if ( i > 0 ) --i ;
      return sizeType(i) ;
    }

    void allocate( valueType const x[], valueType const y[], sizeType n ) ;

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
    build ( VectorOfValues const & x, VectorOfValues const & y ) = 0 ;

    //! Cancel the support points, empty the spline.
    virtual
    void
    clear (void) {
      X . clear() ;
      Y . clear() ;
      lastInterval = 0 ;
      npts         = 0 ;
    }

    //! return x-minumum spline value
    valueType xMin() const { return X . front() ; }

    //! return x-maximum spline value
    valueType xMax() const { return X . back()  ; }

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
    virtual void writeToStream( ostream & s ) const = 0 ;

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
    void allocate( valueType const x[], valueType const y[], sizeType n ) ;

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
    virtual void writeToStream( ostream & s ) const ;

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
     * \param x vector of x-coordinates
     * \param y vector of y-coordinates
     * \param n total number of points
     */
    virtual
    void
    build ( valueType const x[], valueType const y[], sizeType n )
    { allocate( x, y, n ) ; build() ; }

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
     * \param x vector of x-coordinates
     * \param y vector of y-coordinates
     * \param n total number of points
     */
    virtual
    void
    build ( valueType const x[], valueType const y[], sizeType n )
    { allocate( x, y, n ) ; build() ; }

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
     * \param x vector of x-coordinates
     * \param y vector of y-coordinates
     * \param n total number of points
     */
    virtual
    void
    build ( valueType const x[], valueType const y[], sizeType n )
    { allocate( x, y, n ) ; build() ; }

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
     * \param x vector of x-coordinates
     * \param y vector of y-coordinates
     * \param n total number of points
     */
    virtual
    void
    build ( valueType const x[], valueType const y[], sizeType n )
    { allocate( x, y, n ) ; build() ; }

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
    virtual
    void
    build( valueType const x[], valueType const y[], sizeType n ) {
      SPLINE_ASSERT( n > 1,
                     "LinearSpline::build, n =" << n <<
                     " not enought point to define a spline\n") ;
      npts = n ;
      X . resize(n) ;
      Y . resize(n) ;
      std::copy( x, x+n, X . begin() ) ;
      std::copy( y, y+n, Y . begin() ) ;
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
      if ( x < X.front() ) return Y . front() ;
      if ( x > X.back()  ) return Y . back() ;
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
    virtual void writeToStream( ostream & s ) const ;

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
      X . clear() ; X . push_back( x0 ) ;
      Y . clear() ;
      npts = 1 ;
    }

    //! given x and y vectors build a piecewise constants spline
    virtual
    void
    build( valueType const x[], valueType const y[], sizeType n ) {
      SPLINE_ASSERT( n > 1,
                     "ConstantsSpline::build, n =" << n <<
                     " not enought point to define a spline\n") ;
      npts = n ;
      X . resize(n) ;
      Y . resize(n-1) ;
      std::copy( x, x+n,   X . begin() ) ;
      std::copy( y, y+n-1, Y . begin() ) ;
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
      if ( x < X.front() ) return Y . front() ;
      if ( x > X.back()  ) return Y . back() ;
      return Y[search(x)] ;
    }

    //! First derivative
    virtual valueType D( valueType ) const { return 0 ; }
    
    //! Second derivative
    virtual valueType DD( valueType ) const { return 0 ; }

    //! Third derivative
    virtual valueType DDD( valueType ) const { return 0 ; }

    //! Print spline coefficients
    virtual void writeToStream( ostream & ) const ;

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
    void allocate( valueType const x[], valueType const y[], sizeType n ) ;

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
    virtual void writeToStream( ostream & s ) const ;

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
     * \param x vector of x-coordinates
     * \param y vector of y-coordinates
     * \param n total number of points
     */
    void
    build ( valueType const x[], valueType const y[], sizeType n )
    { allocate( x, y, n ) ; build() ; }

    //! Build a Monotone spline.
    /*!
     * \param x vector of x-coordinates
     * \param y vector of y-coordinates
     */
    void
    build ( VectorOfValues const & x, VectorOfValues const & y )
    { build( &x.front(), &y.front(), sizeType(x.size()) ) ; }
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
}

#endif
