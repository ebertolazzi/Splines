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

#include "Malloc.hh"

#ifdef SPLINES_USE_GENERIC_CONTAINER
#include "GenericContainer.hh"
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

//
// file: Splines
//

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
      ost << "In spline: " << name()        \
          << " line: " << __LINE__          \
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
      std::cout << "In spline: " << name() \
                << " line: " << __LINE__   \
                << " file: " << __FILE__   \
                << MSG << '\n' ;           \
    }
#endif

#ifdef DEBUG
  #define SPLINE_CHECK_NAN( PTR, MSG, DIM ) Splines::checkNaN( PTR, MSG, DIM )
#else
  #define SPLINE_CHECK_NAN( PTR, MSG, DIM )
#endif

#ifdef _MSC_VER
  #include <math.h>
  #define SPLINE_USE_ALLOCA
#endif

//! Various kind of splines
namespace Splines {

  using namespace ::std ; // load standard namspace

  typedef double   valueType ; //!< Floating point type for splines
  typedef unsigned sizeType  ; //!< Unsigned integer type for splines
  typedef int      indexType ; //!< Signed integer type for splines

  //! Associate a number for each type of splines implemented
  typedef enum { CONSTANT_TYPE,
                 LINEAR_TYPE,
                 AKIMA_TYPE,
                 BESSEL_TYPE,
                 PCHIP_TYPE,
                 CUBIC_TYPE,
                 QUINTIC_TYPE,
                 SPLINE_SET_TYPE } SplineType ;
  
  extern char const *spline_type[] ;

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
  bilinear3( valueType const p[4],
             valueType const M[4][4],
             valueType const q[4] ) ;

  valueType
  bilinear5( valueType const p[6],
             valueType const M[6][6],
             valueType const q[6] ) ;

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
  
    string _name ;
    bool   _check_range ;
    bool   _repeated_points_continuity ;

    sizeType  npts, npts_reserved ;
    valueType *X ; // allocated in the derived class!
    valueType *Y ; // allocated in the derived class!
    
    // for bbox
    valueType Y_min ;
    valueType Y_max ;

    sizeType search( valueType x ) const ;
    mutable sizeType lastInterval ;

    Spline(Spline const &) ; // block copy constructor
    Spline const & operator = (Spline const &) ; // block copy method

  public:

    //! spline constructor
    Spline( string const & name = "Spline", bool ck = false, bool rp = false )
    : _name(name)
    , _check_range(ck)
    , _repeated_points_continuity(rp)
    , npts(0)
    , npts_reserved(0)
    , X(nullptr)
    , Y(nullptr)
    , lastInterval(0)
    {}

    //! spline destructor
    virtual 
    ~Spline()
    {}
    
    string const & name() const { return _name ; }

    void setCheckRange( bool ck ) { _check_range = ck ; }
    bool getCheckRange() const { return _check_range ; }

    void setRepeatedPointsContinuity( bool rp ) { _repeated_points_continuity = rp ; }
    bool getRepeatedPointsContinuity() const { return _repeated_points_continuity ; }

    //! return the number of support points of the spline.
    sizeType numPoints(void) const { return npts ; }

    //! return the i-th node of the spline (x component).
    valueType xNode( sizeType i ) const { return X[i] ; }

    //! return the i-th node of the spline (y component).
    valueType yNode( sizeType i ) const { return Y[i] ; }

    //! return first node of the spline (x component).
    valueType xBegin() const { return X[0] ; }

    //! return first node of the spline (y component).
    valueType yBegin() const { return Y[0] ; }

    //! return last node of the spline (x component).
    valueType xEnd() const { return X[npts-1] ; }

    //! return last node of the spline (y component).
    valueType yEnd() const { return Y[npts-1] ; }

    //! Allocate memory for `npts` points
    virtual void reserve( sizeType npts ) = 0 ;

    //! Add a support point (x,y) to the spline.
    void pushBack( valueType x, valueType y ) ;

    //! Drop a support point to the spline.
    void dropBack() {
      if ( npts > 0 ) --npts ;
    }

    // must be defined in derived classes
    virtual
    void
    build (void) = 0 ;

    #ifdef SPLINES_USE_GENERIC_CONTAINER
    virtual void build ( GC::GenericContainer const & gc ) = 0 ;
    #endif

    //! Build a spline.
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
            sizeType n ) = 0 ;

    //! Build a spline.
    /*!
     * \param x vector of x-coordinates
     * \param y vector of y-coordinates
     * \param n total number of points
     */
    virtual
    void
    build ( valueType const x[], valueType const y[], sizeType n ) = 0 ;

    //! Build a spline.
    /*!
     * \param x vector of x-coordinates
     * \param y vector of y-coordinates
     */
    virtual
    void
    build ( vector<valueType> const & x, vector<valueType> const & y ) = 0 ;

    //! Cancel the support points, empty the spline.
    virtual
    void
    clear(void) = 0 ;

    //! return x-minumum spline value
    valueType xMin() const { return X[0] ; }

    //! return x-maximum spline value
    valueType xMax() const { return X[npts-1] ; }

    //! return y-minumum spline value
    valueType yMin() const { return Y_min ; }

    //! return y-maximum spline value
    valueType yMax() const { return Y_max ; }

    ///////////////////////////////////////////////////////////////////////////
    //! change X-origin of the spline
    void setOrigin( valueType x0 ) ;

    //! change X-range of the spline
    void setRange( valueType xmin, valueType xmax ) ;

    ///////////////////////////////////////////////////////////////////////////
    //! dump a sample of the spline
    void dump( ostream & s, sizeType nintervals, char const header[] = "x\ty" ) const ;

    void
    dump( char const fname[], sizeType nintervals, char const header[] = "x\ty" ) const
    { ofstream file(fname) ; dump( file, nintervals, header ) ; file.close() ; }

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
    char const * type_name() const { return spline_type[type()] ; }

    //! Return spline type (as number)
    virtual unsigned type() const = 0 ;

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
    Malloc<valueType> baseValue ;
    mutable valueType base[4] ;
    mutable valueType base_D[4] ;
    mutable valueType base_DD[4] ;
    mutable valueType base_DDD[4] ;

    valueType         * Yp ;
    bool              _external_alloc ;

  public:

    using Spline::build ;
  
    //! spline constructor
    CubicSplineBase( string const & name = "Spline", bool ck = false )
    : Spline(name,ck)
    , baseValue(name+"_memory")
    , Yp(nullptr)
    , _external_alloc(false)
    {}
    
    virtual
    ~CubicSplineBase()
    {}

    void copySpline( CubicSplineBase const & S ) ;

    //! return the i-th node of the spline (y' component).
    valueType ypNode( sizeType i ) const { return Yp[i] ; }

    //! change X-range of the spline
    void setRange( valueType xmin, valueType xmax ) ;

    //! Use externally allocated memory for `npts` points
    void reserve_external( sizeType n, valueType *& p_x, valueType *& p_y, valueType *& p_dy ) ;

    // --------------------------- VIRTUALS -----------------------------------

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

    // --------------------------- VIRTUALS -----------------------------------

    //! Allocate memory for `npts` points
    virtual void reserve( sizeType npts ) ;

    // must be defined in derived classes
    virtual
    void
    build (void) = 0 ;

    #ifdef SPLINES_USE_GENERIC_CONTAINER
    virtual void build ( GC::GenericContainer const & gc ) = 0 ;
    #endif

    //! Build a spline.
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
            sizeType n ) ;

    //! Build a spline.
    /*!
     * \param x vector of x-coordinates
     * \param y vector of y-coordinates
     * \param n total number of points
     */
    virtual
    void
    build ( valueType const x[], valueType const y[], sizeType n ) ;

    //! Build a spline.
    /*!
     * \param x vector of x-coordinates
     * \param y vector of y-coordinates
     */
    virtual
    void
    build ( vector<valueType> const & x, vector<valueType> const & y ) ;

    //! Cancel the support points, empty the spline.
    virtual
    void
    clear(void) ;

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
  public:

    using Spline::build ;
    using CubicSplineBase::reserve ;

    //! spline constructor
    AkimaSpline( string const & name = "Spline", bool ck = false )
    : CubicSplineBase( name, ck )
    {}

    //! spline destructor
    virtual
    ~AkimaSpline()
    {}

    //! Return spline type (as number)
    virtual unsigned type() const { return AKIMA_TYPE ; }

    // --------------------------- VIRTUALS -----------------------------------

    //! Build an Akima spline from previously inserted points
    virtual
    void
    build (void) ;

    #ifdef SPLINES_USE_GENERIC_CONTAINER
    virtual void build ( GC::GenericContainer const & gc ) ;
    #endif

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

    using Spline::build ;
    using CubicSplineBase::reserve ;

    //! spline constructor
    BesselSpline( string const & name = "Spline", bool ck = false )
    : CubicSplineBase( name, ck )
    {}

    //! spline destructor
    virtual
    ~BesselSpline()
    {}

    //! Return spline type (as number)
    virtual unsigned type() const { return BESSEL_TYPE ; }

    // --------------------------- VIRTUALS -----------------------------------

    //! Build a Bessel spline from previously inserted points
    virtual
    void
    build (void) ;

    #ifdef SPLINES_USE_GENERIC_CONTAINER
    virtual void build ( GC::GenericContainer const & gc ) ;
    #endif

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

    using Spline::build ;
    using CubicSplineBase::reserve ;

    //! spline constructor
    CubicSpline( string const & name = "CubicSpline", bool ck = false )
    : CubicSplineBase( name, ck )
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

    //! Return spline type (as number)
    virtual unsigned type() const { return CUBIC_TYPE ; }

    // --------------------------- VIRTUALS -----------------------------------

    //! Build a Cubic spline from previously inserted points
    virtual
    void
    build (void) ;

    #ifdef SPLINES_USE_GENERIC_CONTAINER
    virtual void build ( GC::GenericContainer const & gc ) ;
    #endif

  } ;

  /*
  //   ____      _     _      ____        _ _            
  //  |  _ \ ___| |__ (_)_ __/ ___| _ __ | (_)_ __   ___ 
  //  | |_) / __| '_ \| | '_ \___ \| '_ \| | | '_ \ / _ \
  //  |  __/ (__| | | | | |_) |__) | |_) | | | | | |  __/
  //  |_|   \___|_| |_|_| .__/____/| .__/|_|_|_| |_|\___|
  //                    |_|        |_|                   
  */
  void
  pchip( valueType const X[],
         valueType const Y[],
         valueType       Yp[],
         sizeType        n ) ;

  //! Pchip (Piecewise Cubic Hermite Interpolating Polynomial) spline class
  class PchipSpline : public CubicSplineBase {
  public:

    using Spline::build ;
    using CubicSplineBase::reserve ;

    //! spline constructor
    PchipSpline( string const & name = "PchipSpline", bool ck = false )
    : CubicSplineBase( name, ck )
    {}

    //! spline destructor
    virtual
    ~PchipSpline()
    {}

    //! Return spline type (as number)
    virtual unsigned type() const { return PCHIP_TYPE ; }

    // --------------------------- VIRTUALS -----------------------------------

    //! Build a Monotone spline from previously inserted points
    virtual
    void
    build (void) ;

    #ifdef SPLINES_USE_GENERIC_CONTAINER
    virtual void build ( GC::GenericContainer const & gc ) ;
    #endif

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
    Malloc<valueType> baseValue ;
    bool              _external_alloc ;
  public:

    LinearSpline( string const & name = "LinearSpline", bool ck = false )
    : Spline(name,ck)
    , baseValue( name+"_memory")
    , _external_alloc(false)
    {}

    virtual
    ~LinearSpline()
    {}

    //! Use externally allocated memory for `npts` points
    void reserve_external( sizeType n, valueType *& p_x, valueType *& p_y ) ;

    // --------------------------- VIRTUALS -----------------------------------

    //! Evalute spline value at `x`
    virtual
    valueType
    operator () ( valueType x ) const {
      if ( x < X[0]      ) return Y[0] ;
      if ( x > X[npts-1] ) return Y[npts-1] ;
      sizeType i = search(x) ;
      valueType s = (x-X[i])/(X[i+1] - X[i]) ;
      return (1-s)*Y[i] + s * Y[i+1] ;
    }

    //! First derivative
    virtual
    valueType
    D( valueType x ) const {
      if ( x < X[0]      ) return 0 ;
      if ( x > X[npts-1] ) return 0 ;
      sizeType i = search(x) ;
      return ( Y[i+1] - Y[i] ) / ( X[i+1] - X[i] ) ;
    }

    //! Second derivative
    virtual valueType DD( valueType ) const { return 0 ; }

    //! Third derivative
    virtual valueType DDD( valueType ) const { return 0 ; }

    //! Print spline coefficients
    virtual void writeToStream( std::basic_ostream<char> & s ) const ;

    //! Return spline type (as number)
    virtual unsigned type() const { return LINEAR_TYPE ; }

    // --------------------------- VIRTUALS -----------------------------------

    //! Allocate memory for `npts` points
    virtual void reserve( sizeType npts ) ;

    //! added for compatibility with cubic splines
    virtual
    void
    build()
    {}

    #ifdef SPLINES_USE_GENERIC_CONTAINER
    virtual void build ( GC::GenericContainer const & gc ) ;
    #endif

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
            sizeType n ) ;

    //! given x and y vectors build a linear spline
    /*!
     * \param x vector of x-coordinates
     * \param y vector of y-coordinates
     * \param n total number of points
     */
    virtual
    void
    build ( valueType const x[], valueType const y[], sizeType n ) ;

    //! given x and y vectors build a linear spline
    /*!
     * \param x vector of x-coordinates
     * \param y vector of y-coordinates
     */
    virtual
    void
    build ( vector<valueType> const & x, vector<valueType> const & y ) ;

    //! Cancel the support points, empty the spline.
    virtual
    void
    clear(void) ;

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
  class ConstantSpline : public Spline {
    Malloc<valueType> baseValue ;
    bool              _external_alloc ;
  public:

    ConstantSpline( string const & name = "ConstantSpline", bool ck = false )
    : Spline(name,ck)
    , baseValue(name+"_memory")
    , _external_alloc(false)
    {}

    ~ConstantSpline() {}

    //! Use externally allocated memory for `npts` points
    void reserve_external( sizeType n, valueType *& p_x, valueType *& p_y ) ;

    // --------------------------- VIRTUALS -----------------------------------


#if 0
    //! Cancel the support points, empty the spline. `x0` is the initial abscissa of empty spline
    virtual
    void
    clear( valueType x0, sizeType n_reserved = 0 ) {
      Spline::clear() ;
      Spline::reserve( n_reserved ) ;
      Spline::pushBack( x0, 0 ) ;
    }
#endif

    //! Evalute spline value at `x`
    virtual valueType operator () ( valueType x ) const ;

    //! First derivative
    virtual valueType D( valueType ) const { return 0 ; }
    
    //! Second derivative
    virtual valueType DD( valueType ) const { return 0 ; }

    //! Third derivative
    virtual valueType DDD( valueType ) const { return 0 ; }

    //! Print spline coefficients
    virtual void writeToStream( std::basic_ostream<char> & ) const ;

    //! Return spline type (as number)
    virtual unsigned type() const { return CONSTANT_TYPE ; }

    // --------------------------- VIRTUALS -----------------------------------

    //! Allocate memory for `npts` points
    virtual void reserve( sizeType npts ) ;

    //! added for compatibility with cubic splines
    virtual
    void
    build()
    {}

    #ifdef SPLINES_USE_GENERIC_CONTAINER
    virtual void build ( GC::GenericContainer const & gc ) ;
    #endif

    //! given x and y vectors build a piecewise constants spline
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
            sizeType n ) ;

    //! given x and y vectors build a linear spline
    /*!
     * \param x vector of x-coordinates
     * \param y vector of y-coordinates
     * \param n total number of points
     */
    virtual
    void
    build ( valueType const x[], valueType const y[], sizeType n ) ;

    //! given x and y vectors build a linear spline
    /*!
     * \param x vector of x-coordinates
     * \param y vector of y-coordinates
     */
    virtual
    void
    build ( vector<valueType> const & x, vector<valueType> const & y ) ;

    //! Cancel the support points, empty the spline.
    virtual
    void
    clear(void) ;

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
    Malloc<valueType> baseValue ;
    mutable valueType base[6] ;
    mutable valueType base_D[6] ;
    mutable valueType base_DD[6] ;
    mutable valueType base_DDD[6] ;
    mutable valueType base_DDDD[6] ;
    mutable valueType base_DDDDD[6] ;

    valueType * Yp ;
    valueType * Ypp ;
    bool      _external_alloc ;

  public:

    using Spline::build ;

    //! spline constructor
    QuinticSplineBase( string const & name = "Spline", bool ck = false )
    : Spline(name,ck)
    , baseValue(name+"_memeory")
    , Yp(nullptr)
    , Ypp(nullptr)
    , _external_alloc(false)
    {}
    
    virtual
    ~QuinticSplineBase()
    {}

    void copySpline( QuinticSplineBase const & S ) ;

    //! return the i-th node of the spline (y' component).
    valueType ypNode( sizeType i ) const { return Yp[i] ; }

    //! return the i-th node of the spline (y'' component).
    valueType yppNode( sizeType i ) const { return Ypp[i] ; }

    //! change X-range of the spline
    void setRange( valueType xmin, valueType xmax ) ;

    //! Use externally allocated memory for `npts` points
    void reserve_external( sizeType     n,
                           valueType *& p_x,
                           valueType *& p_y,
                           valueType *& p_Yp,
                           valueType *& p_Ypp ) ;

    // --------------------------- VIRTUALS -----------------------------------

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

    //! Return spline type (as number)
    virtual unsigned type() const { return QUINTIC_TYPE ; }


    // --------------------------- VIRTUALS -----------------------------------

    //! Allocate memory for `npts` points
    virtual void reserve( sizeType npts ) ;

    // must be defined in derived classes
    virtual
    void
    build (void) = 0 ;

    //! Build a spline.
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
            sizeType n ) ;

    //! Build a spline.
    /*!
     * \param x vector of x-coordinates
     * \param y vector of y-coordinates
     * \param n total number of points
     */
    virtual
    void
    build ( valueType const x[], valueType const y[], sizeType n ) ;

    //! Build a spline.
    /*!
     * \param x vector of x-coordinates
     * \param y vector of y-coordinates
     */
    virtual
    void
    build ( vector<valueType> const & x, vector<valueType> const & y ) ;

    //! Cancel the support points, empty the spline.
    virtual
    void
    clear(void) ;

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

    using Spline::build ;
    using QuinticSplineBase::reserve ;

    //! spline constructor
    QuinticSpline( string const & name = "Spline", bool ck = false )
    : QuinticSplineBase( name, ck )
    {}

    //! spline destructor
    virtual
    ~QuinticSpline()
    {}

    // --------------------------- VIRTUALS -----------------------------------

    //! Build a Monotone quintic spline from previously inserted points
    virtual
    void
    build (void) ;

    #ifdef SPLINES_USE_GENERIC_CONTAINER
    virtual void build ( GC::GenericContainer const & gc ) ;
    #endif

  } ;

  /*
  //   ____        _ _            ____       _
  //  / ___| _ __ | (_)_ __   ___/ ___|  ___| |_
  //  \___ \| '_ \| | | '_ \ / _ \___ \ / _ \ __|
  //   ___) | |_) | | | | | |  __/___) |  __/ |_
  //  |____/| .__/|_|_|_| |_|\___|____/ \___|\__|
  //        |_|
  */

  //! Splines Management Class
  class SplineSet {

    SplineSet(SplineSet const &) ; // block copy constructor
    SplineSet const & operator = (SplineSet const &) ; // block copy method

  protected:

    string const _name ;
    
    sizeType _npts ;
    sizeType _nspl ;

    Malloc<valueType>  baseValue ;
    Malloc<valueType*> basePointer ;

    valueType *_X ;
    valueType **_Y ;
    valueType **_Yp ;
    valueType **_Ypp ;
    valueType *_Ymin ;
    valueType *_Ymax ;

    mutable sizeType lastInterval ;
    sizeType search( valueType x ) const ;
    
    vector<Spline*> splines ;

  public:

    //! spline constructor
    SplineSet( string const & name = "Splines" )
    : _name(name)
    , baseValue(name+"_values")
    , basePointer(name+"_pointers")
    , _X(nullptr)
    , _Y(nullptr)
    , _Yp(nullptr)
    , _Ypp(nullptr)
    , _Ymin(nullptr)
    , _Ymax(nullptr)
    , lastInterval(0)
    {}

    //! spline destructor
    virtual 
    ~SplineSet()
    { baseValue.free() ; basePointer.free() ; }

    string const & name() const { return _name ; }
    string const & header( sizeType const i ) const { return splines[i]->name() ; }

    //! return the number of support points of the splines
    sizeType numPoints(void) const { return _npts ; }

    //! return the number splines in the spline set
    sizeType numSplines(void) const { return _nspl ; }

    //! return the vector of values of x-nodes
    valueType const * xNodes() const { return _X ; }

    //! return the vector of values of x-nodes
    valueType const * yNodes( sizeType i ) const {
      SPLINE_ASSERT( i >=0 && i < _nspl,
                     "SplineSet::yNodes( " << i << ") argument out of range [0," << _nspl-1 << "]" ) ;
      return _Y[i] ;
    }

    //! return the i-th node of the spline (x component).
    valueType xNode( sizeType npt ) const { return _X[npt] ; }

    //! return the i-th node of the spline (y component).
    valueType yNode( sizeType npt, sizeType spl ) const { return _Y[spl][npt] ; }

    //! return x-minumum spline value
    valueType xMin() const { return _X[0] ; }

    //! return x-maximum spline value
    valueType xMax() const { return _X[_npts-1] ; }

    //! return y-minumum spline value
    valueType yMin( sizeType spl ) const { return _Ymin[spl] ; }

    //! return y-maximum spline value
    valueType yMax( sizeType spl ) const { return _Ymax[spl] ; }

    //! Return pointer to the `i`-th spline
    Spline * getSpline( sizeType i ) const {
      SPLINE_ASSERT( i < _nspl, "SplineSet::getSpline( " << i << ") argument out of range [0," << _nspl-1 << "]" ) ;
      return splines[i] ;
    }

    //! Evalute spline value
    valueType operator () ( valueType x, sizeType spl ) const { return (*getSpline(spl))(x) ; }

    //! First derivative
    valueType D( valueType x, sizeType spl ) const { return getSpline(spl)->D(x) ; }

    //! Second derivative
    valueType DD( valueType x, sizeType spl ) const { return getSpline(spl)->DD(x) ; }

    ///////////////////////////////////////////////////////////////////////////
    /*! Build a set of splines
     * \param nspl       the number of splines
     * \param nspl       the number of splines
     * \param headers    the names of the splines
     * \param stype      the type of each spline
     * \param rp_policy  treatment of repeated points of eachj splines. True = continuous,
     *                   rp_policy = NULL all splines are discontinuous at repeated points
     * \param X          pointer to X independent values
     * \param Y          vector of `nspl` pointers to Y depentendent values.
     */

    void
    build ( indexType  const nspl,
            indexType  const npts,
            char       const *headers[],
            SplineType const stype[],
            valueType  const X[],
            valueType  const *Y[],
            bool       const rp_policy[] = nullptr ) ;

    #ifdef SPLINES_USE_GENERIC_CONTAINER
    void build ( GC::GenericContainer const & gc ) ;
    #endif

    //! Return spline type (as number)
    virtual unsigned type() const { return SPLINE_SET_TYPE ; }

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

    SplineSurf(SplineSurf const &) ; // block copy constructor
    SplineSurf const & operator = (SplineSurf const &) ; // block copy method

  protected:
  
    string const _name ;
    bool         _check_range ;

    vector<valueType> X, Y, Z ;
    
    valueType Z_min, Z_max ;

    mutable sizeType lastInterval_x ;
    sizeType search_x( valueType x ) const ;

    mutable sizeType lastInterval_y ;
    sizeType search_y( valueType y ) const ;

    sizeType ipos_C( sizeType i, sizeType j, sizeType ldZ ) const { return i*ldZ + j ; }
    sizeType ipos_F( sizeType i, sizeType j, sizeType ldZ ) const { return i + ldZ*j ; }

    sizeType ipos_C( sizeType i, sizeType j ) const { return ipos_C(i,j,sizeType(Y.size())) ; }
    sizeType ipos_F( sizeType i, sizeType j ) const { return ipos_F(i,j,sizeType(X.size())) ; }

    virtual void makeSpline() = 0 ;

  public:

    //! spline constructor
    SplineSurf( string const & name = "Spline", bool ck = false )
    : _name(name)
    , _check_range(ck)
    , X()
    , Y()
    , Z()
    , Z_min(0)
    , Z_max(0)
    , lastInterval_x(0)
    , lastInterval_y(0)
    {}

    //! spline destructor
    virtual 
    ~SplineSurf()
    {}

    string const & name() const { return _name ; }

    void setCheckRange( bool ck ) { _check_range = ck ; }
    bool getCheckRange() const { return _check_range ; }

    //! Cancel the support points, empty the spline.
    virtual
    void
    clear (void) {
      X.clear() ;
      Y.clear() ;
      Z.clear() ;
      Z_min = Z_max = 0 ;
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
    valueType zNode( sizeType i, sizeType j ) const { return Z[ipos_C(i,j)] ; }

    //! return x-minumum spline value
    valueType xMin() const { return X.front() ; }

    //! return x-maximum spline value
    valueType xMax() const { return X.back()  ; }

    //! return y-minumum spline value
    valueType yMin() const { return Y.front() ; }

    //! return y-maximum spline value
    valueType yMax() const { return Y.back()  ; }

    //! return z-minumum spline value
    valueType zMin() const { return Z_min ; }

    //! return z-maximum spline value
    valueType zMax() const { return Z_max ; }

    ///////////////////////////////////////////////////////////////////////////
    /*! Build surface spline
     * \param x       vector of x-coordinates
     * \param incx    access elements as x[0], x[incx], x[2*incx],...
     * \param y       vector of y-coordinates
     * \param incy    access elements as y[0], y[incy], x[2*incy],...
     * \param z       matrix of z-values
     * \param ldZ     leading dimension of the matrix. Elements are stored
     *                by row Z(i,j) = z[i*ldZ+j] as C-matrix
     * \param nx      total number of points in direction x
     * \param ny      total number of points in direction y
     * \param fortran_storage if true elements are stored by column
     *                        i.e. Z(i,j) = z[i+j*ldZ] as Fortran-matrix
     * \param transposed      if true matrix Z is stored transposed
     */
    void
    build ( valueType const x[], sizeType incx,
            valueType const y[], sizeType incy,
            valueType const z[], sizeType ldZ,
            sizeType nx, sizeType ny,
            bool fortran_storage = false,
            bool transposed      = false ) ;

    /*! Build surface spline
     * \param x       vector of x-coordinates, nx = x.size()
     * \param y       vector of y-coordinates, ny = y.size()
     * \param z       matrix of z-values. Elements are stored
     *                by row Z(i,j) = z[i*ny+j] as C-matrix
     * \param fortran_storage if true elements are stored by column
     *                        i.e. Z(i,j) = z[i+j*nx] as Fortran-matrix
     * \param transposed      if true matrix Z is stored transposed
     */
    void
    build ( vector<valueType> const & x,
            vector<valueType> const & y,
            vector<valueType> const & z,
            bool fortran_storage = false,
            bool transposed      = false ) {
      sizeType nx = sizeType(x.size()) ;
      sizeType ny = sizeType(y.size()) ;
      valueType const * xx = &x.front() ;
      valueType const * yy = &y.front() ;
      valueType const * zz = &z.front() ;
      if ( fortran_storage ) build ( xx, 1, yy, 1, zz, nx, nx, ny, fortran_storage, transposed ) ;
      else                   build ( xx, 1, yy, 1, zz, ny, nx, ny, fortran_storage, transposed ) ;
    }

    /*! Build surface spline
     * \param z               matrix of z-values. Elements are stored
     *                        by row Z(i,j) = z[i*ny+j] as C-matrix
     * \param ldZ             leading dimension of the matrix. Elements are stored
     *                        by row Z(i,j) = z[i*ldZ+j] as C-matrix
     * \param fortran_storage if true elements are stored by column
     *                        i.e. Z(i,j) = z[i+j*nx] as Fortran-matrix
     * \param transposed      if true matrix Z is stored transposed
     */
    void
    build ( valueType const z[], sizeType ldZ, sizeType nx, sizeType ny,
            bool fortran_storage = false,
            bool transposed      = false ) ;

    /*! Build surface spline
     * \param z               matrix of z-values. Elements are stored
     *                        by row Z(i,j) = z[i*ny+j] as C-matrix.
     *                        ldZ leading dimension of the matrix is ny for C-storage
     *                        and nx for Fortran storage.
     * \param fortran_storage if true elements are stored by column
     *                        i.e. Z(i,j) = z[i+j*nx] as Fortran-matrix
     * \param transposed      if true matrix Z is stored transposed
     */
    void
    build ( vector<valueType> const & z, sizeType nx, sizeType ny,
            bool fortran_storage = false,
            bool transposed      = false ) {
      if ( fortran_storage ) build ( &z.front(), nx, nx, ny, fortran_storage, transposed ) ;
      else                   build ( &z.front(), ny, nx, ny, fortran_storage, transposed ) ;
    }

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
    virtual void makeSpline() {}
  public:
  
    //! spline constructor
    BilinearSpline( string const & name = "Spline", bool ck = false )
    : SplineSurf(name,ck)
    {}
    
    virtual
    ~BilinearSpline()
    {}

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
  //! Bi-cubic spline base class
  class BiCubicSplineBase : public SplineSurf {
  protected:
  
    //! \cond PRIVATE

    vector<valueType> DX, DY, DXY ;
    mutable valueType u[4] ;
    mutable valueType u_D[4] ;
    mutable valueType u_DD[4] ;
    mutable valueType v[4] ;
    mutable valueType v_D[4] ;
    mutable valueType v_DD[4] ;
    mutable valueType bili3[4][4] ;

    void load( sizeType i, sizeType j ) const ;

    //! \endcond

  public:
  
    //! spline constructor
    BiCubicSplineBase( string const & name = "Spline", bool ck = false )
    : SplineSurf( name, ck )
    , DX()
    , DY()
    {}
    
    virtual
    ~BiCubicSplineBase()
    {}

    valueType DxNode ( sizeType i, sizeType j ) const { return DX[ipos_C(i,j)] ; }
    valueType DyNode ( sizeType i, sizeType j ) const { return DY[ipos_C(i,j)] ; }
    valueType DxyNode( sizeType i, sizeType j ) const { return DXY[ipos_C(i,j)] ; }

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
    BiCubicSpline( string const & name = "Spline", bool ck = false )
    : BiCubicSplineBase( name, ck )
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
    Akima2Dspline( string const & name = "Spline", bool ck = false )
    : BiCubicSplineBase( name, ck )
    {}
    
    virtual
    ~Akima2Dspline()
    {}

    //! Print spline coefficients
    virtual void writeToStream( ostream & s ) const ;

    //! Return spline typename
    virtual char const * type_name() const ;

  } ;
  
  /*
  //   ____  _  ___        _       _   _      ____        _ _            ____
  //  | __ )(_)/ _ \ _   _(_)_ __ | |_(_) ___/ ___| _ __ | (_)_ __   ___| __ )  __ _ ___  ___
  //  |  _ \| | | | | | | | | '_ \| __| |/ __\___ \| '_ \| | | '_ \ / _ \  _ \ / _` / __|/ _ \
  //  | |_) | | |_| | |_| | | | | | |_| | (__ ___) | |_) | | | | | |  __/ |_) | (_| \__ \  __/
  //  |____/|_|\__\_\\__,_|_|_| |_|\__|_|\___|____/| .__/|_|_|_| |_|\___|____/ \__,_|___/\___|
  //                                               |_|
  */
  //! Bi-quintic spline base class
  class BiQuinticSplineBase : public SplineSurf {
  protected:
  
    //! \cond PRIVATE

    vector<valueType> DX, DXX, DY, DYY, DXY, DXYY, DXXY, DXXYY ;
    mutable valueType u[6] ;
    mutable valueType u_D[6] ;
    mutable valueType u_DD[6] ;
    mutable valueType v[6] ;
    mutable valueType v_D[6] ;
    mutable valueType v_DD[6] ;
    mutable valueType bili5[6][6] ;

    void load( sizeType i, sizeType j ) const ;

    //! \endcond

  public:
  
    //! spline constructor
    BiQuinticSplineBase( string const & name = "Spline", bool ck = false )
    : SplineSurf( name, ck )
    , DX()
    , DXX()
    , DY()
    , DYY()
    , DXY()
    {}
    
    virtual
    ~BiQuinticSplineBase()
    {}

    valueType DxNode ( sizeType i, sizeType j ) const { return DX[ipos_C(i,j)] ; }
    valueType DyNode ( sizeType i, sizeType j ) const { return DY[ipos_C(i,j)] ; }
    valueType DxxNode( sizeType i, sizeType j ) const { return DXX[ipos_C(i,j)] ; }
    valueType DyyNode( sizeType i, sizeType j ) const { return DYY[ipos_C(i,j)] ; }
    valueType DxyNode( sizeType i, sizeType j ) const { return DXY[ipos_C(i,j)] ; }

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
  } ;

  /*
  //   ____  _  ___        _       _   _      ____        _ _
  //  | __ )(_)/ _ \ _   _(_)_ __ | |_(_) ___/ ___| _ __ | (_)_ __   ___
  //  |  _ \| | | | | | | | | '_ \| __| |/ __\___ \| '_ \| | | '_ \ / _ \ 
  //  | |_) | | |_| | |_| | | | | | |_| | (__ ___) | |_) | | | | | |  __/ 
  //  |____/|_|\__\_\\__,_|_|_| |_|\__|_|\___|____/| .__/|_|_|_| |_|\___|
  //                                               |_|
  */
  //! cubic spline base class
  class BiQuinticSpline : public BiQuinticSplineBase {
    virtual void makeSpline() ;

  public:
  
    //! spline constructor
    BiQuinticSpline( string const & name = "Spline", bool ck = false )
    : BiQuinticSplineBase( name, ck )
    {}
    
    virtual
    ~BiQuinticSpline()
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
  using Splines::ConstantSpline ;
  using Splines::QuinticSpline ;

  using Splines::BilinearSpline ;
  using Splines::BiCubicSpline ;
  using Splines::BiQuinticSpline ;
  using Splines::Akima2Dspline ;

  using Splines::SplineSet ;
  using Splines::SplineType ;

}

#endif
