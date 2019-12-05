
  /*\
   |   ____                  _ _
   |  | __ )       ___ _ __ | (_)_ __   ___
   |  |  _ \ _____/ __| '_ \| | | '_ \ / _ \
   |  | |_) |_____\__ \ |_) | | | | | |  __/
   |  |____/      |___/ .__/|_|_|_| |_|\___|
   |                  |_|
  \*/

  //! B-spline base class
  template <size_t _degree>
  class BSpline : public Spline {
  protected:
    SplineMalloc<real_type> baseValue;
    real_type * knots;
    real_type * yPolygon;
    bool        _external_alloc;

    integer knot_search( real_type x ) const;

    // extrapolation
    real_type s_L,   s_R;
    real_type ds_L,  ds_R;
    real_type dds_L, dds_R;

  public:

    using Spline::build;

    //! spline constructor
    BSpline( string const & name = "BSpline" )
    : Spline(name)
    , baseValue(name+"_memory")
    , knots(nullptr)
    , yPolygon(nullptr)
    , _external_alloc(false)
    {}

    virtual
    ~BSpline() SPLINES_OVERRIDE
    {}

    void
    copySpline( BSpline const & S );

    //! return the i-th node of the spline (y' component).
    real_type
    yPoly( integer i ) const
    { return yPolygon[size_t(i)]; }

    //! change X-range of the spline
    void
    setRange( real_type xmin, real_type xmax );

    //! Use externally allocated memory for `npts` points
    void
    reserve_external( integer      n,
                      real_type *& p_x,
                      real_type *& p_y,
                      real_type *& p_knots,
                      real_type *& p_yPolygon );

    // --------------------------- VIRTUALS -----------------------------------

    //! Return spline type (as number)
    virtual
    unsigned
    type() const SPLINES_OVERRIDE
    { return BSPLINE_TYPE; }

    //! Evaluate spline value
    virtual
    real_type
    operator () ( real_type x ) const SPLINES_OVERRIDE;

    //! First derivative
    virtual
    real_type
    D( real_type x ) const SPLINES_OVERRIDE;

    //! Second derivative
    virtual
    real_type
    DD( real_type x ) const SPLINES_OVERRIDE;

    //! Third derivative
    virtual
    real_type
    DDD( real_type x ) const SPLINES_OVERRIDE;

    //! Print spline coefficients
    virtual
    void
    writeToStream( ostream_type & s ) const SPLINES_OVERRIDE;

    // --------------------------- VIRTUALS -----------------------------------

    //! Allocate memory for `npts` points
    virtual
    void
    reserve( integer npts ) SPLINES_OVERRIDE;

    virtual
    void
    build(void) SPLINES_OVERRIDE;

    //! Cancel the support points, empty the spline.
    virtual
    void
    clear(void) SPLINES_OVERRIDE;

    //! get the piecewise polinomials of the spline
    virtual
    integer // order
    coeffs( real_type cfs[],
            real_type nodes[],
            bool transpose = false ) const SPLINES_OVERRIDE;

    virtual
    integer // order
    order() const SPLINES_OVERRIDE
    { return _degree+1; }

    // ---------------------------  UTILS  -----------------------------------

    //! Evaluate spline value
    static
    real_type
    eval( real_type       x,
          integer         n,
          real_type const Knots[],
          real_type const yPolygon[] );

    //! First derivative
    static
    real_type
    eval_D( real_type       x,
            integer         n,
            real_type const Knots[],
            real_type const yPolygon[] );

    //! Second derivative
    static
    real_type
    eval_DD( real_type       x,
             integer         n,
             real_type const Knots[],
             real_type const yPolygon[] );

    //! Third derivative
    static
    real_type
    eval_DDD( real_type       x,
              integer         n,
              real_type const Knots[],
              real_type const yPolygon[] );

    //! B-spline bases
    void bases( real_type x, real_type vals[] ) const;
    void bases_D( real_type x, real_type vals[] ) const;
    void bases_DD( real_type x, real_type vals[] ) const;
    void bases_DDD( real_type x, real_type vals[] ) const;

    integer bases_nz( real_type x, real_type vals[] ) const;
    integer bases_D_nz( real_type x, real_type vals[] ) const;
    integer bases_DD_nz( real_type x, real_type vals[] ) const;
    integer bases_DDD_nz( real_type x, real_type vals[] ) const;

    // Utilities
    static
    void
    knots_sequence( integer         n,
                    real_type const X[],
                    real_type     * Knots );

    static
    void
    sample_bases( integer             nx, // number of sample points
                  real_type const     X[],
                  integer             nb, // number of bases
                  real_type const     Knots[],
                  vector<integer>   & II, // GCC on linux bugged for I
                  vector<integer>   & JJ,
                  vector<real_type> & vals );

  };



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
#include <iomanip>

#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat"
#endif

/**
 * 
 */

namespace Splines {

  using namespace std; // load standard namespace


  template <size_t degree>
  class BSplineBase {
  public:

    /*
    // dati i nodi knot[0]...knot[2*degree+1] calcola le basi
    // B[0]...B[degree] che hanno supporto knot[j]..knot[j+degree+1]
    // vale la condizione knot[degree] <= x <= knot[degree+1]
    */
    static
    void
    eval( real_type       x,
          real_type const knot[],
          real_type       Bbase[] ) {
      real_type B[2*degree+2];
      std::fill( B, B+2*degree+2, 0.0 );
      B[degree] = 1;
      for ( size_t r = 1; r <= size_t(degree); ++r ) {
        for ( size_t j = 0; j <= size_t(2*degree-r); ++j ) {
          if ( knot[j] <= x && x <= knot[j+r+1] ) {
            real_type oma = 0;
            real_type omb = 0;
            if ( knot[j+r]   > knot[j]   ) oma = (x-knot[j])/(knot[j+r]-knot[j]);
            if ( knot[j+r+1] > knot[j+1] ) omb = (knot[j+r+1]-x)/(knot[j+r+1]-knot[j+1]);
            B[j] = oma*B[j] + omb*B[j+1];
          }
        }
      }
      std::copy( B, B+degree+1, Bbase );
    }

    static
    void
    eval_D( real_type       x,
            real_type const knot[],
            real_type       B_D[] ) {
      real_type B[degree+2];
      B[0] = B[degree+1] = 0;
      BSplineBase<degree-1>::eval( x, knot+1, B+1 );
      for ( size_t j = 1; j <= size_t(degree); ++j ) B[j] *= degree/(knot[j+degree]-knot[j]);
      for ( size_t j = 0; j <= size_t(degree); ++j ) B_D[j] = B[j]-B[j+1];
    }

    static
    void
    eval_DD( real_type       x,
             real_type const knot[],
             real_type       B_DD[] ) {
      real_type B_D[degree+2];
      B_D[0] = B_D[degree+1] = 0;
      BSplineBase<degree-1>::eval_D( x, knot+1, B_D+1 );
      for ( size_t j = 1; j <= size_t(degree); ++j ) B_D[j] *= degree/(knot[j+degree]-knot[j]);
      for ( size_t j = 0; j <= size_t(degree); ++j ) B_DD[j] = B_D[j]-B_D[j+1];
    }

    static
    void
    eval_DDD ( real_type       x,
               real_type const knot[],
               real_type       B_DDD[] ) {
      real_type B_DD[degree+2];
      B_DD[0] = B_DD[degree+1] = 0;
      BSplineBase<degree-1>::eval_DD( x, knot+1, B_DD+1 );
      for ( size_t j = 1; j <= size_t(degree); ++j ) B_DD[j] *= degree/(knot[j+degree]-knot[j]);
      for ( size_t j = 0; j <= size_t(degree); ++j ) B_DDD[j] = B_DD[j]-B_DD[j+1];
    }
  };


  template <>
  class BSplineBase<0> {
  public:

    static void
    eval( real_type       x,
          real_type const knot[],
          real_type       Bbase[] ) {
      Bbase[0] = knot[0] <= x && x <= knot[1] ? 1 : 0;
    }

    static void
    eval_D( real_type, real_type const [], real_type B_D[] ) { B_D[0] = 0; }

    static void
    eval_DD( real_type, real_type const [], real_type B_DD[] ) { B_DD[0] = 0; }

    static void
    eval_DDD( real_type, real_type const [], real_type B_DDD[] ) { B_DDD[0] = 0; }
  };

  template <size_t degree>
  class BSplineEval {
  public:

    /*
    //  Calcola il valore di
    //
    //  y[0] * B[0](x) + y[1] * B[1](x) + ... + y[degree] * B[degree](x)
    //
    //  Usando usando la ricorsione e l'algoritmo di De Boor.
    */
    static
    void
    eval_levels( real_type x, real_type const knot[], real_type y[] ) {
      size_t j = degree;
      do {
        if ( knot[j+degree] > knot[j] ) {
          real_type omega = (x-knot[j])/(knot[j+degree]-knot[j]);
          y[j] = (1-omega)*y[j-1]+omega*y[j];
        } else {
          y[j] = y[j-1];
        }
      } while ( --j > 0 );
      BSplineEval<degree-1>::eval_levels( x, knot+1, y+1 );
    }

    /*
    //  Calcola il valore di
    //
    //  y[0] * B[0](x) + y[1] * B[1](x) + ... + y[degree] * B[degree](x)
    //
    //  Usando usando la ricorsione e l'algoritmo di De Boor.
    //  offs e' tale che knot[degree-offs] <= x <= knot[degree-offs+1]
    */
    static
    real_type
    eval( real_type       x,
          real_type const knot[],
          real_type const y[] ) {
      real_type c[degree+1], kn[2*degree+2];
      std::copy( y,    y+degree+1,      c  );
      std::copy( knot, knot+2*degree+2, kn );
      BSplineEval<degree>::eval_levels( x, kn, c );
      return c[degree];
    }

    static
    real_type
    eval_D( real_type       x,
            real_type const knot[],
            real_type const y[] ) {
      real_type d[degree]; // poligono derivata
      for ( size_t j = 0; j < size_t(degree); ++j )
        d[j] = degree*(y[j+1]-y[j])/(knot[j+degree+1] - knot[j+1]);
      return BSplineEval<degree-1>::eval( x, knot+1, d );
    }

    static
    real_type
    eval_DD( real_type       x,
             real_type const knot[],
             real_type const y[] ) {
      real_type d[degree]; // poligono derivata
      for ( size_t j = 0; j < size_t(degree); ++j )
        d[j] = degree*(y[j+1]-y[j])/(knot[j+degree+1] - knot[j+1]);
      return BSplineEval<degree-1>::eval_D( x, knot+1, d );
    }

    static
    real_type
    eval_DDD( real_type       x,
              real_type const knot[],
              real_type const y[] ) {
      real_type d[degree]; // poligono derivata
      for ( size_t j = 0; j < size_t(degree); ++j )
        d[j] = degree*(y[j+1]-y[j])/(knot[j+degree+1] - knot[j+1]);
      return BSplineEval<degree-1>::eval_DD( x, knot+1, d );
    }
  };

  template <>
  class BSplineEval<1> {
  public:

    static void
    eval_B( real_type       x,
            real_type const knot[],
            real_type       Bbase[] ) {
      real_type B[4];
      std::fill( B, B+4, 0.0 );
      B[1] = 1;
      for ( size_t j = 0; j <= 1; ++j ) {
        if ( knot[j] <= x && x <= knot[j+2] ) {
          real_type oma = 0;
          real_type omb = 0;
          if ( knot[j+1] > knot[j]   ) oma = (x-knot[j])/(knot[j+1]-knot[j]);
          if ( knot[j+2] > knot[j+1] ) omb = (knot[j+2]-x)/(knot[j+2]-knot[j+1]);
          B[j] = oma*B[j] + omb*B[j+1];
        }
      }
      std::copy( B, B+2, Bbase );
    }

    static void
    eval_B_D( real_type, real_type const [], real_type B_D[] ) { B_D[0] = 0; }

    static void
    eval_B_DD( real_type, real_type const [], real_type B_DD[] ) { B_DD[0] = 0; }

    static void
    eval_B_DDD( real_type, real_type const [], real_type B_DDD[] ) { B_DDD[0] = 0; }

    static
    void
    eval_levels( real_type       x,
                 real_type const knot[],
                 real_type       y[] ) {
      real_type omega = (x-knot[1])/(knot[2]-knot[1]);
      y[1] = (1-omega)*y[0]+omega*y[1];
    }

    static
    real_type
    eval( real_type       x,
          real_type const knot[],
          real_type const y[] ) {
      real_type omega = (x-knot[1])/(knot[2]-knot[1]);
      return (1-omega)*y[0]+omega*y[1];
    }

    static
    real_type
    eval_D( real_type             ,
            real_type const knot[],
            real_type const y[] ) {
      return (y[1]-y[0])/(knot[2] - knot[1]);
    }

    static
    real_type
    eval_DD( real_type         ,
             real_type const [],
             real_type const [] ) {
      return 0;
    }

    static
    real_type
    eval_DDD( real_type         ,
              real_type const [],
              real_type const [] ) {
      return 0;
    }

  };
  
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // standard average nodes
  template <size_t _degree>
  void
  BSpline<_degree>::knots_sequence( integer         n,
                                    real_type const X[],
                                    real_type     * Knots ) {
    size_t nn = size_t(n) - _degree - 1;
    for ( size_t j = 0; j <= _degree; ++j ) Knots[j] = X[0];
    Knots += _degree+1;
    for ( size_t j = 0; j < nn; ++j ) {
      real_type tmp = 0;
      for ( size_t k = 1; k <= _degree; ++k ) tmp += X[j+k];
      Knots[j] = tmp / _degree;
    }
    Knots += nn;
    for ( size_t j = 0; j <= _degree; ++j ) Knots[j] = X[n-1];
  }

  template <size_t _degree>
  void
  BSpline<_degree>::sample_bases( integer             nx,
                                  real_type const     X[],
                                  integer             nb,
                                  real_type const     Knots[],
                                  vector<integer>   & II,
                                  vector<integer>   & JJ,
                                  vector<real_type> & vals ) {
    real_type row[_degree+1];
    II.clear();   II.reserve( size_t(nx) );
    JJ.clear();   JJ.reserve( size_t(nx) );
    vals.clear(); vals.reserve( size_t(nx) );
    for ( integer i = 0; i < nx; ++i ) {
      integer ii = integer(lower_bound( Knots+_degree, Knots+nb+_degree, X[i] ) - Knots );
      if ( ii > integer(_degree) ) --ii;
      ii -= _degree;
      BSplineBase<_degree>::eval( X[i], Knots+ii, row );
      for ( integer j = 0; j <= integer(_degree); ++j ) {
        II.push_back( i );
        JJ.push_back( ii+j );
        vals.push_back( row[size_t(j)] );
      }
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <size_t _degree>
  integer
  BSpline<_degree>::knot_search( real_type x ) const {
    SPLINE_ASSERT( npts > 0, "\nknot_search(" << x << ") empty spline");
    if ( x < knots[lastInterval] || knots[lastInterval+1] < x ) {
      if ( _check_range ) {
        SPLINE_ASSERT(
          x >= X[0] && x <= X[npts-1],
          "method search( " << x << " ) out of range: [" <<
          X[0] << ", " << X[npts-1] << "]"
        );
      }
      // 0 1 2 3
      lastInterval = integer(lower_bound( knots+_degree, knots+npts+1, x ) - knots);
      if ( lastInterval > integer(_degree) ) --lastInterval;
    }
    return lastInterval - integer(_degree);
  }

  /*
  // Solve banded linear system
  */
  static
  void
  solveBanded( real_type * rows,
               real_type * rhs,
               integer     n,
               integer     ndiag ) {
    // forward
    integer rsize = 1+2*ndiag;
    integer i = 0;
    real_type * rowsi = rows + ndiag;
    do {
      // pivot
      real_type pivot = rowsi[0]; // elemento sulla diagonale
      
      // scala equazione
      for ( integer j = 0; j <= ndiag; ++j ) rowsi[j] /= pivot;
      rhs[i] /= pivot;

      // azzera colonna
      integer nr = i+ndiag >= n ? n-i-1 : ndiag;
      for ( integer k = 1; k <= nr; ++k ) {
        real_type * rowsk = rowsi + k * (rsize-1);
        real_type tmp = rowsk[0];
        rowsk[0] = 0;
        for ( integer j = 1; j <= nr; ++j ) rowsk[j] -= tmp*rowsi[j];
        rhs[i+k] -= tmp*rhs[i];
      }
      rowsi += rsize;
    } while ( ++i < n );
    // backward
    while ( i > 0 ) {
      --i;
      rowsi -= rsize;
      integer nr = i+ndiag >= n ? n-i-1 : ndiag;
      for ( integer j = 1; j <= nr; ++j ) rhs[i] -= rhs[i+j]*rowsi[j];
    }
  }

  template <size_t _degree>
  void
  BSpline<_degree>::build(void) {
    std::vector<real_type> band;
    band.resize( size_t(npts*(2*integer(_degree)+1)) );
    
    std::fill( band.begin(), band.end(), 0 );
    knots_sequence( npts, X, knots );
    
    // costruzione sistema lineare

    // calcola il valore delle basi non zero quando
    // knot[degree] <= x <= knot[degree+1]
    integer nr = integer(2*_degree+1);
    for ( integer i = 0; i < npts; ++i ) {
      real_type * rowi = &band[size_t(nr * i)];
      integer ii = knot_search( X[i] );
      BSplineBase<_degree>::eval( X[i], knots+ii, (rowi + ii + _degree) - i );
      yPolygon[i] = Y[i];
    }
    solveBanded( &band.front(), yPolygon, npts, _degree );
    // extrapolation
    real_type const * knots_R    = knots    + npts - _degree - 1;
    real_type const * yPolygon_R = yPolygon + npts - _degree - 1;
    real_type x_L = X[0];
    real_type x_R = X[npts-1];
    s_L   = BSplineEval<_degree>::eval(x_L,knots,yPolygon);
    s_R   = BSplineEval<_degree>::eval(x_R,knots_R,yPolygon_R);
    ds_L  = BSplineEval<_degree>::eval_D(x_L,knots,yPolygon);
    ds_R  = BSplineEval<_degree>::eval_D(x_R,knots_R,yPolygon_R);
    dds_L = BSplineEval<_degree>::eval_DD(x_L,knots,yPolygon);
    dds_R = BSplineEval<_degree>::eval_DD(x_R,knots_R,yPolygon_R);
  }
  
  template <size_t _degree>
  void
  BSpline<_degree>::clear(void) {
    if ( !_external_alloc ) baseValue.free();
    npts = npts_reserved = 0;
    _external_alloc = false;
    X = Y = knots = yPolygon = nullptr;
  }

  template <size_t _degree>
  void
  BSpline<_degree>::reserve( integer n ) {
    if ( _external_alloc && n <= npts_reserved ) {
      // nothing to do!, already allocated
    } else {
      npts_reserved = n;
      baseValue.allocate( size_t(4*n+integer(_degree)+1) );
      X        = baseValue( size_t(n) );
      Y        = baseValue( size_t(n) );
      knots    = baseValue( size_t(n+integer(_degree)+1) );
      yPolygon = baseValue( size_t(n) );
      _external_alloc = false;
    }
    npts         = 0;
    lastInterval = _degree;
  }

  template <size_t _degree>
  void
  BSpline<_degree>::reserve_external( integer      n,
                                      real_type *& p_x,
                                      real_type *& p_y,
                                      real_type *& p_knots,
                                      real_type *& p_yPolygon ) {
    npts_reserved   = n;
    X               = p_x;
    Y               = p_y;
    knots           = p_knots;
    yPolygon        = p_yPolygon;
    npts            = 0;
    lastInterval    = _degree;
    _external_alloc = true;
  }

  template <size_t _degree>
  real_type
  BSpline<_degree>::eval( real_type       x,
                          integer         n,
                          real_type const Knots[],
                          real_type const yPoly[] ) {
    integer idx = integer(lower_bound( Knots+_degree, Knots+n+1, x ) - Knots);
    if ( idx > integer(_degree) ) --idx;
    idx -= _degree;
    return BSplineEval<_degree>::eval(x,Knots+idx,yPoly+idx);
  }

  template <size_t _degree>
  real_type
  BSpline<_degree>::eval_D( real_type       x,
                            integer         n,
                            real_type const Knots[],
                            real_type const yPoly[] ) {
    integer idx = integer(lower_bound( Knots+_degree, Knots+n+1, x ) - Knots);
    if ( idx > integer(_degree) ) --idx;
    idx -= _degree;
    return BSplineEval<_degree>::eval_D(x,Knots+idx,yPoly+idx);
  }

  template <size_t _degree>
  real_type
  BSpline<_degree>::eval_DD( real_type       x,
                             integer         n,
                             real_type const Knots[],
                             real_type const yPoly[] ) {
    integer idx = integer(lower_bound( Knots+_degree, Knots+n+1, x ) - Knots);
    if ( idx > integer(_degree) ) --idx;
    idx -= _degree;
    return BSplineEval<_degree>::eval_DD(x,Knots+idx,yPoly+idx);
  }

  template <size_t _degree>
  real_type
  BSpline<_degree>::eval_DDD( real_type       x,
                              integer         n,
                              real_type const Knots[],
                              real_type const yPoly[] ) {
    integer idx = integer(lower_bound( Knots+_degree, Knots+n+1, x ) - Knots);
    if ( idx > integer(_degree) ) --idx;
    idx -= _degree;
    return BSplineEval<_degree>::eval_DDD(x,Knots+idx,yPoly+idx);
  }

  template <size_t _degree>
  real_type
  BSpline<_degree>::operator () ( real_type x ) const {
    if ( x >= X[npts-1] ) {
      real_type dx = x - X[npts-1];
      return s_R + dx * ( ds_R + 0.5 * dds_R * dx );
    } else if ( x <= X[0] ) {
      real_type dx = x - X[0];
      return s_L + dx * ( ds_L + 0.5 * dds_L * dx );
    } else {
      integer idx = knot_search( x );
      return BSplineEval<_degree>::eval(x,knots+idx,yPolygon+idx);
    }
  }

  template <size_t _degree>
  real_type
  BSpline<_degree>::D( real_type x ) const {
    if ( x >= X[npts-1] ) {
      real_type dx = x - X[npts-1];
      return ds_R + dds_R * dx;
    } else if ( x <= X[0] ) {
      real_type dx = x - X[0];
      return ds_L + dds_L * dx;
    } else {
      integer idx = knot_search( x );
      return BSplineEval<_degree>::eval_D(x,knots+idx,yPolygon+idx);
    }
  }

  template <size_t _degree>
  real_type
  BSpline<_degree>::DD( real_type x ) const {
    if ( x >= X[npts-1] ) {
      return dds_R;
    } else if ( x <= X[0] ) {
      return dds_L;
    } else {
      integer idx = knot_search( x );
      return BSplineEval<_degree>::eval_DD(x,knots+idx,yPolygon+idx);
    }
  }

  template <size_t _degree>
  real_type
  BSpline<_degree>::DDD( real_type x ) const {
    if ( x >= X[npts-1] || x <= X[0] ) {
      return 0;
    } else {
      integer idx = knot_search( x );
      return BSplineEval<_degree>::eval_DDD(x,knots+idx,yPolygon+idx);
    }
  }

  template <size_t _degree>
  integer // order
  BSpline<_degree>::coeffs( real_type cfs[],
                            real_type nodes[],
                            bool transpose ) const {
    SPLINE_ASSERT( false, "BSpline<_degree>::coeffs not yet implemented" );
#if 0
    size_t n = size_t(npts > 0 ? npts-1 : 0);
    for ( size_t i = 0; i < n; ++i ) {
      nodes[i] = X[i];
      real_type H  = X[i+1]-X[i];
      real_type DY = (Y[i+1]-Y[i])/H;
      real_type a = Y[i];
      real_type b = Yp[i];
      real_type c = (3*DY-2*Yp[i]-Yp[i+1])/H;
      real_type d = (Yp[i+1]+Yp[i]-2*DY)/(H*H);
      if ( transpose ) {
        cfs[4*i+3] = a;
        cfs[4*i+2] = b;
        cfs[4*i+1] = c;
        cfs[4*i+0] = d;
      } else {
        cfs[i+3*n] = a;
        cfs[i+2*n] = b;
        cfs[i+1*n] = c;
        cfs[i+0*n] = d;
      }
    }
#endif
    return _degree+1;
  }

  // Implementation
  template <size_t _degree>
  void
  BSpline<_degree>::copySpline( BSpline const & S ) {
    BSpline::reserve(S.npts);
    npts = S.npts;
    std::copy( S.X,        S.X+npts,               X );
    std::copy( S.Y,        S.Y+npts,               Y );
    std::copy( S.knots,    S.knots+npts+_degree+1, knots );
    std::copy( S.yPolygon, S.yPolygon+npts,        yPolygon );
  }

  //! change X-range of the spline
  template <size_t _degree>
  void
  BSpline<_degree>::setRange( real_type xmin, real_type xmax ) {
    Spline::setRange( xmin, xmax );
    real_type recS = ( X[npts-1] - X[0] ) / (xmax - xmin);
    real_type * iy = Y;
    while ( iy < Y + npts ) *iy++ *= recS;
  }

  template <size_t _degree>
  void
  BSpline<_degree>::writeToStream( ostream_type & s ) const {
    size_t nseg = size_t(npts > 0 ? npts - 1 : 0);
    for ( size_t i = 0; i < nseg; ++i )
      s << "segment N." << setw(4) << i
        << " X:[" << X[i] << ", " << X[i+1]
        << "] Y:[" << Y[i] << ", " << Y[i+1]
        << "] slope: " << (Y[i+1]-Y[i])/(X[i+1]-X[i])
        << '\n';
  }

  template <size_t _degree>
  void
  BSpline<_degree>::bases( real_type x, real_type vals[] ) const {
    std::fill( vals, vals+npts, 0 );
    if ( x >= X[0] && x <= X[npts-1] ) {
      integer idx = knot_search( x );
      BSplineBase<_degree>::eval( x, knots+idx, vals+idx );
    }
  }

  template <size_t _degree>
  void
  BSpline<_degree>::bases_D( real_type x, real_type vals[] ) const {
    std::fill( vals, vals+npts, 0 );
    if ( x >= X[0] && x <= X[npts-1] ) {
      integer idx = knot_search( x );
      BSplineBase<_degree>::eval_D( x, knots+idx, vals+idx );
    }
  }

  template <size_t _degree>
  void
  BSpline<_degree>::bases_DD( real_type x, real_type vals[] ) const {
    std::fill( vals, vals+npts, 0 );
    if ( x >= X[0] && x <= X[npts-1] ) {
      integer idx = knot_search( x );
      BSplineBase<_degree>::eval_DD( x, knots+idx, vals+idx );
    }
  }

  template <size_t _degree>
  void
  BSpline<_degree>::bases_DDD( real_type x, real_type vals[] ) const {
    std::fill( vals, vals+npts, 0 );
    if ( x >= X[0] && x <= X[npts-1] ) {
      integer idx = knot_search( x );
      BSplineBase<_degree>::eval_DDD( x, knots+idx, vals+idx );
    }
  }

  template <size_t _degree>
  integer
  BSpline<_degree>::bases_nz( real_type x, real_type vals[] ) const {
    integer idx = knot_search( x );
    if ( x >= X[0] && x <= X[npts-1] ) BSplineBase<_degree>::eval( x, knots+idx, vals );
    return idx;
  }

  template <size_t _degree>
  integer
  BSpline<_degree>::bases_D_nz( real_type x, real_type vals[] ) const {
    integer idx = knot_search( x );
    if ( x >= X[0] && x <= X[npts-1] )
      BSplineBase<_degree>::eval_D( x, knots+idx, vals );
    return idx;
  }

  template <size_t _degree>
  integer
  BSpline<_degree>::bases_DD_nz( real_type x, real_type vals[] ) const {
    integer idx = knot_search( x );
    if ( x >= X[0] && x <= X[npts-1] )
      BSplineBase<_degree>::eval_DD( x, knots+idx, vals );
    return idx;
  }

  template <size_t _degree>
  integer
  BSpline<_degree>::bases_DDD_nz( real_type x, real_type vals[] ) const {
    integer idx = knot_search( x );
    if ( x >= X[0] && x <= X[npts-1] )
      BSplineBase<_degree>::eval_DDD( x, knots+idx, vals );
    return idx;
  }

  #ifdef __clang__
    #pragma clang diagnostic ignored "-Wweak-template-vtables"
  #endif

  template class BSpline<1>;
  template class BSpline<2>;
  template class BSpline<3>;
  template class BSpline<4>;
  template class BSpline<5>;
  template class BSpline<6>;
  template class BSpline<7>;
  template class BSpline<8>;
  template class BSpline<9>;
  template class BSpline<10>;

}
