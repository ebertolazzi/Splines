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

/*\
 |    ___        _       _   _      ____        _ _            ____
 |   / _ \ _   _(_)_ __ | |_(_) ___/ ___| _ __ | (_)_ __   ___| __ )  __ _ ___  ___
 |  | | | | | | | | '_ \| __| |/ __\___ \| '_ \| | | '_ \ / _ \  _ \ / _` / __|/ _ \
 |  | |_| | |_| | | | | | |_| | (__ ___) | |_) | | | | | |  __/ |_) | (_| \__ \  __/
 |   \__\_\\__,_|_|_| |_|\__|_|\___|____/| .__/|_|_|_| |_|\___|____/ \__,_|___/\___|
 |                                       |_|
 |
\*/

namespace Splines {

  //! cubic quintic base class
  class QuinticSplineBase : public Spline {
  protected:
    Utils::Malloc<real_type> m_baseValue;

    real_type * m_Yp;
    real_type * m_Ypp;
    bool        m_external_alloc;

  public:

    using Spline::build;

    //! spline constructor
    QuinticSplineBase( string const & name = "Spline" )
    : Spline(name)
    , m_baseValue(name+"_memeory")
    , m_Yp(nullptr)
    , m_Ypp(nullptr)
    , m_external_alloc(false)
    {}

    virtual
    ~QuinticSplineBase() override
    {}

    void
    copySpline( QuinticSplineBase const & S );

    //! return the i-th node of the spline (y' component).
    real_type
    ypNode( integer i ) const
    { return m_Yp[size_t(i)]; }

    //! return the i-th node of the spline (y'' component).
    real_type
    yppNode( integer i ) const
    { return m_Ypp[size_t(i)]; }

    //! change X-range of the spline
    void
    setRange( real_type xmin, real_type xmax );

    //! Use externally allocated memory for `npts` points
    void
    reserve_external(
      integer       n,
      real_type * & p_x,
      real_type * & p_y,
      real_type * & p_Yp,
      real_type * & p_Ypp
    );

    // --------------------------- VIRTUALS -----------------------------------

    //! Evaluate spline value
    virtual
    real_type
    operator () ( real_type x ) const override;

    //! First derivative
    virtual
    real_type
    D( real_type x ) const override;

    //! Second derivative
    virtual
    real_type
    DD( real_type x ) const override;

    //! Third derivative
    virtual
    real_type
    DDD( real_type x ) const override;

    //! Fourth derivative
    virtual
    real_type
    DDDD( real_type x ) const override;

    //! Fifth derivative
    virtual
    real_type
    DDDDD( real_type x ) const override;

    //! Evaluate spline value knowing interval
    virtual
    real_type
    id_eval( integer ni, real_type x ) const override;

    //! First derivative
    virtual
    real_type
    id_D( integer ni, real_type x ) const override;

    //! Second derivative
    virtual
    real_type
    id_DD( integer ni, real_type x ) const override;

    //! Third derivative
    virtual
    real_type
    id_DDD( integer ni, real_type x ) const override;

    //! Fourth derivative
    virtual
    real_type
    id_DDDD( integer ni, real_type x ) const override;

    //! Fifth derivative
    virtual
    real_type
    id_DDDDD( integer ni, real_type x ) const override;

    //! Print spline coefficients
    virtual
    void
    writeToStream( ostream_type & s ) const override;

    //! Return spline type (as number)
    virtual
    unsigned
    type() const override
    { return QUINTIC_TYPE; }

    //! Allocate memory for `npts` points
    virtual
    void
    reserve( integer npts ) override;

    //! Cancel the support points, empty the spline.
    virtual
    void
    clear() override;

    //! get the piecewise polinomials of the spline
    virtual
    integer // order
    coeffs(
      real_type cfs[],
      real_type nodes[],
      bool      transpose = false
    ) const override;

    virtual
    integer // order
    order() const override;

  };

}

// EOF: SplineQuinticBase.hxx
