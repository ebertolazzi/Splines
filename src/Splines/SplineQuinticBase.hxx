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

  //!
  //! Quintic spline base class
  //!
  class QuinticSplineBase : public Spline {
  protected:
    Malloc_real m_baseValue;
    real_type * m_Yp;
    real_type * m_Ypp;
    bool        m_external_alloc;

  public:

    //!
    //! \name Constructors
    //!
    ///@{

    #ifndef DOXYGEN_SHOULD_SKIP_THIS
    using Spline::build;
    #endif

    //!
    //! Spline constructor
    //!
    QuinticSplineBase( string const & name = "Spline" )
    : Spline(name)
    , m_baseValue(name+"_memeory")
    , m_Yp(nullptr)
    , m_Ypp(nullptr)
    , m_external_alloc(false)
    {}

    ~QuinticSplineBase() override {}

    ///@}

    void copy_spline( QuinticSplineBase const & S );
    #ifndef SPLINES_NO_COMPATIBILITY
    void copySpline( QuinticSplineBase const & S ) { this->copy_spline(S); }
    #endif

    //!
    //! \name Info
    //!
    ///@{

    //!
    //! Return the i-th node of the spline (y' component).
    //!
    real_type yp_node( integer i ) const { return m_Yp[size_t(i)]; }
    #ifndef SPLINES_NO_COMPATIBILITY
    real_type ypNode( integer i ) const { return this->yp_node(i); }
    #endif

    //!
    //! Return the i-th node of the spline (y'' component).
    //!
    real_type ypp_node( integer i ) const { return m_Ypp[size_t(i)]; }
    #ifndef SPLINES_NO_COMPATIBILITY
    real_type yppNode( integer i ) const { return this->ypp_node(i); }
    #endif

    void
    y_min_max(
      integer   & i_min_pos,
      real_type & x_min_pos,
      real_type & y_min,
      integer   & i_max_pos,
      real_type & x_max_pos,
      real_type & y_max
    ) const override;

    void
    y_min_max(
      vector<integer>   & i_min_pos,
      vector<real_type> & x_min_pos,
      vector<real_type> & y_min,
      vector<integer>   & i_max_pos,
      vector<real_type> & x_max_pos,
      vector<real_type> & y_max
    ) const override;

    void write_to_stream( ostream_type & s ) const override;

    SplineType1D type() const override { return SplineType1D::QUINTIC; }

    ///@}

    //!
    //! Change X-range of the spline
    //!
    void set_range( real_type xmin, real_type xmax );
    #ifndef SPLINES_NO_COMPATIBILITY
    void setRange( real_type xmin, real_type xmax ) { this->set_range( xmin, xmax ); }
    #endif

    //!
    //! Use externally allocated memory for `npts` points
    //!
    void
    reserve_external(
      integer       n,
      real_type * & p_x,
      real_type * & p_y,
      real_type * & p_Yp,
      real_type * & p_Ypp
    );

    // --------------------------- VIRTUALS -----------------------------------

    //!
    //! \name Evaluation Aliases
    //!
    ///@{
    real_type operator () ( real_type x ) const override;
    real_type D( real_type x ) const override;
    real_type DD( real_type x ) const override;
    real_type DDD( real_type x ) const override;
    real_type DDDD( real_type x ) const override;
    real_type DDDDD( real_type x ) const override;
    ///@}

    //!
    //! \name Evaluation when segment is known
    ///@{
    real_type id_eval( integer ni, real_type x ) const override;
    real_type id_D( integer ni, real_type x ) const override;
    real_type id_DD( integer ni, real_type x ) const override;
    real_type id_DDD( integer ni, real_type x ) const override;
    real_type id_DDDD( integer ni, real_type x ) const override;
    real_type id_DDDDD( integer ni, real_type x ) const override;
    ///@}

    void reserve( integer npts ) override;
    void clear() override;

    //!
    //! Get the piecewise polinomials of the spline
    //!
    integer // order
    coeffs(
      real_type * const cfs,
      real_type * const nodes,
      bool              transpose = false
    ) const override;

    integer order() const override;

  };

}

// EOF: SplineQuinticBase.hxx
