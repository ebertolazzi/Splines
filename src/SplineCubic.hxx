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
 |    ____      _     _      ____        _ _
 |   / ___|   _| |__ (_) ___/ ___| _ __ | (_)_ __   ___
 |  | |  | | | | '_ \| |/ __\___ \| '_ \| | | '_ \ / _ \
 |  | |__| |_| | |_) | | (__ ___) | |_) | | | | | |  __/
 |   \____\__,_|_.__/|_|\___|____/| .__/|_|_|_| |_|\___|
 |                                |_|
\*/

namespace Splines {

  typedef enum {
    EXTRAPOLATE_BC = 0,
    NATURAL_BC,
    PARABOLIC_RUNOUT_BC,
    NOT_A_KNOT
  } CUBIC_SPLINE_TYPE_BC;

  #ifndef DOXYGEN_SHOULD_SKIP_THIS

  void
  CubicSpline_build(
    real_type const * X,
    real_type const * Y,
    real_type       * Yp,
    integer           npts,
    CUBIC_SPLINE_TYPE_BC bc0,
    CUBIC_SPLINE_TYPE_BC bcn
  );

  void
  CubicSpline_build(
    real_type const * X,
    real_type const * Y,
    real_type       * Yp,
    real_type       * Ypp,
    real_type       * L,
    real_type       * D,
    real_type       * U,
    integer           npts,
    CUBIC_SPLINE_TYPE_BC bc0,
    CUBIC_SPLINE_TYPE_BC bcn
  );

  #endif

  //! Cubic Spline Management Class
  class CubicSpline : public CubicSplineBase {
  private:
    CUBIC_SPLINE_TYPE_BC m_bc0, m_bcn;
  public:

    //! \name Constructor
    ///@{

    #ifndef DOXYGEN_SHOULD_SKIP_THIS
    using CubicSplineBase::build;
    using CubicSplineBase::reserve;
    #endif

    //!
    //! spline constructor
    //!
    CubicSpline( string const & name = "CubicSpline" )
    : CubicSplineBase( name )
    , m_bc0( EXTRAPOLATE_BC )
    , m_bcn( EXTRAPOLATE_BC )
    {}

    //!
    //! spline destructor
    //!
    ~CubicSpline() override {}

    ///@}

    //! \name Setup
    ///@{

    //!
    //! Set the boudary consition for initial point
    //! \param[in] bc0  initial boundary condition.
    //!
    void
    setInitialBC( CUBIC_SPLINE_TYPE_BC bc0 )
    { m_bc0 = bc0; }

    //!
    //! Set the boudary consition for final point
    //! \param[in] bcn final boundary condition.
    //!
    void
    setFinalBC( CUBIC_SPLINE_TYPE_BC bcn )
    { m_bcn = bcn; }

    // --------------------------- VIRTUALS -----------------------------------

    void build() override;
    void setup( GenericContainer const & gc ) override;

    ///@}

    //!
    //! Return spline type (as number)
    //!
    unsigned type() const override { return CUBIC_TYPE; }
  };

}

// EOF: SplineCubic.hxx
