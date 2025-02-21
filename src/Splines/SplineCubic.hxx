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

  #ifndef DOXYGEN_SHOULD_SKIP_THIS

  using CubicSpline_BC = enum class CubicSpline_BC : integer {
    EXTRAPOLATE      = 0,
    NATURAL          = 1,
    PARABOLIC_RUNOUT = 2,
    NOT_A_KNOT       = 3
  };

  void
  CubicSpline_build(
    real_type const X[],
    real_type const Y[],
    real_type       Yp[],
    integer         npts,
    CubicSpline_BC  bc0,
    CubicSpline_BC  bcn
  );

  void
  CubicSpline_build(
    real_type const X[],
    real_type const Y[],
    real_type       Yp[],
    real_type       Ypp[],
    real_type       L[],
    real_type       D[],
    real_type       U[],
    integer         npts,
    CubicSpline_BC  bc0,
    CubicSpline_BC  bcn
  );

  #endif

  //!
  //! Cubic Spline Management Class
  //!
  class CubicSpline : public CubicSplineBase {
  private:
    CubicSpline_BC m_bc0{CubicSpline_BC::EXTRAPOLATE};
    CubicSpline_BC m_bcn{CubicSpline_BC::EXTRAPOLATE};
  public:
    //!
    //! \name Constructors
    //!
    ///@{

    #ifndef DOXYGEN_SHOULD_SKIP_THIS
    using CubicSplineBase::build;
    using CubicSplineBase::reserve;
    #endif

    //!
    //! Build an empty spline of `CubicSpline` type
    //!
    //! \param name the name of the spline
    //!
    explicit
    CubicSpline( string_view name = "CubicSpline" )
    : CubicSplineBase( name )
    {}

    //!
    //! Spline destructor.
    //!
    ~CubicSpline() override {}

    ///@}

    //!
    //! \name Setup
    //!
    ///@{

    //!
    //! Set the boudary consition for initial point
    //! \param[in] bc0  initial boundary condition.
    //!
    void
    set_initial_BC( CubicSpline_BC bc0 )
    { m_bc0 = bc0; }

    //!
    //! Set the boudary consition for final point
    //! \param[in] bcn final boundary condition.
    //!
    void
    set_final_BC( CubicSpline_BC bcn )
    { m_bcn = bcn; }

    // --------------------------- VIRTUALS -----------------------------------

    void build() override;
    void setup( GenericContainer const & gc ) override;

    ///@}

    //!
    //! Return spline type (as number)
    //!
    SplineType1D type() const override { return SplineType1D::CUBIC; }

    #ifdef SPLINES_BACK_COMPATIBILITY
    void setInitialBC( CubicSpline_BC bc0 ) { m_bc0 = bc0; }
    void setFinalBC( CubicSpline_BC bcn ) { m_bcn = bcn; }
    #endif

  };

}

// EOF: SplineCubic.hxx
