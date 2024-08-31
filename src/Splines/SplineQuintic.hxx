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
 |    ___        _       _   _      ____        _ _
 |   / _ \ _   _(_)_ __ | |_(_) ___/ ___| _ __ | (_)_ __   ___
 |  | | | | | | | | '_ \| __| |/ __\___ \| '_ \| | | '_ \ / _ \
 |  | |_| | |_| | | | | | |_| | (__ ___) | |_) | | | | | |  __/
 |   \__\_\\__,_|_|_| |_|\__|_|\___|____/| .__/|_|_|_| |_|\___|
 |                                       |_|
 |
\*/

namespace Splines {

  #ifndef DOXYGEN_SHOULD_SKIP_THIS

  using QuinticSpline_sub_type = enum class QuinticSpline_sub_type : integer {
    CUBIC  = 0,
    PCHIP  = 1,
    AKIMA  = 2,
    BESSEL = 3
  };

  #endif

  //! Quintic spline class
  class QuinticSpline : public QuinticSplineBase {
    QuinticSpline_sub_type m_q_sub_type{QuinticSpline_sub_type::CUBIC};
  public:

    //!
    //! \name Constructors
    //!
    ///@{

    #ifndef DOXYGEN_SHOULD_SKIP_THIS
    using QuinticSplineBase::build;
    using QuinticSplineBase::reserve;
    #endif

    //!
    //! Build an empty spline of `QuinticSpline` type
    //!
    //! \param name the name of the spline
    //!
    QuinticSpline( string const & name = "Spline" )
    : QuinticSplineBase( name )
    {}

    //!
    //! Spline destructor.
    //!
    ~QuinticSpline() override {}

    ///@}

    //!
    //! Set spline type
    //!
    //! - CUBIC
    //! - PCHIP
    //! - AKIMA
    //! - BESSEL
    //!
    void
    setQuinticType( QuinticSpline_sub_type qt )
    { m_q_sub_type = qt; }

    // --------------------------- VIRTUALS -----------------------------------
    //! Build a Monotone quintic spline from previously inserted points
    void build() override;
    //! Build a Monotone quintic spline from data from `gc`
    void setup( GenericContainer const & gc ) override;
  };

}

// EOF: SplineQuintic.hxx
