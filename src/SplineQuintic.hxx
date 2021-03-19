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

  typedef enum {
    CUBIC_QUINTIC = 0,
    PCHIP_QUINTIC,
    AKIMA_QUINTIC,
    BESSEL_QUINTIC
  } QUINTIC_SPLINE_TYPE;

  //! Quintic spline class
  class QuinticSpline : public QuinticSplineBase {
    QUINTIC_SPLINE_TYPE m_q_sub_type;
  public:

    using QuinticSplineBase::build;
    using QuinticSplineBase::reserve;

    //! spline constructor
    QuinticSpline( string const & name = "Spline" )
    : QuinticSplineBase( name )
    , m_q_sub_type(CUBIC_QUINTIC)
    {}

    //! spline destructor
    virtual
    ~QuinticSpline() override
    {}

    void
    setQuinticType( QUINTIC_SPLINE_TYPE qt )
    { m_q_sub_type = qt; }

    // --------------------------- VIRTUALS -----------------------------------
    //! Build a Monotone quintic spline from previously inserted points
    virtual
    void
    internal_build() override;

    virtual
    void
    setup( GenericContainer const & gc ) override;

  };

}

// EOF: SplineQuintic.hxx
