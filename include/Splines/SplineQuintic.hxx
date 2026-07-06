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
 |      Università degli Studi di Trento                                    |
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

namespace Splines
{

  //! Quintic spline class
  class QuinticSpline final : public QuinticSplineBase
  {
    using QuinticSplineBase::m_sub_type;

  public:
    //!
    //! \name Constructors
    //!
    ///@{

    using QuinticSplineBase::build;
    using QuinticSplineBase::reserve;

    //!
    //! Build an empty spline of `QuinticSpline` type
    //!
    //! \param type spline type
    //! \param name the name of the spline
    //!
    explicit QuinticSpline( Spline_sub_type type = Spline_sub_type::PCHIP, string_view name = "Spline" )
      : QuinticSplineBase( type, name )
    {
    }

    //!
    //! Build an empty spline of `QuinticSpline` type
    //!
    //! \param name the name of the spline
    //!
    explicit QuinticSpline( string_view name ) : QuinticSplineBase( Spline_sub_type::PCHIP, name ) {}

    //!
    //! Spline destructor.
    //!
    ~QuinticSpline() override = default;

    ///@}

    // --------------------------- VIRTUALS -----------------------------------
    //! Build a Monotone quintic spline from previously inserted points
    void build() override;

    //! Build a Monotone quintic spline from data from `gc`
    void setup( GenericContainer const & gc ) override;
  };

}  // namespace Splines

// EOF: SplineQuintic.hxx
