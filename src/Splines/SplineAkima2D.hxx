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
 |      _    _    _                 ____  ____            _ _
 |     / \  | | _(_)_ __ ___   __ _|___ \|  _ \ ___ _ __ | (_)_ __   ___
 |    / _ \ | |/ / | '_ ` _ \ / _` | __) | | | / __| '_ \| | | '_ \ / _ \
 |   / ___ \|   <| | | | | | | (_| |/ __/| |_| \__ \ |_) | | | | | |  __/
 |  /_/   \_\_|\_\_|_| |_| |_|\__,_|_____|____/|___/ .__/|_|_|_| |_|\___|
 |                                                 |_|
\*/

namespace Splines {

  //!
  //! Smooth Curve Fitting Based on Local Procedures
  //!
  //! *Reference*
  //!
  //! - *Hiroshi Akima*, A Method of Bivariate Interpolation and
  //!   Smooth Surface Fitting for Irregularly Distributed Data Points.
  //!   ACM Transactions on Mathematical Software, Vol.4, 148-164, 1978.
  //!
  class Akima2Dspline : public BiCubicSplineBase {
    void make_spline() override;

  public:

    //!
    //! Build an empty spline of `Akima2Dspline` type
    //!
    //! \param name the name of the spline
    //!
    Akima2Dspline( string const & name = "Spline" )
    : BiCubicSplineBase( name )
    {}

    //!
    //! Spline destructor.
    //!
    ~Akima2Dspline() override {}

    void write_to_stream( ostream_type & s ) const override;
    char const * type_name() const override;

  };
}

// EOF: SplineAkima2D.hxx

