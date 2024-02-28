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
 |   ____  _ _ _                       ____        _ _
 |  | __ )(_) (_)_ __   ___  __ _ _ __/ ___| _ __ | (_)_ __   ___
 |  |  _ \| | | | '_ \ / _ \/ _` | '__\___ \| '_ \| | | '_ \ / _ \
 |  | |_) | | | | | | |  __/ (_| | |   ___) | |_) | | | | | |  __/
 |  |____/|_|_|_|_| |_|\___|\__,_|_|  |____/| .__/|_|_|_| |_|\___|
 |                                          |_|
\*/

namespace Splines {

  //! bilinear spline base class
  class BilinearSpline : public SplineSurf {

    void make_spline() override {}

    using SplineSurf::m_nx;
    using SplineSurf::m_ny;

    using SplineSurf::m_X;
    using SplineSurf::m_Y;
    using SplineSurf::m_Z;

  public:

    //! spline constructor
    BilinearSpline( string const & name = "Spline" )
    : SplineSurf(name)
    {}

    ~BilinearSpline() override {}

    real_type eval( real_type x, real_type y ) const override;

    void D( real_type x, real_type y, real_type d[3] ) const override;
    real_type Dx( real_type x, real_type y ) const override;
    real_type Dy( real_type x, real_type y ) const override;

    void DD( real_type x, real_type y, real_type dd[6] ) const override;
    real_type Dxx( real_type , real_type ) const override { return 0; }
    real_type Dxy( real_type , real_type ) const override { return 0; }
    real_type Dyy( real_type , real_type ) const override { return 0; }

    void write_to_stream( ostream_type & s ) const override;
    char const * type_name() const override;

  };

}

// EOF: SplineBilinear.hxx
