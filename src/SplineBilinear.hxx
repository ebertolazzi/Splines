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
    virtual void makeSpline() override {}
  public:

    using SplineSurf::m_nx;
    using SplineSurf::m_ny;

    using SplineSurf::m_X;
    using SplineSurf::m_Y;
    using SplineSurf::m_Z;

    //! spline constructor
    BilinearSpline( string const & name = "Spline" )
    : SplineSurf(name)
    {}

    virtual
    ~BilinearSpline() override
    {}

    //! Evaluate spline value
    virtual
    real_type
    operator () ( real_type x, real_type y ) const override;

    //! First derivative
    virtual
    void
    D( real_type x, real_type y, real_type d[3] ) const override;

    virtual
    real_type
    Dx( real_type x, real_type y ) const override;

    virtual
    real_type
    Dy( real_type x, real_type y ) const override;

    //! Second derivative
    virtual
    void
    DD( real_type x, real_type y, real_type dd[6] ) const override;

    virtual
    real_type
    Dxx( real_type , real_type ) const override
    { return 0; }

    virtual
    real_type
    Dxy( real_type , real_type ) const override
    { return 0; }

    virtual
    real_type
    Dyy( real_type , real_type ) const override
    { return 0; }

    //! Print spline coefficients
    virtual
    void
    writeToStream( ostream_type & s ) const override;

    //! Return spline typename
    virtual
    char const *
    type_name() const override;

  };

}

// EOF: SplineBilinear.hxx
