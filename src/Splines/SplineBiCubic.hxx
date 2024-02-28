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

namespace Splines {

  /*\
   |   ____  _  ____      _     _      ____        _ _            ____
   |  | __ )(_)/ ___|   _| |__ (_) ___/ ___| _ __ | (_)_ __   ___| __ )  __ _ ___  ___
   |  |  _ \| | |  | | | | '_ \| |/ __\___ \| '_ \| | | '_ \ / _ \  _ \ / _` / __|/ _ \
   |  | |_) | | |__| |_| | |_) | | (__ ___) | |_) | | | | | |  __/ |_) | (_| \__ \  __/
   |  |____/|_|\____\__,_|_.__/|_|\___|____/| .__/|_|_|_| |_|\___|____/ \__,_|___/\___|
   |                                        |_|
  \*/

  //!
  //! Bi-cubic spline base class
  //!
  class BiCubicSplineBase : public SplineSurf {
  protected:

    Malloc_real m_mem_bicubic;

    real_type * m_DX{nullptr};
    real_type * m_DY{nullptr};
    real_type * m_DXY{nullptr};

    using SplineSurf::m_nx;
    using SplineSurf::m_ny;

    using SplineSurf::m_X;
    using SplineSurf::m_Y;
    using SplineSurf::m_Z;

    void load( integer i, integer j, real_type bili3[4][4] ) const;

  public:

    //! spline constructor
    BiCubicSplineBase( string const & name = "Spline" )
    : SplineSurf( name )
    , m_mem_bicubic("BiCubicSplineBase")
    {}

    ~BiCubicSplineBase() override {}

    //!
    //! \name Estimated derivatives at interpolation nodes
    //!
    ///@{

    real_type
    Dx_node( integer i, integer j ) const
    { return m_DX[size_t(this->ipos_C(i,j))]; }

    real_type
    Dy_node( integer i, integer j ) const
    { return m_DY[size_t(this->ipos_C(i,j))]; }

    real_type
    Dxy_node( integer i, integer j ) const
    { return m_DXY[size_t(this->ipos_C(i,j))]; }

    ///@}

    //!
    //! \name Evaluate
    //!
    ///@{
    real_type eval( real_type x, real_type y ) const override;

    void D( real_type x, real_type y, real_type d[3] ) const override;
    real_type Dx( real_type x, real_type y ) const override;
    real_type Dy( real_type x, real_type y ) const override;

    void DD( real_type x, real_type y, real_type dd[6] ) const override;
    real_type Dxx( real_type x, real_type y ) const override;
    real_type Dxy( real_type x, real_type y ) const override;
    real_type Dyy( real_type x, real_type y ) const override;
    ///@}
  };

  /*\
   |   ____  _  ____      _     _      ____        _ _
   |  | __ )(_)/ ___|   _| |__ (_) ___/ ___| _ __ | (_)_ __   ___
   |  |  _ \| | |  | | | | '_ \| |/ __\___ \| '_ \| | | '_ \ / _ \
   |  | |_) | | |__| |_| | |_) | | (__ ___) | |_) | | | | | |  __/
   |  |____/|_|\____\__,_|_.__/|_|\___|____/| .__/|_|_|_| |_|\___|
   |                                        |_|
  \*/
  //!
  //! Cubic spline base class
  //!
  class BiCubicSpline : public BiCubicSplineBase {
    void make_spline() override;

    using BiCubicSplineBase::m_mem_bicubic;
    using BiCubicSplineBase::m_DX;
    using BiCubicSplineBase::m_DY;
    using BiCubicSplineBase::m_DXY;

  public:

    //!
    //! Spline constructor
    //!
    BiCubicSpline( string const & name = "Spline" )
    : BiCubicSplineBase( name )
    {}

    ~BiCubicSpline() override {}

    void write_to_stream( ostream_type & s ) const override;
    char const * type_name() const override;

  };

}

// EOF: SplineBiCubic.hxx
