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
 |   ____  _  ___        _       _   _      ____        _ _            ____
 |  | __ )(_)/ _ \ _   _(_)_ __ | |_(_) ___/ ___| _ __ | (_)_ __   ___| __ )  __ _ ___  ___
 |  |  _ \| | | | | | | | | '_ \| __| |/ __\___ \| '_ \| | | '_ \ / _ \  _ \ / _` / __|/ _ \
 |  | |_) | | |_| | |_| | | | | | |_| | (__ ___) | |_) | | | | | |  __/ |_) | (_| \__ \  __/
 |  |____/|_|\__\_\\__,_|_|_| |_|\__|_|\___|____/| .__/|_|_|_| |_|\___|____/ \__,_|___/\___|
 |                                               |_|
\*/

namespace Splines {

  //! Bi-quintic spline base class
  class BiQuinticSplineBase : public SplineSurf {
  protected:

    Malloc_real mem;

    real_type * m_DX{nullptr};
    real_type * m_DXX{nullptr};
    real_type * m_DY{nullptr};
    real_type * m_DYY{nullptr};
    real_type * m_DXY{nullptr};
    real_type * m_DXYY{nullptr};
    real_type * m_DXXY{nullptr};
    real_type * m_DXXYY{nullptr};
    void load( integer i, integer j, real_type bili5[6][6] ) const;

  public:

    //! spline constructor
    BiQuinticSplineBase( string const & name = "Spline" )
    : SplineSurf( name )
    , mem("BiQuinticSplineBase")
    {}

    ~BiQuinticSplineBase() override
    { mem.free(); }

    real_type
    Dx_node( integer i, integer j ) const
    { return m_DX[size_t(this->ipos_C(i,j))]; }

    real_type
    Dy_node( integer i, integer j ) const
    { return m_DY[size_t(this->ipos_C(i,j))]; }

    real_type
    Dxx_node( integer i, integer j ) const
    { return m_DXX[size_t(this->ipos_C(i,j))]; }

    real_type
    Dyy_node( integer i, integer j ) const
    { return m_DYY[size_t(this->ipos_C(i,j))]; }

    real_type
    Dxy_node( integer i, integer j ) const
    { return m_DXY[size_t(this->ipos_C(i,j))]; }

    real_type eval( real_type x, real_type y ) const override;

    void D( real_type x, real_type y, real_type d[3] ) const override;
    real_type Dx( real_type x, real_type y ) const override;
    real_type Dy( real_type x, real_type y ) const override;

    void DD( real_type x, real_type y, real_type dd[6] ) const override;
    real_type Dxx( real_type x, real_type y ) const override;
    real_type Dxy( real_type x, real_type y ) const override;
    real_type Dyy( real_type x, real_type y ) const override;
  };

  /*\
   |   ____  _  ___        _       _   _      ____        _ _
   |  | __ )(_)/ _ \ _   _(_)_ __ | |_(_) ___/ ___| _ __ | (_)_ __   ___
   |  |  _ \| | | | | | | | | '_ \| __| |/ __\___ \| '_ \| | | '_ \ / _ \
   |  | |_) | | |_| | |_| | | | | | |_| | (__ ___) | |_) | | | | | |  __/
   |  |____/|_|\__\_\\__,_|_|_| |_|\__|_|\___|____/| .__/|_|_|_| |_|\___|
   |                                               |_|
  \*/
  //! cubic spline base class
  class BiQuinticSpline : public BiQuinticSplineBase {
    void make_spline() override;
  public:

    //! spline constructor
    BiQuinticSpline( string const & name = "BiQuinticSpline" )
    : BiQuinticSplineBase( name )
    {}

    ~BiQuinticSpline() override {}

    void write_to_stream( ostream_type & s ) const override;
    char const * type_name() const override;

  };

}

// EOF: SplineBiQuintic.hxx
