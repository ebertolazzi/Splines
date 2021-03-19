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

    Utils::Malloc<real_type> mem;

    real_type * m_DX;
    real_type * m_DXX;
    real_type * m_DY;
    real_type * m_DYY;
    real_type * m_DXY;
    real_type * m_DXYY;
    real_type * m_DXXY;
    real_type * m_DXXYY;
    void load( integer i, integer j, real_type bili5[6][6] ) const;

  public:

    //! spline constructor
    BiQuinticSplineBase( string const & name = "Spline" )
    : SplineSurf( name )
    , mem("BiQuinticSplineBase")
    , m_DX(nullptr)
    , m_DXX(nullptr)
    , m_DY(nullptr)
    , m_DYY(nullptr)
    , m_DXY(nullptr)
    , m_DXYY(nullptr)
    , m_DXXY(nullptr)
    {}

    virtual
    ~BiQuinticSplineBase() override
    { mem.free(); }

    real_type
    DxNode( integer i, integer j ) const
    { return m_DX[size_t(this->ipos_C(i,j))]; }

    real_type
    DyNode( integer i, integer j ) const
    { return m_DY[size_t(this->ipos_C(i,j))]; }

    real_type
    DxxNode( integer i, integer j ) const
    { return m_DXX[size_t(this->ipos_C(i,j))]; }

    real_type
    DyyNode( integer i, integer j ) const
    { return m_DYY[size_t(this->ipos_C(i,j))]; }

    real_type
    DxyNode( integer i, integer j ) const
    { return m_DXY[size_t(this->ipos_C(i,j))]; }

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
    Dxx( real_type x, real_type y ) const override;

    virtual
    real_type
    Dxy( real_type x, real_type y ) const override;

    virtual
    real_type
    Dyy( real_type x, real_type y ) const override;
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
    virtual void makeSpline() override;
  public:

    //! spline constructor
    BiQuinticSpline( string const & name = "Spline" )
    : BiQuinticSplineBase( name )
    {}

    virtual
    ~BiQuinticSpline() override
    {}

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

// EOF: SplineBiQuintic.hxx
