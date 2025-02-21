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
    #ifndef DOXYGEN_SHOULD_SKIP_THIS

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

    #endif

  public:

    //! spline constructor
    explicit
    BiQuinticSplineBase( string_view name = "Spline" )
    : SplineSurf( name )
    , mem("BiQuinticSplineBase")
    {}

    ~BiQuinticSplineBase() override
    { mem.free(); }

    //!
    //! \name Estimated derivatives at interpolation nodes
    //!
    ///@{

    //!
    //! Estimated `x` derivatives at node `(i,j)`
    //!
    real_type
    Dx_node( integer i, integer j ) const
    { return m_DX[ this->ipos_C(i,j) ]; }

    //!
    //! Estimated `y` derivatives at node `(i,j)`
    //!
    real_type
    Dy_node( integer i, integer j ) const
    { return m_DY[ this->ipos_C(i,j) ]; }

    //!
    //! Estimated `x` second derivatives at node `(i,j)`
    //!
    real_type
    Dxx_node( integer i, integer j ) const
    { return m_DXX[ this->ipos_C(i,j) ]; }

    //!
    //! Estimated `y` derivatives at node `(i,j)`
    //!
    real_type
    Dyy_node( integer i, integer j ) const
    { return m_DYY[ this->ipos_C(i,j) ]; }

    //!
    //! Estimated mixed `xy` derivatives at node `(i,j)`
    //!
    real_type
    Dxy_node( integer i, integer j ) const
    { return m_DXY[ this->ipos_C(i,j) ]; }
    ///@}

    //!
    //! \name Evaluate
    //!
    ///@{
    //!
    //! Evaluate spline at point \f$ (x,y) \f$
    //!
    real_type eval( real_type x, real_type y ) const override;

    //!
    //! Evaluate spline with derivative at point \f$ (x,y) \f$
    //!
    //! - `d[0]` the value of the spline
    //! - `d[1]` the value of the spline `x` derivative
    //! - `d[2]` the value of the spline `y` derivative
    //!
    void D( real_type x, real_type y, real_type d[3] ) const override;
    //!
    //! Evaluate spline `x`  derivative at point \f$ (x,y) \f$
    //!
    real_type Dx( real_type x, real_type y ) const override;
    //!
    //! Evaluate spline `y`  derivative at point \f$ (x,y) \f$
    //!
    real_type Dy( real_type x, real_type y ) const override;

    //!
    //! Evaluate spline with derivative at point \f$ (x,y) \f$
    //!
    //! - `d[0]` the value of the spline
    //! - `d[1]` the value of the spline `x` derivative
    //! - `d[2]` the value of the spline `y` derivative
    //! - `d[3]` the value of the spline `x` second derivative
    //! - `d[4]` the value of the spline `y` second derivative
    //! - `d[5]` the value of the spline `xy` mixed derivative
    //!
    void DD( real_type x, real_type y, real_type dd[6] ) const override;
    //!
    //! Evaluate spline `x` second derivative at point \f$ (x,y) \f$
    //!
    real_type Dxx( real_type x, real_type y ) const override;
    //!
    //! Evaluate spline `xy` mixed derivative at point \f$ (x,y) \f$
    //!
    real_type Dxy( real_type x, real_type y ) const override;
    //!
    //! Evaluate spline `y` second derivative at point \f$ (x,y) \f$
    //!
    real_type Dyy( real_type x, real_type y ) const override;

    ///@}
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

    //!
    //! Build an empty spline of `BiQuinticSpline` type
    //!
    //! \param name the name of the spline
    //!
    explicit
    BiQuinticSpline( string_view name = "BiQuinticSpline" )
    : BiQuinticSplineBase( name )
    {}

    //!
    //! Spline destructor.
    //!
    ~BiQuinticSpline() override {}

    void write_to_stream( ostream_type & s ) const override;
    char const * type_name() const override;

  };

}

// EOF: SplineBiQuintic.hxx
