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

    #ifndef DOXYGEN_SHOULD_SKIP_THIS

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

    real_type & Dx_node_ref  ( integer i, integer j ) { return m_DX  [ this->ipos_C(i,j) ]; }
    real_type & Dy_node_ref  ( integer i, integer j ) { return m_DY  [ this->ipos_C(i,j) ]; }
    real_type & Dxy_node_ref ( integer i, integer j ) { return m_DXY [ this->ipos_C(i,j) ]; }

    #endif

  public:

    //! spline constructor
    explicit
    BiCubicSplineBase( string_view name = "BiCubicSplineBase" );

    ~BiCubicSplineBase() override {}

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
    //! Build an empty spline of `BiCubicSpline` type
    //!
    //! \param name the name of the spline
    //!
    explicit
    BiCubicSpline( string_view name = "BiCubicSpline" )
    : BiCubicSplineBase( name )
    {}

    //!
    //! Spline destructor.
    //!
    ~BiCubicSpline() override {}

    void write_to_stream( ostream_type & s ) const override;
    char const * type_name() const override;

  };

}

// EOF: SplineBiCubic.hxx
