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

    Malloc_real m_mem_biquintic{"BiQuinticSplineBase"};

    real_type * m_DX{nullptr};
    real_type * m_DY{nullptr};

    real_type * m_DXX{nullptr};
    real_type * m_DYY{nullptr};
    real_type * m_DXY{nullptr};

    real_type * m_DXYY{nullptr};
    real_type * m_DXXY{nullptr};
    real_type * m_DXXYY{nullptr};

    using SplineSurf::m_nx;
    using SplineSurf::m_ny;

    using SplineSurf::m_X;
    using SplineSurf::m_Y;
    using SplineSurf::m_Z;

    real_type & Dx_node_ref    ( integer const i, integer const j ) { return m_DX    [ ipos_C(i,j) ]; }
    real_type & Dy_node_ref    ( integer const i, integer const j ) { return m_DY    [ ipos_C(i,j) ]; }

    real_type & Dxx_node_ref   ( integer const i, integer const j ) { return m_DXX   [ ipos_C(i,j) ]; }
    real_type & Dyy_node_ref   ( integer const i, integer const j ) { return m_DYY   [ ipos_C(i,j) ]; }
    real_type & Dxy_node_ref   ( integer const i, integer const j ) { return m_DXY   [ ipos_C(i,j) ]; }

    real_type & Dxyy_node_ref  ( integer const i, integer const j ) { return m_DXYY  [ ipos_C(i,j) ]; }
    real_type & Dxxy_node_ref  ( integer const i, integer const j ) { return m_DXXY  [ ipos_C(i,j) ]; }
    real_type & Dxxyy_node_ref ( integer const i, integer const j ) { return m_DXXYY [ ipos_C(i,j) ]; }

    void load( integer const i, integer const j, real_type bili5[6][6] ) const;

    #endif

  public:

    using SplineSurf::eval;

    //! spline constructor
    explicit
    BiQuinticSplineBase( string_view name = "Spline" )
    : SplineSurf( name )
    {}

    ~BiQuinticSplineBase() override
    { m_mem_biquintic.free(); }

    //!
    //! \name Estimated derivatives at interpolation nodes
    //!
    ///@{

    //!
    //! Estimated `x` derivatives at node `(i,j)`
    //!
    real_type
    Dx_node( integer const i, integer const j ) const
    { return m_DX[ ipos_C(i,j) ]; }

    //!
    //! Estimated `y` derivatives at node `(i,j)`
    //!
    real_type
    Dy_node( integer const i, integer const j ) const
    { return m_DY[ ipos_C(i,j) ]; }

    //!
    //! Estimated second `x` second derivatives at node `(i,j)`
    //!
    real_type
    Dxx_node( integer const i, integer const j ) const
    { return m_DXX[ ipos_C(i,j) ]; }

    //!
    //! Estimated second`y` derivatives at node `(i,j)`
    //!
    real_type
    Dyy_node( integer const i, integer const j ) const
    { return m_DYY[ ipos_C(i,j) ]; }

    //!
    //! Estimated mixed `xy` derivatives at node `(i,j)`
    //!
    real_type
    Dxy_node( integer const i, integer const j ) const
    { return m_DXY[ ipos_C(i,j) ]; }

    //!
    //! Estimated `xxy` second derivatives at node `(i,j)`
    //!
    real_type
    Dxxy_node( integer const i, integer const j ) const
    { return m_DXXY[ ipos_C(i,j) ]; }

    //!
    //! Estimated `xyy` derivatives at node `(i,j)`
    //!
    real_type
    Dxyy_node( integer const i, integer const j ) const
    { return m_DXYY[ ipos_C(i,j) ]; }

    //!
    //! Estimated mixed `xxyy` derivatives at node `(i,j)`
    //!
    real_type
    Dxxyy_node( integer const i, integer const j ) const
    { return m_DXXYY[ ipos_C(i,j) ]; }

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
    void D( real_type const x, real_type const y, real_type d[3] ) const override;
    //!
    //! Evaluate spline `x`  derivative at point \f$ (x,y) \f$
    //!
    real_type Dx( real_type const x, real_type const y ) const override;
    //!
    //! Evaluate spline `y`  derivative at point \f$ (x,y) \f$
    //!
    real_type Dy( real_type const x, real_type const y ) const override;

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
    void DD( real_type const x, real_type const y, real_type dd[6] ) const override;
    //!
    //! Evaluate spline `x` second derivative at point \f$ (x,y) \f$
    //!
    real_type Dxx( real_type const x, real_type const y ) const override;
    //!
    //! Evaluate spline `xy` mixed derivative at point \f$ (x,y) \f$
    //!
    real_type Dxy( real_type const x, real_type const y ) const override;
    //!
    //! Evaluate spline `y` second derivative at point \f$ (x,y) \f$
    //!
    real_type Dyy( real_type const x, real_type const y ) const override;

    ///@}

    #ifdef AUTODIFF_SUPPORT
    //!
    //! \name Autodiff
    //!
    ///@{
    autodiff::dual1st eval( autodiff::dual1st const & x, autodiff::dual1st const & y ) const;
    autodiff::dual2nd eval( autodiff::dual2nd const & x, autodiff::dual2nd const & y ) const;

    template <typename T1, typename T2>
    autodiff::HigherOrderDual<autodiff::detail::DualOrder<T1,T2>::value,real_type>
    eval( T1 const & x, T2 const & y ) const {
      autodiff::HigherOrderDual<autodiff::detail::DualOrder<T1,T2>::value,real_type> X{x}, Y{y};
      return eval( X, Y );
    }
    ///@}
    #endif

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
