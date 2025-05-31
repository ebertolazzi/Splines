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

    void make_spline() override { m_search_x.must_reset(); m_search_y.must_reset(); }

    using SplineSurf::m_nx;
    using SplineSurf::m_ny;

    using SplineSurf::m_X;
    using SplineSurf::m_Y;
    using SplineSurf::m_Z;

  public:

    using SplineSurf::eval;

    //!
    //! Build an empty spline of `BilinearSpline` type
    //!
    //! \param name the name of the spline
    //!
    explicit
    BilinearSpline( string_view name = "BilinearSpline" )
    : SplineSurf(name)
    {}

    //!
    //! Spline destructor.
    //!
    ~BilinearSpline() override {}

    real_type eval( real_type const x, real_type const y ) const override;

    void D( real_type const x, real_type const y, real_type d[3] ) const override;
    real_type Dx( real_type const x, real_type const y ) const override;
    real_type Dy( real_type const x, real_type const y ) const override;

    void DD( real_type const x, real_type const y, real_type dd[6] ) const override;
    real_type Dxx( real_type const, real_type const ) const override { return 0; }
    real_type Dxy( real_type const, real_type const ) const override { return 0; }
    real_type Dyy( real_type const, real_type const ) const override { return 0; }

    #ifdef AUTIDIFF_SUPPORT
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

    void write_to_stream( ostream_type & s ) const override;
    char const * type_name() const override;

  };

}

// EOF: SplineBilinear.hxx
