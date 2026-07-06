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

#pragma once

#ifndef SPLINE_SURF_HH
#define SPLINE_SURF_HH

namespace Splines
{

  /*\
   |   ____        _ _            ____              __
   |  / ___| _ __ | (_)_ __   ___/ ___| _   _ _ __ / _|
   |  \___ \| '_ \| | | '_ \ / _ \___ \| | | | '__| |_
   |   ___) | |_) | | | | | |  __/___) | |_| | |  |  _|
   |  |____/| .__/|_|_|_| |_|\___|____/ \__,_|_|  |_|
   |        |_|
  \*/

  //!
  //! Spline Management Class
  //!
  class SplineSurf
  {
  protected:
    string const m_name;
    bool         m_x_closed     = false;
    bool         m_y_closed     = false;
    bool         m_x_can_extend = true;
    bool         m_y_can_extend = true;

    integer m_nx = 0;
    integer m_ny = 0;

    real_type * m_X_ptr = nullptr;
    real_type * m_Y_ptr = nullptr;
    Vec         mX;
    Vec         mY;
    MatC        mZ;

    real_type m_Z_min = 0;
    real_type m_Z_max = 0;

    Utils::SearchInterval<real_type, integer> m_search_x;
    Utils::SearchInterval<real_type, integer> m_search_y;

    //!
    //! Find interval and calculate local delta/intervals
    //! Returns: { i, j, dx, dy, DX, DY }
    //!
    std::tuple<integer, integer, real_type, real_type, real_type, real_type> find_patch(
      real_type const x,
      real_type const y ) const
    {
      std::pair<integer, real_type> X( 0, x ), Y( 0, y );
      m_search_x.find( X );
      m_search_y.find( Y );

      integer const i = X.first;
      integer const j = Y.first;

      real_type const dx = X.second - mX.coeff( i );
      real_type const dy = Y.second - mY.coeff( j );

      // Assumo mX e mY siano i vettori Eigen corretti (ho uniformato m_X -> mX)
      real_type const DX = mX.coeff( i + 1 ) - mX.coeff( i );
      real_type const DY = mY.coeff( j + 1 ) - mY.coeff( j );

      return { i, j, dx, dy, DX, DY };
    }

    real_type & z_node_ref( integer const i, integer const j ) { return mZ.coeffRef( i, j ); }

    template <typename Derived> void load_Z( const Eigen::ArrayBase<Derived> & Z, bool transposed )
    {
      if ( transposed )
      {
        UTILS_ASSERT(
          Z.rows() >= m_ny && Z.cols() >= m_nx,
          "SplineSurf::load_Z( Z, transposed={} ) bad dimension found {} x {} expected {} x {}",
          transposed,
          Z.rows(),
          Z.cols(),
          m_ny,
          m_nx );
        for ( integer ix = 0; ix < m_nx; ++ix )
          for ( integer iy = 0; iy < m_ny; ++iy ) z_node_ref( ix, iy ) = Z( iy, ix );
      }
      else
      {
        UTILS_ASSERT(
          Z.rows() >= m_nx && Z.cols() >= m_ny,
          "SplineSurf::load_Z( Z, transposed={} ) bad dimension found {} x {} expected {} x {}",
          transposed,
          Z.rows(),
          Z.cols(),
          m_nx,
          m_ny );
        for ( integer ix = 0; ix < m_nx; ++ix )
          for ( integer iy = 0; iy < m_ny; ++iy ) z_node_ref( ix, iy ) = Z( ix, iy );
      }
      m_Z_max = Z.maxCoeff();
      m_Z_min = Z.minCoeff();
    }

    void load_Z( real_type const z[], integer const ldZ, bool fortran_storage, bool transposed );

    virtual void make_spline() = 0;

    template <typename Derived>
    void make_derivative_x( CubicSplineBase * S, Eigen::ArrayBase<Derived> const & Z, Eigen::ArrayBase<Derived> & DX )
    {
      for ( integer j = 0; j < m_ny; ++j )
      {
        S->build( mX, Z.col( j ) );
        for ( integer i = 0; i < m_nx; ++i ) DX( i, j ) = S->yp_node( i );
      }
    }

    template <typename Derived>
    void make_derivative_y( CubicSplineBase * S, Eigen::ArrayBase<Derived> const & Z, Eigen::ArrayBase<Derived> & DY )
    {
      for ( integer i = 0; i < m_nx; ++i )
      {
        S->build( mX, Z.row( i ) );
        for ( integer j = 0; j < m_ny; ++j ) DY( i, j ) = S->yp_node( j );
      }
    }

    template <typename Derived> void make_derivative_xy(
      CubicSplineBase *                 S,
      Eigen::ArrayBase<Derived> const & DX,
      Eigen::ArrayBase<Derived> const & DY,
      Eigen::ArrayBase<Derived> &       DXY )
    {
      auto minmod = []( real_type a, real_type b ) -> real_type
      {
        if ( a * b <= 0 ) return 0;
        if ( a > 0 ) return std::min( a, b );
        return std::max( a, b );
      };

      for ( integer j = 0; j < m_ny; ++j )
      {
        S->build( mX, DY.col( j ) );
        for ( integer i = 0; i < m_nx; ++i ) DXY( i, j ) = S->yp_node( i );
      }

      for ( integer i = 0; i < m_nx; ++i )
      {
        S->build( mY, DX.row( i ) );
        for ( integer j = 0; j < m_ny; ++j ) DXY( i, j ) = minmod( DXY( i, j ), S->yp_node( j ) );
      }
    }

    void resize( integer const nx, integer const ny );

  public:
    SplineSurf( SplineSurf const & )                   = delete;  // block copy constructor
    SplineSurf const & operator=( SplineSurf const & ) = delete;  // block copy method

    //!
    //! Spline constructor
    //!
    explicit SplineSurf( string_view name = "Spline" ) : m_name( name )
    {
      m_search_x.setup( &m_name, &m_nx, &m_X_ptr, &m_x_closed, &m_x_can_extend );
      m_search_y.setup( &m_name, &m_ny, &m_Y_ptr, &m_y_closed, &m_y_can_extend );
    }

    //!
    //! Spline destructor
    //!
    virtual ~SplineSurf() = default;

    //!
    //! \name Open/Close
    //!
    ///@{

    //!
    //! Return `true` if the surface is assumed closed in the `x` direction.
    //!
    [[nodiscard]] bool is_x_closed() const noexcept { return m_x_closed; }

    //!
    //! Setup the surface as closed in the `x` direction.
    //!
    void make_x_closed() { m_x_closed = true; }

    //!
    //! Setup the surface as open in the `x` direction.
    //!
    void make_x_opened() { m_x_closed = false; }

    //!
    //! Return `true` if the surface is assumed closed in the `y` direction.
    //!
    [[nodiscard]] bool is_y_closed() const noexcept { return m_y_closed; }

    //!
    //! Setup the surface as closed in the `y` direction.
    //!
    void make_y_closed() { m_y_closed = true; }

    //!
    //! Setup the surface as open in the `y` direction.
    //!
    void make_y_opened() { m_y_closed = false; }

    //!
    //! Return `true` if the parameter `x` assumed bounded.
    //! If false the spline is estrapolated for `x` values
    //! outside the range.
    //!
    [[nodiscard]] bool is_x_bounded() const noexcept { return !m_x_can_extend; }

    //!
    //! Make the spline surface unbounded in the `x` direction.
    //!
    void make_x_unbounded() { m_x_can_extend = true; }

    //!
    //! Make the spline surface bounded in the `x` direction.
    //!
    void make_x_bounded() { m_x_can_extend = false; }

    //!
    //! Return `true` if the parameter `y` assumed bounded.
    //! If false the spline is extrapolated for `y` values
    //! outside the range.
    //!
    [[nodiscard]] bool is_y_bounded() const noexcept { return !m_y_can_extend; }

    //!
    //! Make the spline surface unbounded in the `y` direction
    //!
    void make_y_unbounded() { m_y_can_extend = true; }

    //!
    //! Make the spline surface bounded in the `x` direction.
    //!
    void make_y_bounded() { m_y_can_extend = false; }

    ///@}

    //!
    //! Cancel the support points, empty the spline.
    //!
    void clear();

    //!
    //! \name Info
    //!
    ///@{

    //!
    //! \return string with the name of the spline
    //!
    [[nodiscard]] string_view name() const noexcept { return m_name; }

    //!
    //! Return the number of support points of the spline along x direction.
    //!
    [[nodiscard]] integer num_point_x() const noexcept { return m_nx; }

    //!
    //! Return the number of support points of the spline along y direction.
    //!
    [[nodiscard]] integer num_point_y() const noexcept { return m_ny; }

    //!
    //! Return the i-th node of the spline (x component).
    //!
    [[nodiscard]] real_type x_node( integer const i ) const { return mX.coeff( i ); }

    //!
    //! Return the i-th node of the spline (y component).
    //!
    [[nodiscard]] real_type y_node( integer const i ) const { return mY.coeff( i ); }

    //!
    //! Return the i-th node of the spline (y component).
    //!
    [[nodiscard]] real_type z_node( integer const i, integer const j ) const { return mZ.coeff( i, j ); }

    //!
    //! Return x-minumum spline value.
    //!
    [[nodiscard]] real_type x_min() const { return mX.coeff( 0 ); }

    //!
    //! Return x-maximum spline value.
    //!
    [[nodiscard]] real_type x_max() const { return mX.coeff( mX.size() - 1 ); }

    //!
    //! Return y-minumum spline value.
    //!
    [[nodiscard]] real_type y_min() const { return mY.coeff( 0 ); }

    //!
    //! Return y-maximum spline value.
    //!
    [[nodiscard]] real_type y_max() const { return mY.coeff( mY.size() - 1 ); }

    //!
    //! Return z-minumum spline value.
    //!
    [[nodiscard]] real_type z_min() const noexcept { return m_Z_min; }

    //!
    //! Return z-maximum spline value.
    //!
    [[nodiscard]] real_type z_max() const noexcept { return m_Z_max; }

    ///@}

    //!
    //! \name Build Spline
    //!
    ///@{


    ///
    /// Build the surface spline from generic coordinate/value accessors.
    ///
    /// The arguments `X`, `Y`, and `Z` must be callable objects.
    ///
    /// Required interface:
    /// - `X(ix)` returns the x-coordinate of node `ix`, for `0 <= ix < nx`
    /// - `Y(iy)` returns the y-coordinate of node `iy`, for `0 <= iy < ny`
    /// - `Z(ix, iy)` returns the scalar value associated with `X(ix), Y(iy)`
    ///
    /// The callable objects may be:
    /// - lambdas
    /// - functors
    /// - wrappers around vectors/matrices
    /// - any object accepted by `std::invoke`
    ///
    /// \tparam XFunc  callable type for x-coordinates
    /// \tparam YFunc  callable type for y-coordinates
    /// \tparam ZFunc  callable type for surface values
    ///
    /// \param X   accessor for x-coordinates
    /// \param Y   accessor for y-coordinates
    /// \param Z   accessor for surface values
    /// \param nx  number of grid points along the x-direction
    /// \param ny  number of grid points along the y-direction
    ///
    template <typename XFunc, typename YFunc, typename ZFunc>
    void build( XFunc && X, YFunc && Y, ZFunc && Z, integer const nx, integer const ny )
    {
      static_assert(
        std::is_invocable_r_v<real_type, XFunc, integer>,
        "build(): X(i) must be callable and return real_type" );

      static_assert(
        std::is_invocable_r_v<real_type, YFunc, integer>,
        "build(): Y(j) must be callable and return real_type" );

      static_assert(
        std::is_invocable_r_v<real_type, ZFunc, integer, integer>,
        "build(): Z(ix, iy) must be callable and return real_type" );

      resize( nx, ny );

      for ( integer i = 0; i < nx; ++i ) { mX.coeffRef( i ) = std::invoke( std::forward<XFunc>( X ), i ); }

      for ( integer j = 0; j < ny; ++j ) { mY.coeffRef( j ) = std::invoke( std::forward<YFunc>( Y ), j ); }

      for ( integer ix = 0; ix < m_nx; ++ix )
      {
        for ( integer iy = 0; iy < m_ny; ++iy )
        {
          z_node_ref( ix, iy ) = std::invoke( std::forward<ZFunc>( Z ), ix, iy );
        }
      }

      make_spline();
    }

    //!
    //! Build surface spline
    //!
    //! \param x               vector of `x`-coordinates
    //! \param incx            access elements as `x[0]`, `x[incx]`, `x[2*incx]`,...
    //! \param y               vector of `y`-coordinates
    //! \param incy            access elements as `y[0]`, `y[incy]`, `y[2*incy]`,...
    //! \param z               matrix of `z`-values. Elements are stored
    //!                        by row Z(i,j) = z[i*ny+j] as C-matrix
    //! \param ldZ             leading dimension of `z`
    //! \param nx              number of points in `x` direction
    //! \param ny              number of points in `y` direction
    //! \param fortran_storage if true elements are stored by column
    //!                        i.e. Z(i,j) = z[i+j*nx] as Fortran-matrix
    //! \param transposed      if true matrix Z is stored transposed
    //!
    void build(
      real_type const x[],
      integer const   incx,
      real_type const y[],
      integer const   incy,
      real_type const z[],
      integer const   ldZ,
      integer const   nx,
      integer const   ny,
      bool            fortran_storage = false,
      bool            transposed      = false );

    //!
    //! Build surface spline
    //!
    //! \param x               vector of `x`-coordinates
    //! \param y               vector of `y`-coordinates
    //! \param z               matrix of `z`-values. Elements are stored
    //!                        by row Z(i,j) = z[i*ny+j] as C-matrix
    //! \param nx              number of points in `x` direction
    //! \param ny              number of points in `y` direction
    //! \param fortran_storage if true elements are stored by column
    //!                        i.e. Z(i,j) = z[i+j*nx] as Fortran-matrix
    //! \param transposed      if true matrix Z is stored transposed
    //!
    void build(
      real_type const x[],
      real_type const y[],
      real_type const z[],
      integer const   nx,
      integer const   ny,
      bool            fortran_storage = false,
      bool            transposed      = false )
    {
      integer ldZ = fortran_storage ? nx : ny;
      build( x, 1, y, 1, z, ldZ, nx, ny, fortran_storage, transposed );
    }

    //!
    //! Build surface spline
    //!
    //! \param x               vector of x-coordinates, nx = x.size()
    //! \param y               vector of y-coordinates, ny = y.size()
    //! \param Z               matrix of z-values. Elements are stored
    //!                        by row Z(i,j) = z[i*ny+j] as C-matrix
    //! \param transposed      if true matrix Z is stored transposed
    //!
    template <typename Derived> void build(
      Eigen::Ref<const Vec>             x,
      Eigen::Ref<const Vec>             y,
      Eigen::ArrayBase<Derived> const & Z,
      bool const                        transposed )
    {
      resize( x.size(), y.size() );
      mX = x;
      mY = y;
      load_Z( Z, transposed );
      make_spline();
    }

    //!
    //! Build surface spline
    //!
    //! \param x               vector of x-coordinates, nx = x.size()
    //! \param y               vector of y-coordinates, ny = y.size()
    //! \param z               matrix of z-values. Elements are stored
    //!                        by row Z(i,j) = z[i*ny+j] as C-matrix
    //! \param fortran_storage if true elements are stored by column
    //!                        i.e. Z(i,j) = z[i+j*nx] as Fortran-matrix
    //! \param transposed      if true matrix Z is stored transposed
    //!
    void build(
      vector<real_type> const & x,
      vector<real_type> const & y,
      vector<real_type> const & z,
      bool                      fortran_storage = false,
      bool                      transposed      = false )
    {
      integer nx  = static_cast<integer>( x.size() );
      integer ny  = static_cast<integer>( y.size() );
      size_t  nz  = static_cast<size_t>( nx ) * static_cast<size_t>( ny );
      UTILS_ASSERT(
        z.size() == nz,
        "SplineSurf::build( x, y, z, ... ) bad z size found {} expected {} = {} x {}",
        z.size(),
        nz,
        nx,
        ny );
      integer ldZ = fortran_storage ? nx : ny;
      build( x.data(), 1, y.data(), 1, z.data(), ldZ, nx, ny, fortran_storage, transposed );
    }

    void build(
      real_type const z[],
      integer         ldZ,
      integer         nx,
      integer         ny,
      bool            fortran_storage = false,
      bool            transposed      = false );

    //!
    //! Build surface spline
    //!
    //! \param z               matrix of z-values. Elements are stored
    //!                        by row Z(i,j) = z[i*ny+j] as C-matrix.
    //!                        ldZ leading dimension of the matrix is ny for C-storage
    //!                        and nx for Fortran storage.
    //! \param nx              x-dimension
    //! \param ny              y-dimension
    //! \param fortran_storage if true elements are stored by column
    //!                        i.e. Z(i,j) = z[i+j*nx] as Fortran-matrix
    //! \param transposed      if true matrix Z is stored transposed
    //!
    void build(
      vector<real_type> const & z,
      integer const             nx,
      integer const             ny,
      bool                      fortran_storage = false,
      bool                      transposed      = false )
    {
      size_t nz = static_cast<size_t>( nx ) * static_cast<size_t>( ny );
      UTILS_ASSERT(
        z.size() == nz,
        "SplineSurf::build( z, nx, ny, ... ) bad z size found {} expected {} = {} x {}",
        z.size(),
        nz,
        nx,
        ny );
      integer ldZ = fortran_storage ? nx : ny;
      this->build( z.data(), ldZ, nx, ny, fortran_storage, transposed );
    }

    //!
    //! Build spline using data in `gc`
    //!
    void setup( GenericContainer const & gc );

    //!
    //! Setup a spline using a `GenericContainer` readed from file
    //!
    //! - file_name file name of the file with data
    //!
    //! - `.json`
    //! - `.yaml` or `.yml`
    //! - `.toml`
    //!
    void setup( string const & file_name )
    {
      GenericContainer gc;
      UTILS_ASSERT( gc.from_file( file_name ), "Spline::setup( '{}' ) failed to read\n", file_name );
      setup( gc );
    }

    //!
    //! Build a spline using data in `GenericContainer`
    //!
    void build( GenericContainer const & gc ) { setup( gc ); }

    //!
    //! Build a spline using data in a file `file_name`
    //!
    void build( string const & file_name ) { setup( file_name ); }

    ///@}

    //!
    //! \name Evaluate
    //!
    ///@{

    //!
    //! Evaluate spline value at point \f$ (x,y) \f$.
    //!
    [[nodiscard]] virtual real_type eval( real_type const x, real_type const y ) const = 0;

#ifdef AUTODIFF_SUPPORT
    // Metodi base per dual1st e dual2nd
    [[nodiscard]] virtual autodiff::dual1st eval( autodiff::dual1st const & x, autodiff::dual1st const & y ) const = 0;
    [[nodiscard]] virtual autodiff::dual2nd eval( autodiff::dual2nd const & x, autodiff::dual2nd const & y ) const = 0;

    // Template per due parametri (x, y) - per SplineSurf
    // Promuove automaticamente a double, dual1st o dual2nd in base ai tipi di input
    template <typename T1, typename T2> [[nodiscard]] auto eval( T1 const & x, T2 const & y ) const
    {
      if constexpr ( std::is_arithmetic<T1>::value && std::is_arithmetic<T2>::value )
      {
        // Entrambi numerici: ritorna real_type
        return eval( static_cast<real_type>( x ), static_cast<real_type>( y ) );
      }
      else
      {
        // Almeno uno è duale: determina il tipo duale appropriato
        constexpr int order1    = std::is_arithmetic<T1>::value ? 0 : autodiff::detail::DualOrder<T1>::value;
        constexpr int order2    = std::is_arithmetic<T2>::value ? 0 : autodiff::detail::DualOrder<T2>::value;
        constexpr int max_order = ( order1 > order2 ) ? order1 : order2;

        if constexpr ( max_order == 1 )
        {
          // Promuovi a dual1st
          autodiff::dual1st X = std::is_arithmetic<T1>::value ? autodiff::dual1st{ static_cast<real_type>( x ) }
                                                              : autodiff::dual1st{ x };
          autodiff::dual1st Y = std::is_arithmetic<T2>::value ? autodiff::dual1st{ static_cast<real_type>( y ) }
                                                              : autodiff::dual1st{ y };
          return eval( X, Y );
        }
        else
        {
          // Promuovi a dual2nd
          autodiff::dual2nd X = std::is_arithmetic<T1>::value ? autodiff::dual2nd{ static_cast<real_type>( x ) }
                                                              : autodiff::dual2nd{ x };
          autodiff::dual2nd Y = std::is_arithmetic<T2>::value ? autodiff::dual2nd{ static_cast<real_type>( y ) }
                                                              : autodiff::dual2nd{ y };
          return eval( X, Y );
        }
      }
    }

    // Operator() per due parametri
    template <typename T1, typename T2>
    [[nodiscard]] auto operator()( T1 const & x, T2 const & y ) const -> decltype( eval( x, y ) )
    { return this->eval( x, y ); }
#endif

    //!
    //! Value and first derivatives at point \f$ (x,y) \f$:
    //!
    //! - d[0] value of the spline \f$ S(x,y) \f$
    //! - d[1] derivative respect to \f$ x \f$ of the spline: \f$ S_x(x,y) \f$
    //! - d[2] derivative respect to \f$ y \f$ of the spline: \f$ S_y(x,y) \f$
    //!
    virtual void D( real_type const x, real_type const y, real_type d[3] ) const = 0;

    //!
    //! First derivatives respect to \f$ x \f$ at point \f$ (x,y) \f$
    //! of the spline: \f$ S_x(x,y) \f$.
    //!
    [[nodiscard]] virtual real_type Dx( real_type const x, real_type const y ) const = 0;

    //!
    //! First derivatives respect to \f$ y \f$ at point \f$ (x,y) \f$
    //! of the spline: \f$ S_y(x,y) \f$.
    //!
    [[nodiscard]] virtual real_type Dy( real_type const x, real_type const y ) const = 0;

    //!
    //! Value, first and second derivatives at point \f$ (x,y) \f$:
    //!
    //! - dd[0] value of the spline \f$ S(x,y) \f$
    //! - dd[1] derivative respect to \f$ x \f$ of the spline: \f$ S_x(x,y) \f$
    //! - dd[2] derivative respect to \f$ y \f$ of the spline: \f$ S_y(x,y) \f$
    //! - dd[3] second derivative respect to \f$ x \f$ of the spline: \f$ S_{xx}(x,y) \f$
    //! - dd[4] mixed second derivative: \f$ S_{xy}(x,y) \f$
    //! - dd[5] second derivative respect to \f$ y \f$ of the spline: \f$ S_{yy}(x,y) \f$
    //!
    virtual void DD( real_type const x, real_type const y, real_type dd[6] ) const = 0;

    //!
    //! Second derivatives respect to \f$ x \f$ at point \f$ (x,y) \f$
    //! of the spline: \f$ S_{xx}(x,y) \f$.
    //!
    [[nodiscard]] virtual real_type Dxx( real_type const x, real_type const y ) const = 0;

    //!
    //! Mixed second derivatives: \f$ S_{xy}(x,y) \f$.
    //!
    [[nodiscard]] virtual real_type Dxy( real_type const x, real_type const y ) const = 0;

    //!
    //! Second derivatives respect to \f$ y \f$ at point \f$ (x,y) \f$
    //! of the spline: \f$ S_{yy}(x,y) \f$.
    //!
    [[nodiscard]] virtual real_type Dyy( real_type const x, real_type const y ) const = 0;

    //!
    //! Evaluate spline value at point \f$ (x,y) \f$.
    //!
    [[nodiscard]] real_type operator()( real_type const x, real_type const y ) const { return this->eval( x, y ); }

    //!
    //! Alias for `Dx(x,y)`
    //!
    [[nodiscard]] real_type eval_D_1( real_type const x, real_type const y ) const { return this->Dx( x, y ); }

    //!
    //! Alias for `Dy(x,y)`
    //!
    [[nodiscard]] real_type eval_D_2( real_type const x, real_type const y ) const { return this->Dy( x, y ); }

    //!
    //! Alias for `Dxx(x,y)`
    //!
    [[nodiscard]] real_type eval_D_1_1( real_type const x, real_type const y ) const { return this->Dxx( x, y ); }

    //!
    //! Alias for `Dxy(x,y)`
    //!
    [[nodiscard]] real_type eval_D_1_2( real_type const x, real_type const y ) const { return this->Dxy( x, y ); }

    //!
    //! Alias for `Dyy(x,y)`
    //!
    [[nodiscard]] real_type eval_D_2_2( real_type const x, real_type const y ) const { return this->Dyy( x, y ); }

    ///@}

    //!
    //! Print spline coefficients.
    //!
    virtual void write_to_stream( ostream_type & s ) const = 0;

    //!
    //! Return spline type as a string pointer.
    //!
    [[nodiscard]] virtual char const * type_name() const = 0;

    //!
    //! String information of the kind and order of the spline
    //!
    [[nodiscard]] virtual string info() const { return fmt::format( "Bivariate spline [{}] of type = {}", name(), type_name() ); }

    //!
    //! Print information of the kind and order of the spline
    //!
    void info( ostream_type & stream ) const { stream << this->info() << '\n'; }

    //!
    //! Print stored data x, y, and matrix z.
    //!
    void dump_data( ostream_type & s ) const
    {
      s << "X = [ " << x_node( 0 );
      for ( integer i = 1; i < m_nx; ++i ) s << ", " << x_node( i );
      s << " ]\nY = [ " << y_node( 0 );
      for ( integer j = 1; j < m_ny; ++j ) s << ", " << y_node( j );
      s << " ]\nZ = [\n";
      for ( integer j = 0; j < m_ny; ++j )
      {
        s << "  [ " << z_node( 0, j );
        for ( integer i = 1; i < m_nx; ++i ) s << ", " << z_node( i, j );
        s << " ]\n";
      }
      s << "\n];\n";
    }

#ifdef SPLINES_BACK_COMPATIBILITY
    [[nodiscard]] integer   numPointX() const noexcept { return m_nx; }
    [[nodiscard]] integer   numPointY() const noexcept { return m_ny; }
    [[nodiscard]] real_type xNode( integer i ) const { return m_X[i]; }
    [[nodiscard]] real_type yNode( integer i ) const { return m_Y[i]; }
    [[nodiscard]] real_type zNode( integer i, integer j ) const { return z_node( i, j ); }
    [[nodiscard]] real_type xMin() const { return this->x_min(); }
    [[nodiscard]] real_type xMax() const { return this->x_max(); }
    [[nodiscard]] real_type yMin() const { return this->y_min(); }
    [[nodiscard]] real_type yMax() const { return this->y_max(); }
    [[nodiscard]] real_type zMin() const noexcept { return m_Z_min; }
    [[nodiscard]] real_type zMax() const noexcept { return m_Z_max; }
    void                    writeToStream( ostream_type & s ) const { write_to_stream( s ); }
#endif
  };

}  // namespace Splines

#endif

// EOF: SplineSurf.hxx
