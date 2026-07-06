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

#ifndef SPLINES2D_HXX
#define SPLINES2D_HXX

namespace Splines
{

  /*\
   |   ____        _ _            ____  ____
   |  / ___| _ __ | (_)_ __   ___|___ \|  _ \
   |  \___ \| '_ \| | | '_ \ / _ \ __) | | | |
   |   ___) | |_) | | | | | |  __// __/| |_| |
   |  |____/| .__/|_|_|_| |_|\___|_____|____/
   |        |_|
  \*/

  //!
  //! \class Spline2D
  //! \brief A wrapper class for 2D spline surfaces (bilinear, bicubic, biquintic, etc.)
  //!
  //! This class provides a unified interface for various types of 2D spline surfaces.
  //! It manages the construction, evaluation, and querying of spline surfaces defined
  //! over a rectangular grid in the x-y plane.
  //!
  //! The surface is defined by:
  //!   - A set of nodal points (x_i, y_j) forming a rectangular grid.
  //!   - A matrix of z-values Z(i,j) defined at each grid point.
  //!
  //! The Z matrix can be stored in different memory layouts (row-major or column-major)
  //! and can optionally be transposed. See the build() methods for detailed storage descriptions.
  //!
  class Spline2D final
  {
  protected:
    std::string  m_name;                 ///< Name identifier for the spline surface
    SplineSurf * m_spline_2D = nullptr;  ///< Pointer to the actual spline implementation

    //! \brief Internal method to create a new spline of the specified type
    //! \param[in] tp The type of spline surface to create (e.g., bilinear, bicubic, biquintic)
    void new_spline( SplineType2D tp )
    {
      if ( m_spline_2D != nullptr )
      {
        delete m_spline_2D;
        m_spline_2D = nullptr;
      }
      switch ( tp )
      {
        case SplineType2D::BILINEAR: m_spline_2D = new BilinearSpline( m_name ); break;
        case SplineType2D::BICUBIC_CUBIC: m_spline_2D = new BiCubicSpline( Spline_sub_type::CUBIC, m_name ); break;
        case SplineType2D::BICUBIC_AKIMA: m_spline_2D = new BiCubicSpline( Spline_sub_type::AKIMA, m_name ); break;
        case SplineType2D::BICUBIC_VANLEER: m_spline_2D = new BiCubicSpline( Spline_sub_type::VANLEER, m_name ); break;
        case SplineType2D::BICUBIC_PCHIP: m_spline_2D = new BiCubicSpline( Spline_sub_type::PCHIP, m_name ); break;
        case SplineType2D::BIQUINTIC_CUBIC: m_spline_2D = new BiQuinticSpline( Spline_sub_type::CUBIC, m_name ); break;
        case SplineType2D::BIQUINTIC_AKIMA: m_spline_2D = new BiQuinticSpline( Spline_sub_type::AKIMA, m_name ); break;
        case SplineType2D::BIQUINTIC_VANLEER:
          m_spline_2D = new BiQuinticSpline( Spline_sub_type::VANLEER, m_name );
          break;
        case SplineType2D::BIQUINTIC_PCHIP: m_spline_2D = new BiQuinticSpline( Spline_sub_type::PCHIP, m_name ); break;
      }
    }

  public:
    //! \name Constructors and Destructor
    ///@{

    //! \brief Constructs an empty spline surface with the given name.
    //! \param[in] name Optional name identifier for the spline (default: "Spline2D")
    explicit Spline2D( string_view name = "Spline2D" ) : m_name( name ) {}

    //! \brief Destructor. Safely deletes the internal spline object.
    ~Spline2D()
    {
      if ( m_spline_2D != nullptr )
      {
        delete m_spline_2D;
        m_spline_2D = nullptr;
      }
    }

    ///@}

    //! \name Open/Close/Boundary Conditions
    //! Methods to control periodicity and extrapolation behavior.
    ///@{

    //! \brief Returns true if the surface is closed (periodic) in the x-direction.
    [[nodiscard]] bool is_x_closed() const { return m_spline_2D->is_x_closed(); }

    //! \brief Makes the surface closed (periodic) in the x-direction.
    void make_x_closed() { m_spline_2D->make_x_closed(); }

    //! \brief Makes the surface open (non‑periodic) in the x‑direction.
    void make_x_opened() { m_spline_2D->make_x_opened(); }

    //! \brief Returns true if the surface is closed (periodic) in the y-direction.
    [[nodiscard]] bool is_y_closed() const { return m_spline_2D->is_y_closed(); }

    //! \brief Makes the surface closed (periodic) in the y-direction.
    void make_y_closed() { m_spline_2D->make_y_closed(); }

    //! \brief Makes the surface open (non‑periodic) in the y‑direction.
    void make_y_opened() { m_spline_2D->make_y_opened(); }

    //! \brief Returns true if x-values are bounded (no extrapolation).
    //! If false, the spline will extrapolate for x outside the defined range.
    [[nodiscard]] bool is_x_bounded() const { return m_spline_2D->is_x_bounded(); }

    //! \brief Allows extrapolation in the x‑direction (unbounded).
    void make_x_unbounded() { m_spline_2D->make_x_unbounded(); }

    //! \brief Restricts evaluation to the defined x‑range (bounded).
    void make_x_bounded() { m_spline_2D->make_x_bounded(); }

    //! \brief Returns true if y-values are bounded (no extrapolation).
    //! If false, the spline will extrapolate for y outside the defined range.
    [[nodiscard]] bool is_y_bounded() const { return m_spline_2D->is_y_bounded(); }

    //! \brief Allows extrapolation in the y‑direction (unbounded).
    void make_y_unbounded() { m_spline_2D->make_y_unbounded(); }

    //! \brief Restricts evaluation to the defined y‑range (bounded).
    void make_y_bounded() { m_spline_2D->make_y_bounded(); }

    ///@}

    //! \name Information and Accessors
    //! Methods to query spline properties and nodal data.
    ///@{

    //! \brief Returns the name of the spline surface.
    [[nodiscard]] string_view name() const { return m_spline_2D->name(); }

    //! \brief Returns the number of grid points in the x-direction.
    [[nodiscard]] integer num_point_x() const { return m_spline_2D->num_point_x(); }

    //! \brief Returns the number of grid points in the y-direction.
    [[nodiscard]] integer num_point_y() const { return m_spline_2D->num_point_y(); }

    //! \brief Returns the x-coordinate of the i-th grid node.
    //! \param[in] i Index of the node in the x-direction (0‑based).
    [[nodiscard]] real_type x_node( integer const i ) const { return m_spline_2D->x_node( i ); }

    //! \brief Returns the y-coordinate of the i-th grid node.
    //! \param[in] i Index of the node in the y-direction (0‑based).
    [[nodiscard]] real_type y_node( integer const i ) const { return m_spline_2D->y_node( i ); }

    //! \brief Returns the z-value at grid node (i, j).
    //! \param[in] i Index in the x-direction (0‑based).
    //! \param[in] j Index in the y-direction (0‑based).
    [[nodiscard]] real_type z_node( integer const i, integer const j ) const { return m_spline_2D->z_node( i, j ); }

    ///@}

    //! \brief Clears all support points, leaving the spline empty.
    void clear() { m_spline_2D->clear(); }

    //! \name Bounding Box Queries
    //! Methods to get the min/max ranges of the spline domain and values.
    ///@{

    //! \brief Minimum x-coordinate of the spline domain.
    [[nodiscard]] real_type x_min() const { return m_spline_2D->x_min(); }

    //! \brief Maximum x-coordinate of the spline domain.
    [[nodiscard]] real_type x_max() const { return m_spline_2D->x_max(); }

    //! \brief Minimum y-coordinate of the spline domain.
    [[nodiscard]] real_type y_min() const { return m_spline_2D->y_min(); }

    //! \brief Maximum y-coordinate of the spline domain.
    [[nodiscard]] real_type y_max() const { return m_spline_2D->y_max(); }

    //! \brief Minimum z-value over all grid points.
    [[nodiscard]] real_type z_min() const { return m_spline_2D->z_min(); }

    //! \brief Maximum z-value over all grid points.
    [[nodiscard]] real_type z_max() const { return m_spline_2D->z_max(); }

    ///@}

    //! \name Construction Methods
    //! Methods to build the spline surface from data.
    //! The Z matrix can be stored in several layouts:
    //!   - Row‑major (C style): Z(i,j) = z[i * ny + j] with ldZ = ny.
    //!   - Column‑major (Fortran style): Z(i,j) = z[i + j * nx] with ldZ = nx.
    //!   - Transposed: swaps the roles of rows and columns.
    //! See each overload for specific parameter descriptions.
    ///@{

    //! \brief Builds a spline surface from raw arrays.
    //!
    //! \param[in] tp              Type of spline surface (e.g., bilinear, bicubic).
    //! \param[in] x               Array of x-coordinates (length nx).
    //! \param[in] incx            Stride between consecutive x elements (typically 1).
    //! \param[in] y               Array of y-coordinates (length ny).
    //! \param[in] incy            Stride between consecutive y elements (typically 1).
    //! \param[in] z               Flat array containing the matrix of z‑values.
    //! \param[in] ldZ             Leading dimension of the Z matrix in memory.
    //!                            For C storage (row‑major), ldZ is the number of columns (ny).
    //!                            For Fortran storage (column‑major), ldZ is the number of rows (nx).
    //! \param[in] nx              Number of grid points in the x‑direction.
    //! \param[in] ny              Number of grid points in the y‑direction.
    //! \param[in] fortran_storage If true, z is stored column‑major (Fortran style).
    //!                            If false, z is stored row‑major (C style).
    //! \param[in] transposed      If true, the Z matrix is stored transposed.
    //!                            When transposed, the indexing changes:
    //!                            - C storage: Z(i,j) = z[j * ldZ + i] (with ldZ = nx).
    //!                            - Fortran storage: Z(i,j) = z[j + i * ldZ] (with ldZ = ny).
    //!
    //! Example (C storage, not transposed):
    //!   Given nx=3, ny=2, ldZ=ny=2, the flat array z[6] contains:
    //!     [ Z(0,0), Z(0,1), Z(1,0), Z(1,1), Z(2,0), Z(2,1) ]
    //!   Access: Z(i,j) = z[i*2 + j].
    //!
    void build(
      SplineType2D    tp,
      real_type const x[],
      integer         incx,
      real_type const y[],
      integer         incy,
      real_type const z[],
      integer         ldZ,
      integer const   nx,
      integer const   ny,
      bool            fortran_storage = false,
      bool            transposed      = false )
    {
      new_spline( tp );
      m_spline_2D->build( x, incx, y, incy, z, ldZ, nx, ny, fortran_storage, transposed );
    }

    //! \brief Builds a spline surface from std::vector containers.
    //!
    //! \param[in] tp              Type of spline surface.
    //! \param[in] x               Vector of x-coordinates (length nx).
    //! \param[in] y               Vector of y-coordinates (length ny).
    //! \param[in] z               Vector containing the matrix of z‑values, stored as a flat array.
    //! \param[in] fortran_storage If true, z is stored column‑major (Fortran style).
    //!                            If false, z is stored row‑major (C style).
    //! \param[in] transposed      If true, the Z matrix is stored transposed.
    //!
    //! The vector z must have size nx * ny. Storage and transposition behave as in the raw‑array version.
    void build(
      SplineType2D              tp,
      vector<real_type> const & x,
      vector<real_type> const & y,
      vector<real_type> const & z,
      bool                      fortran_storage = false,
      bool                      transposed      = false )
    {
      new_spline( tp );
      m_spline_2D->build( x, y, z, fortran_storage, transposed );
    }

    //! \brief Builds a spline surface with uniformly spaced x and y grids.
    //!
    //! The x-grid is assumed to be 0, 1, ..., nx-1.
    //! The y-grid is assumed to be 0, 1, ..., ny-1.
    //!
    //! \param[in] tp              Type of spline surface.
    //! \param[in] z               Flat array containing the matrix of z‑values.
    //! \param[in] ldZ             Leading dimension of Z in memory.
    //!                            For C storage (row‑major): ldZ = ny.
    //!                            For Fortran storage (column‑major): ldZ = nx.
    //! \param[in] nx              Number of grid points in the x‑direction.
    //! \param[in] ny              Number of grid points in the y‑direction.
    //! \param[in] fortran_storage If true, z is stored column‑major.
    //! \param[in] transposed      If true, the Z matrix is stored transposed.
    void build(
      SplineType2D    tp,
      real_type const z[],
      integer const   ldZ,
      integer const   nx,
      integer const   ny,
      bool            fortran_storage = false,
      bool            transposed      = false )
    {
      new_spline( tp );
      m_spline_2D->build( z, ldZ, nx, ny, fortran_storage, transposed );
    }

    //! \brief Builds a spline surface with uniformly spaced grids using a vector for Z.
    //!
    //! The x-grid is 0, 1, ..., nx-1; the y-grid is 0, 1, ..., ny-1.
    //!
    //! \param[in] tp              Type of spline surface.
    //! \param[in] z               Vector containing the matrix of z‑values (size must be nx * ny).
    //! \param[in] nx              Number of grid points in the x‑direction.
    //! \param[in] ny              Number of grid points in the y‑direction.
    //! \param[in] fortran_storage If true, z is stored column‑major.
    //! \param[in] transposed      If true, the Z matrix is stored transposed.
    void build(
      SplineType2D              tp,
      vector<real_type> const & z,
      integer const             nx,
      integer const             ny,
      bool                      fortran_storage = false,
      bool                      transposed      = false )
    {
      new_spline( tp );
      m_spline_2D->build( z, nx, ny, fortran_storage, transposed );
    }

    //! \brief Configures the spline from a GenericContainer.
    //!
    //! The GenericContainer must contain a key "spline_type" with a string value:
    //!   - "bilinear"  : bilinear spline surface.
    //!   - "bicubic"   : bicubic spline surface.
    //!   - "biquintic" : biquintic spline surface.
    //!   - "Akima" or "akima": cubic spline using Akima algorithm.
    //!
    //! Other parameters (grid and Z data) are also read from the container.
    //! \param[in] gc The GenericContainer holding the spline configuration.
    void setup( GenericContainer const & gc )
    {
      string const   where = fmt::format( "Spline2D[{}]::setup( gc ):", m_name );
      string const & type  = gc.get_map_string( "spline_type", where );
      new_spline( string_to_splineType2D( type ) );
      m_spline_2D->setup( gc );
    }

    //! \brief Alias for setup().
    void build( GenericContainer const & gc ) { setup( gc ); }

    ///@}

    //! \name Evaluation
    //! Methods to evaluate the spline surface and its derivatives.
    ///@{

    //! \brief Evaluates the spline at (x, y).
    //! \param[in] x x-coordinate.
    //! \param[in] y y-coordinate.
    //! \return The interpolated value S(x, y).
    real_type operator()( real_type const x, real_type const y ) const { return m_spline_2D->eval( x, y ); }

    //! \brief Evaluates the spline at (x, y).
    //! \param[in] x x-coordinate.
    //! \param[in] y y-coordinate.
    //! \return The interpolated value S(x, y).
    real_type eval( real_type const x, real_type const y ) const { return m_spline_2D->eval( x, y ); }

#ifdef AUTODIFF_SUPPORT

    //! \brief Evaluates the spline using first‑order automatic differentiation.
    autodiff::dual1st eval( autodiff::dual1st const & x, autodiff::dual1st const & y ) const
    { return m_spline_2D->eval( x, y ); }
    //! \brief Evaluates the spline using second‑order automatic differentiation.
    autodiff::dual2nd eval( autodiff::dual2nd const & x, autodiff::dual2nd const & y ) const
    { return m_spline_2D->eval( x, y ); }

    //! \brief Evaluates the spline with higher‑order automatic differentiation.
    template <typename T1, typename T2> auto eval( T1 const & x, T2 const & y ) const
      -> decltype( m_spline_2D->eval( x, y ) )
    { return m_spline_2D->eval( x, y ); }

    // Operator() per due parametri
    template <typename T1, typename T2> auto operator()( T1 const & x, T2 const & y ) const
      -> decltype( m_spline_2D->eval( x, y ) )
    { return m_spline_2D->eval( x, y ); }

#endif

    ///@}

    //! \name First Derivatives
    //! Methods to compute first derivatives of the spline surface.
    ///@{

    //! \brief Computes value and first derivatives at (x, y).
    //! \param[in]  x  x-coordinate.
    //! \param[in]  y  y-coordinate.
    //! \param[out] d  Array of size 3 where:
    //!                - d[0] = S(x, y)
    //!                - d[1] = ∂S/∂x (x, y)
    //!                - d[2] = ∂S/∂y (x, y)
    void D( real_type const x, real_type const y, real_type d[3] ) const { return m_spline_2D->D( x, y, d ); }

    //! \brief Returns ∂S/∂x at (x, y).
    real_type Dx( real_type const x, real_type const y ) const { return m_spline_2D->Dx( x, y ); }

    //! \brief Returns ∂S/∂y at (x, y).
    real_type Dy( real_type const x, real_type const y ) const { return m_spline_2D->Dy( x, y ); }

    //! \brief Alias for Dx().
    real_type eval_D_1( real_type const x, real_type const y ) const { return this->Dx( x, y ); }

    //! \brief Alias for Dy().
    real_type eval_D_2( real_type const x, real_type const y ) const { return this->Dy( x, y ); }

    ///@}

    //! \name Second Derivatives
    //! Methods to compute second derivatives of the spline surface.
    ///@{

    //! \brief Computes value, first and second derivatives at (x, y).
    //! \param[in]  x   x-coordinate.
    //! \param[in]  y   y-coordinate.
    //! \param[out] dd  Array of size 6 where:
    //!                 - dd[0] = S(x, y)
    //!                 - dd[1] = ∂S/∂x
    //!                 - dd[2] = ∂S/∂y
    //!                 - dd[3] = ∂²S/∂x²
    //!                 - dd[4] = ∂²S/∂x∂y
    //!                 - dd[5] = ∂²S/∂y²
    void DD( real_type const x, real_type const y, real_type dd[6] ) const { return m_spline_2D->DD( x, y, dd ); }

    //! \brief Returns ∂²S/∂x² at (x, y).
    real_type Dxx( real_type const x, real_type const y ) const { return m_spline_2D->Dxx( x, y ); }

    //! \brief Returns ∂²S/∂x∂y at (x, y).
    real_type Dxy( real_type const x, real_type const y ) const { return m_spline_2D->Dxy( x, y ); }

    //! \brief Returns ∂²S/∂y² at (x, y).
    real_type Dyy( real_type const x, real_type const y ) const { return m_spline_2D->Dyy( x, y ); }

    //! \brief Alias for Dxx().
    real_type eval_D_1_1( real_type const x, real_type const y ) const { return this->Dxx( x, y ); }

    //! \brief Alias for Dxy().
    real_type eval_D_1_2( real_type const x, real_type const y ) const { return this->Dxy( x, y ); }

    //! \brief Alias for Dyy().
    real_type eval_D_2_2( real_type const x, real_type const y ) const { return this->Dyy( x, y ); }
    ///@}

    //! \name Output and Debugging
    //! Methods for printing spline information and data.
    ///@{

    //! \brief Writes a human‑readable representation of the spline to a stream.
    void write_to_stream( ostream_type & s ) const { return m_spline_2D->write_to_stream( s ); }

    //! \brief Returns the spline type name as a C‑string.
    char const * type_name() const { return m_spline_2D->type_name(); }

    //! \brief Returns a string describing the spline kind and order.
    string info() const { return m_spline_2D->info(); }

    //! \brief Prints the spline information to a stream.
    void info( ostream_type & stream ) const { m_spline_2D->info( stream ); }

    //! \brief Dumps all spline data (grid and coefficients) to a stream.
    void dump_data( ostream_type & stream ) const { m_spline_2D->dump_data( stream ); }

#ifdef SPLINES_BACK_COMPATIBILITY
    // Legacy method names (camelCase) for backward compatibility.
    integer   numPointX() const { return m_spline_2D->num_point_x(); }
    integer   numPointY() const { return m_spline_2D->num_point_y(); }
    real_type xNode( integer i ) const { return this->x_node( i ); }
    real_type yNode( integer i ) const { return this->y_node( i ); }
    real_type zNode( integer i, integer j ) const { return this->z_node( i, j ); }
    real_type xMin() const { return this->x_min(); }
    real_type xMax() const { return this->x_max(); }
    real_type yMin() const { return this->y_min(); }
    real_type yMax() const { return this->y_max(); }
    real_type zMin() const { return this->z_min(); }
    real_type zMax() const { return this->z_max(); }
#endif
  };

}  // namespace Splines

// EOF Splines2D.hxx
#endif
