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
 |   ____        _ _            ____  ____
 |  / ___| _ __ | (_)_ __   ___|___ \|  _ \
 |  \___ \| '_ \| | | '_ \ / _ \ __) | | | |
 |   ___) | |_) | | | | | |  __// __/| |_| |
 |  |____/| .__/|_|_|_| |_|\___|_____|____/
 |        |_|
\*/

namespace Splines {

  //!
  //! Bi-quintic spline base class
  //!
  class Spline2D {
  protected:

    #ifndef DOXYGEN_SHOULD_SKIP_THIS

    std::string  m_name;
    SplineSurf * m_spline_2D{nullptr};

    void new_spline( SplineType2D tp );

    #endif

  public:

    //! \name Constructors
    ///@{

    //!
    //! Build an empty spline of `Spline2D` type
    //!
    //! \param name the name of the spline
    //!
    explicit
    Spline2D( string_view name = "Spline2D" )
    : m_name(name)
    {}

    //!
    //! Spline destructor.
    //!
    virtual
    ~Spline2D() {
      if ( m_spline_2D != nullptr ) {
        delete m_spline_2D;
        m_spline_2D = nullptr;
      }
    }

    ///@}

    //!
    //! \name Open/Close/Bound
    //!
    ///@{

    //!
    //! Return `true` if the surface is assumed closed in the `x` direction
    //!
    bool is_x_closed() const { return m_spline_2D->is_x_closed(); }

    //!
    //! Setup the surface as closed in the `x` direction.
    //!
    void make_x_closed() { m_spline_2D->make_x_closed(); }

    //!
    //! Setup the surface as open in the `x` direction.
    //!
    void make_x_opened() { m_spline_2D->make_x_opened(); }

    //!
    //! Return `true` if the surface is assumed closed in the `y` direction.
    //!
    bool is_y_closed() const { return m_spline_2D->is_y_closed(); }

    //!
    //! Setup the surface as closed in the `y` direction.
    //!
    void make_y_closed() { m_spline_2D->make_y_closed(); }

    //!
    //! Setup the surface as open in the `y` direction.
    //!
    void make_y_opened() { m_spline_2D->make_y_opened(); }

    //!
    //! Return `true` if the parameter `x` assumed bounded.
    //! If false the spline is estrapolated for `x` values
    //! outside the range.
    //!
    bool is_x_bounded() const { return m_spline_2D->is_x_bounded(); }

    //!
    //! Make the spline surface unbounded in the `x` direction.
    //!
    void make_x_unbounded() { m_spline_2D->make_x_unbounded(); }

    //!
    //! Make the spline surface bounded in the `x` direction.
    //!
    void make_x_bounded() { m_spline_2D->make_x_bounded(); }

    //!
    //! Return `true` if the parameter `y` assumed bounded.
    //! If false the spline is estrapolated for `y` values
    //! outside the range.
    //!
    bool is_y_bounded() const { return m_spline_2D->is_y_bounded(); }

    //!
    //! Make the spline surface unbounded in the `y` direction
    //!
    void make_y_unbounded() { m_spline_2D->make_y_unbounded(); }

    //!
    //! Make the spline surface bounded in the `y` direction
    //!
    void make_y_bounded() { m_spline_2D->make_y_bounded(); }

    ///@}

    //!
    //! \name Info
    //!
    ///@{

    //!
    //! \return string with the name of the spline
    //!
    string_view name() const { return m_spline_2D->name(); }

    //!
    //! Return the number of support points of the spline along x direction.
    //!
    integer num_point_x() const { return m_spline_2D->num_point_x(); }

    //!
    //! Return the number of support points of the spline along y direction.
    //!
    integer num_point_y() const { return m_spline_2D->num_point_y(); }

    //!
    //! Return the i-th node of the spline (x component).
    //!
    real_type x_node( integer i ) const { return m_spline_2D->x_node(i); }

    //!
    //! Return the i-th node of the spline (y component).
    //!
    real_type y_node( integer i ) const { return m_spline_2D->y_node(i); }

    //!
    //! Return the i-th node of the spline (y component).
    //!
    real_type z_node( integer i, integer j ) const { return m_spline_2D->z_node(i,j); }

    ///@}

    //!
    //! Cancel the support points, empty the spline.
    //!
    void clear() { m_spline_2D->clear(); }

    //!
    //! \name Get bounds
    //!
    ///@{

    //!
    //! Return x-minumum spline value.
    //!
    real_type x_min() const { return m_spline_2D->x_min(); }

    //!
    //! Return x-maximum spline value.
    //!
    real_type x_max() const { return m_spline_2D->x_max(); }

    //!
    //! Return y-minumum spline value.
    //!
    real_type y_min() const { return m_spline_2D->y_min(); }

    //!
    //! Return y-maximum spline value
    //!
    real_type y_max() const { return m_spline_2D->y_max(); }

    //!
    //! Return z-minumum spline value
    //!
    real_type z_min() const { return m_spline_2D->z_min(); }

    //!
    //! Return z-maximum spline value
    //!
    real_type z_max() const { return m_spline_2D->z_max(); }

    ///@}

    //!
    //! \name Constructors
    //!
    ///@{

    //!
    //! Build surface spline
    //!
    //! \param tp              spline type
    //! \param x               vector of x-coordinates
    //! \param incx            access elements as `x[0]`, `x[incx]`, `x[2*incx]`,...
    //! \param y               vector of y-coordinates
    //! \param incy            access elements as `y[0]`, `y[incx]`, `y[2*incx]`,...
    //! \param z               matrix of z-values. Elements are stored
    //!                        by row Z(i,j) = z[i*ny+j] as C-matrix
    //! \param ldZ             leading dimension of `z`
    //! \param nx              number of points in `x` direction
    //! \param ny              number of points in `y` direction
    //! \param fortran_storage if true elements are stored by column
    //!                        i.e. Z(i,j) = z[i+j*nx] as Fortran-matrix
    //! \param transposed      if true matrix Z is stored transposed
    //!
    void
    build(
      SplineType2D    tp,
      real_type const x[], integer incx,
      real_type const y[], integer incy,
      real_type const z[], integer ldZ,
      integer         nx,
      integer         ny,
      bool            fortran_storage = false,
      bool            transposed      = false
    ) {
      new_spline( tp );
      m_spline_2D->build( x, incx, y, incy, z, ldZ, nx, ny, fortran_storage, transposed );
    }

    //!
    //! Build surface spline
    //!
    //! \param tp              spline type
    //! \param x               vector of x-coordinates, nx = x.size()
    //! \param y               vector of y-coordinates, ny = y.size()
    //! \param z               matrix of z-values. Elements are stored
    //!                        by row Z(i,j) = z[i*ny+j] as C-matrix
    //! \param fortran_storage if true elements are stored by column
    //!                        i.e. Z(i,j) = z[i+j*nx] as Fortran-matrix
    //! \param transposed      if true matrix Z is stored transposed
    //!
    void
    build(
      SplineType2D              tp,
      vector<real_type> const & x,
      vector<real_type> const & y,
      vector<real_type> const & z,
      bool fortran_storage = false,
      bool transposed      = false
    ) {
      new_spline( tp );
      m_spline_2D->build( x, y, z, fortran_storage, transposed );
    }

    //!
    //! Build surface spline
    //!
    //! \param tp              spline type
    //! \param z               matrix of z-values. Elements are stored
    //!                        by row Z(i,j) = z[i*ny+j] as C-matrix
    //! \param ldZ             leading dimension of the matrix. Elements are stored
    //!                        by row Z(i,j) = z[i*ldZ+j] as C-matrix
    //! \param nx              x-dimension
    //! \param ny              y-dimension
    //! \param fortran_storage if true elements are stored by column
    //!                        i.e. Z(i,j) = z[i+j*nx] as Fortran-matrix
    //! \param transposed      if true matrix Z is stored transposed
    //!
    void
    build(
      SplineType2D    tp,
      real_type const z[],
      integer         ldZ,
      integer         nx,
      integer         ny,
      bool            fortran_storage = false,
      bool            transposed      = false
    ) {
      new_spline( tp );
      m_spline_2D->build( z, ldZ, nx, ny, fortran_storage, transposed );
    }

    //!
    //! Build surface spline
    //!
    //! \param tp              spline type
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
    void
    build(
      SplineType2D              tp,
      vector<real_type> const & z,
      integer                   nx,
      integer                   ny,
      bool fortran_storage = false,
      bool transposed      = false
    ) {
      new_spline( tp );
      m_spline_2D->build( z, nx, ny, fortran_storage, transposed );
    }

    //!
    //! Build spline surface using `gc`
    //!
    //! - gc("spline_type")
    //!     - "bilinear" build a bilinear spline surface
    //!     - "bicubic" build a spline surface with cubic spline
    //!     - "biquintic" build a spline surface with quintic spline
    //!     - "Akima" or "akima "build a spline surface with cubic spline
    //!        using Akima algorithm to avoid obscillation
    //!
    void
    setup( GenericContainer const & gc );

    //!
    //! Build a spline using data in `GenericContainer`
    //!
    void
    build( GenericContainer const & gc )
    { setup(gc); }

    ///@}

    //!
    //! \name Evaluate
    //!
    ///@{

    //!
    //! Evaluate spline value at `(x,y)`.
    //!
    real_type
    operator () ( real_type x, real_type y ) const
    { return m_spline_2D->eval( x, y ); }

    //!
    //! Evaluate spline value at `(x,y)`.
    //!
    real_type
    eval( real_type x, real_type y ) const
    { return m_spline_2D->eval( x, y ); }

    ///@}

    //! \name First derivatives:
    ///@{
    //!
    //! Value and first derivatives at point \f$ (x,y) \f$:
    //!
    //! - d[0] value of the spline \f$ S(x,y) \f$
    //! - d[1] derivative respect to \f$ x \f$ of the spline: \f$ S_x(x,y) \f$
    //! - d[2] derivative respect to \f$ y \f$ of the spline: \f$ S_y(x,y) \f$
    //!
    void
    D( real_type x, real_type y, real_type d[3] ) const
    { return m_spline_2D->D( x, y, d ); }

    //!
    //! First derivatives respect to \f$ x \f$ at point \f$ (x,y) \f$
    //! of the spline: \f$ S_x(x,y) \f$.
    //!
    real_type
    Dx( real_type x, real_type y ) const
    { return m_spline_2D->Dx( x, y ); }

    //!
    //! First derivatives respect to \f$ y \f$ at point \f$ (x,y) \f$
    //! of the spline: \f$ S_x(x,y) \f$.
    //!
    real_type
    Dy( real_type x, real_type y ) const
    { return m_spline_2D->Dy( x, y ); }

    //!
    //! First derivatives respect to \f$ x \f$ at point \f$ (x,y) \f$
    //! of the spline: \f$ S_x(x,y) \f$.
    //!
    real_type
    eval_D_1( real_type x, real_type y ) const
    { return this->Dx(x,y); }

    //!
    //! First derivatives respect to \f$ y \f$ at point \f$ (x,y) \f$
    //! of the spline: \f$ S_x(x,y) \f$.
    //!
    real_type
    eval_D_2( real_type x, real_type y ) const
    { return this->Dy(x,y); }

    ///@}

    //!
    //! \name Second derivatives:
    //!
    ///@{
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
    void
    DD( real_type x, real_type y, real_type dd[6] ) const
    { return m_spline_2D->DD( x, y, dd ); }

    //!
    //! Second derivatives respect to \f$ x \f$ at point \f$ (x,y) \f$
    //! of the spline: \f$ S_{xx}(x,y) \f$.
    //!
    real_type
    Dxx( real_type x, real_type y ) const
    { return m_spline_2D->Dxx( x, y ); }

    //!
    //! Mixed second derivatives: \f$ S_{xy}(x,y) \f$.
    //!
    real_type
    Dxy( real_type x, real_type y ) const
    { return m_spline_2D->Dxy( x, y ); }

    //!
    //! Second derivatives respect to \f$ y \f$ at point \f$ (x,y) \f$
    //! of the spline: \f$ S_{xx}(x,y) \f$.
    //!
    real_type
    Dyy( real_type x, real_type y ) const
    { return m_spline_2D->Dyy( x, y ); }

    //!
    //! Second derivatives respect to \f$ x \f$ at point \f$ (x,y) \f$
    //! of the spline: \f$ S_{xx}(x,y) \f$.
    //!
    real_type
    eval_D_1_1( real_type x, real_type y ) const
    { return this->Dxx(x,y); }

    //!
    //! Mixed second derivatives: \f$ S_{xy}(x,y) \f$.
    //!
    real_type
    eval_D_1_2( real_type x, real_type y ) const
    { return this->Dxy(x,y); }

    //!
    //! Second derivatives respect to \f$ y \f$ at point \f$ (x,y) \f$
    //! of the spline: \f$ S_{xx}(x,y) \f$.
    //!
    real_type
    eval_D_2_2( real_type x, real_type y ) const
    { return this->Dyy(x,y); }
    ///@}

    //!
    //! Print spline coefficients.
    //!
    void
    write_to_stream( ostream_type & s ) const
    { return m_spline_2D->write_to_stream( s ); }

    //!
    //! Return spline typename
    //!
    char const * type_name() const { return m_spline_2D->type_name(); }

    //!
    //! String information of the kind and order of the spline
    //!
    string
    info() const
    { return m_spline_2D->info(); }

    //!
    //! Print information of the kind and order of the spline
    //!
    void
    info( ostream_type & stream ) const
    { m_spline_2D->info( stream ); }

    //!
    //! Dump spline values on the streams
    //!
    void
    dump_data( ostream_type & stream ) const
    { m_spline_2D->dump_data( stream ); }

    #ifdef SPLINES_BACK_COMPATIBILITY
    integer numPointX() const { return m_spline_2D->num_point_x(); }
    integer numPointY() const { return m_spline_2D->num_point_y(); }
    real_type xNode( integer i ) const { return this->x_node(i); }
    real_type yNode( integer i ) const { return this->y_node(i); }
    real_type zNode( integer i, integer j ) const { return this->z_node(i,j); }
    real_type xMin() const { return this->x_min(); }
    real_type xMax() const { return this->x_max(); }
    real_type yMin() const { return this->y_min(); }
    real_type yMax() const { return this->y_max(); }
    real_type zMin() const { return this->z_min(); }
    real_type zMax() const { return this->z_max(); }
    #endif

  };

}

// EOF Splines2D.hxx
