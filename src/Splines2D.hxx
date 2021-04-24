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

  //! Bi-quintic spline base class
  class Spline2D {
  protected:
    std::string  m_name;
    SplineSurf * m_spline_2D;

    void new_spline( SplineType2D tp );

  public:

    //! \name Constructors
    ///@{

    //!
    //! spline constructor
    //!
    Spline2D( string const & name = "Spline2D" )
    : m_name(name)
    , m_spline_2D( nullptr )
    {}

    virtual
    ~Spline2D() {
      if ( m_spline_2D != nullptr ) {
        delete m_spline_2D;
        m_spline_2D = nullptr;
      }
    }

    ///@}

    //! true if the surface is assumed closed in the `x` direction
    bool is_x_closed() const { return m_spline_2D->is_x_closed(); }
    //! setup the surface as closed in the `x` direction
    void make_x_closed() { m_spline_2D->make_x_closed(); }
    //! setup the surface as open in the `x` direction
    void make_x_opened() { m_spline_2D->make_x_opened(); }

    //! true if the surface is assumed closed in the `y` direction
    bool is_y_closed() const { return m_spline_2D->is_y_closed(); }
    //! setup the surface as closed in the `y` direction
    void make_y_closed() { m_spline_2D->make_y_closed(); }
    //! setup the surface as open in the `y` direction
    void make_y_opened() { m_spline_2D->make_y_opened(); }

    //! 
    //! return true if the parameter `x` assumed bounded.
    //! If false the spline is estrapolated for `x` values
    //! outside the range.
    //! 
    bool is_x_bounded() const { return m_spline_2D->is_x_bounded(); }
    //! make the spline surface unbounded in the `x` direction
    void make_x_unbounded() { m_spline_2D->make_x_unbounded(); }
    //! make the spline surface bounded in the `x` direction
    void make_x_bounded() { m_spline_2D->make_x_bounded(); }

    //! 
    //! return true if the parameter `y` assumed bounded.
    //! If false the spline is estrapolated for `y` values
    //! outside the range.
    //! 
    bool is_y_bounded() const { return m_spline_2D->is_y_bounded(); }
    //! make the spline surface unbounded in the `y` direction
    void make_y_unbounded() { m_spline_2D->make_y_unbounded(); }
    //! make the spline surface bounded in the `y` direction
    void make_y_bounded() { m_spline_2D->make_y_bounded(); }

    //! return the name of the spline surface assigned if was constructed
    string const & name() const { return m_spline_2D->name(); }

    //! Cancel the support points, empty the spline.
    void clear() { m_spline_2D->clear(); }

    //! return the number of support points of the spline along x direction
    integer
    numPointX() const { return m_spline_2D->numPointX(); }

    //! return the number of support points of the spline along y direction
    integer
    numPointY() const { return m_spline_2D->numPointY(); }

    //! return the i-th node of the spline (x component).
    real_type
    xNode( integer i ) const { return m_spline_2D->xNode(i); }

    //! return the i-th node of the spline (y component).
    real_type
    yNode( integer i ) const { return m_spline_2D->yNode(i); }

    //! return the i-th node of the spline (y component).
    real_type
    zNode( integer i, integer j ) const { return m_spline_2D->zNode(i,j); }

    //!
    //! \name Get bounds
    //!
    ///@{

    //! return x-minumum spline value
    real_type xMin() const { return m_spline_2D->xMin(); }

    //! return x-maximum spline value
    real_type xMax() const { return m_spline_2D->xMax(); }

    //! return y-minumum spline value
    real_type yMin() const { return m_spline_2D->yMin(); }

    //! return y-maximum spline value
    real_type yMax() const { return m_spline_2D->yMax(); }

    //! return z-minumum spline value
    real_type zMin() const { return m_spline_2D->zMin(); }

    //! return z-maximum spline value
    real_type zMax() const { return m_spline_2D->zMax(); }

    ///@}

    //!
    //! \name Constructors
    //!
    ///@{

    void
    build(
      SplineType2D      tp,
      real_type const * x, integer incx,
      real_type const * y, integer incy,
      real_type const * z, integer ldZ,
      integer           nx,
      integer           ny,
      bool              fortran_storage = false,
      bool              transposed      = false
    ) {
      new_spline( tp );
      m_spline_2D->build(
        x, incx, y, incy, z, ldZ, nx, ny, fortran_storage, transposed
      );
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
      SplineType2D      tp,
      real_type const * z,
      integer           ldZ,
      integer           nx,
      integer           ny,
      bool              fortran_storage = false,
      bool              transposed      = false
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
    //! Build spline surface using a `GenericContainer`
    //! 
    //! - gc("spline_type")
    //!     - "bilinear" build a bilinear spline surface
    //!     - "bicubic" build a spline surface with cubic spline
    //!     - "biquintic" build a spline surface with quintic spline
    //!     - "Akima" or "akima "build a spline surface with cubic spline
    //!        using Akima algorithm to avoid obscillation
    //! 
    void
    setup( GenericContainer const & gc ) {
      string msg = fmt::format("Spline2D[{}]::setup( gc ):", m_name );
      UTILS_ASSERT(
        gc.exists("spline_type"), "{}, missing `spline_type` field!\n", msg
      );
      string type = gc("spline_type").get_string();
      new_spline( string_to_splineType2D( type ) );

      m_spline_2D->setup( gc );
    }

    void
    build( GenericContainer const & gc )
    { setup(gc); }

    ///@}

    //!
    //! \name Evaluate
    //!
    ///@{

    //! Evaluate spline value
    real_type
    operator () ( real_type x, real_type y ) const
    { return (*m_spline_2D)( x, y ); }

    //! First derivative
    void
    D( real_type x, real_type y, real_type d[3] ) const
    { return m_spline_2D->D( x, y, d ); }

    real_type
    Dx( real_type x, real_type y ) const
    { return m_spline_2D->Dx( x, y ); }

    real_type
    Dy( real_type x, real_type y ) const
    { return m_spline_2D->Dy( x, y ); }

    //! Second derivative
    void
    DD( real_type x, real_type y, real_type dd[6] ) const
    { return m_spline_2D->DD( x, y, dd ); }

    real_type
    Dxx( real_type x, real_type y ) const
    { return m_spline_2D->Dxx( x, y ); }

    real_type
    Dxy( real_type x, real_type y ) const
    { return m_spline_2D->Dxy( x, y ); }

    real_type
    Dyy( real_type x, real_type y ) const
    { return m_spline_2D->Dyy( x, y ); }

    //! Evaluate spline value
    real_type
    eval( real_type x, real_type y ) const
    { return (*this)(x,y); }

    //! First derivative
    real_type
    eval_D_1( real_type x, real_type y ) const
    { return this->Dx(x,y); }

    real_type
    eval_D_2( real_type x, real_type y ) const
    { return this->Dy(x,y); }

    //! Second derivative
    real_type
    eval_D_1_1( real_type x, real_type y ) const
    { return this->Dxx(x,y); }

    real_type
    eval_D_1_2( real_type x, real_type y ) const
    { return this->Dxy(x,y); }

    real_type
    eval_D_2_2( real_type x, real_type y ) const
    { return this->Dyy(x,y); }

    ///@}

    //! Print spline coefficients
    void
    writeToStream( ostream_type & s ) const
    { return m_spline_2D->writeToStream( s ); }

    //! Return spline typename
    char const * type_name() const { return m_spline_2D->type_name(); }

    string
    info() const
    { return m_spline_2D->info(); }

    void
    info( ostream_type & stream ) const
    { m_spline_2D->info( stream ); }

    void
    dump_data( ostream_type & stream ) const
    { m_spline_2D->dump_data( stream ); }

  };

}

// EOF Splines2D.hxx
