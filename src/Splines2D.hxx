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
    SplineSurf * m_pSpline2D;
  public:

    //! spline constructor
    Spline2D( string const & name = "Spline2D" )
    : m_name(name)
    , m_pSpline2D( nullptr )
    {}

    ~Spline2D()
    {}

    bool is_x_closed() const { return m_pSpline2D->is_x_closed(); }
    void make_x_closed()     { m_pSpline2D->make_x_closed(); }
    void make_x_opened()     { m_pSpline2D->make_x_opened(); }

    bool is_y_closed() const { return m_pSpline2D->is_y_closed(); }
    void make_y_closed()     { m_pSpline2D->make_y_closed(); }
    void make_y_opened()     { m_pSpline2D->make_y_opened(); }

    bool is_x_bounded() const { return m_pSpline2D->is_x_bounded(); }
    void make_x_unbounded()   { m_pSpline2D->make_x_unbounded(); }
    void make_x_bounded()     { m_pSpline2D->make_x_bounded(); }

    bool is_y_bounded() const { return m_pSpline2D->is_y_bounded(); }
    void make_y_unbounded()   { m_pSpline2D->make_y_unbounded(); }
    void make_y_bounded()     { m_pSpline2D->make_y_bounded(); }

    string const & name() const { return m_pSpline2D->name(); }

    //! Cancel the support points, empty the spline.
    void clear() { m_pSpline2D->clear(); }

    //! return the number of support points of the spline along x direction
    integer
    numPointX() const { return m_pSpline2D->numPointX(); }

    //! return the number of support points of the spline along y direction
    integer
    numPointY() const { return m_pSpline2D->numPointY(); }

    //! return the i-th node of the spline (x component).
    real_type
    xNode( integer i ) const { return m_pSpline2D->xNode(i); }

    //! return the i-th node of the spline (y component).
    real_type
    yNode( integer i ) const { return m_pSpline2D->yNode(i); }

    //! return the i-th node of the spline (y component).
    real_type
    zNode( integer i, integer j ) const { return m_pSpline2D->zNode(i,j); }

    //! return x-minumum spline value
    real_type
    xMin() const { return m_pSpline2D->xMin(); }

    //! return x-maximum spline value
    real_type
    xMax() const { return m_pSpline2D->xMax(); }

    //! return y-minumum spline value
    real_type
    yMin() const { return m_pSpline2D->yMin(); }

    //! return y-maximum spline value
    real_type
    yMax() const { return m_pSpline2D->yMax(); }

    //! return z-minumum spline value
    real_type
    zMin() const { return m_pSpline2D->zMin(); }

    //! return z-maximum spline value
    real_type
    zMax() const { return m_pSpline2D->zMax(); }

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
    );

    /*!
     * \brief Build surface spline
     * 
     * \param tp              spline type
     * \param x                vector of x-coordinates, nx = x.size()
     * \param y                vector of y-coordinates, ny = y.size()
     * \param z                matrix of z-values. Elements are stored
     *                         by row Z(i,j) = z[i*ny+j] as C-matrix
     * \param fortran_storage if true elements are stored by column
     *                        i.e. Z(i,j) = z[i+j*nx] as Fortran-matrix
     * \param transposed      if true matrix Z is stored transposed
     */
    void
    build(
      SplineType2D              tp,
      vector<real_type> const & x,
      vector<real_type> const & y,
      vector<real_type> const & z,
      bool fortran_storage = false,
      bool transposed      = false
    );

    /*!
     * \brief Build surface spline
     * 
     * \param tp              spline type
     * \param z               matrix of z-values. Elements are stored
     *                        by row Z(i,j) = z[i*ny+j] as C-matrix
     * \param ldZ             leading dimension of the matrix. Elements are stored
     *                        by row Z(i,j) = z[i*ldZ+j] as C-matrix
     * \param nx              x-dimension
     * \param ny              y-dimension
     * \param fortran_storage if true elements are stored by column
     *                        i.e. Z(i,j) = z[i+j*nx] as Fortran-matrix
     * \param transposed      if true matrix Z is stored transposed
     */
    void
    build(
      SplineType2D      tp,
      real_type const * z,
      integer           ldZ,
      integer           nx,
      integer           ny,
      bool              fortran_storage = false,
      bool              transposed      = false
    );

    /*!
     * \brief Build surface spline
     * 
     * \param tp              spline type
     * \param z               matrix of z-values. Elements are stored
     *                        by row Z(i,j) = z[i*ny+j] as C-matrix.
     *                        ldZ leading dimension of the matrix is ny for C-storage
     *                        and nx for Fortran storage.
     * \param nx              x-dimension
     * \param ny              y-dimension
     * \param fortran_storage if true elements are stored by column
     *                        i.e. Z(i,j) = z[i+j*nx] as Fortran-matrix
     * \param transposed      if true matrix Z is stored transposed
     */
    void
    build(
      SplineType2D              tp,
      vector<real_type> const & z,
      integer                   nx,
      integer                   ny,
      bool fortran_storage = false,
      bool transposed      = false
    );

    void
    setup( GenericContainer const & gc )
    { build(gc); }

    void
    build( GenericContainer const & gc );

    //! Evaluate spline value
    real_type
    operator () ( real_type x, real_type y ) const
    { return (*m_pSpline2D)( x, y ); }

    //! First derivative
    void
    D( real_type x, real_type y, real_type d[3] ) const
    { return m_pSpline2D->D( x, y, d ); }

    real_type
    Dx( real_type x, real_type y ) const
    { return m_pSpline2D->Dx( x, y ); }

    real_type
    Dy( real_type x, real_type y ) const
    { return m_pSpline2D->Dy( x, y ); }

    //! Second derivative
    void
    DD( real_type x, real_type y, real_type dd[6] ) const
    { return m_pSpline2D->DD( x, y, dd ); }

    real_type
    Dxx( real_type x, real_type y ) const
    { return m_pSpline2D->Dxx( x, y ); }

    real_type
    Dxy( real_type x, real_type y ) const
    { return m_pSpline2D->Dxy( x, y ); }

    real_type
    Dyy( real_type x, real_type y ) const
    { return m_pSpline2D->Dyy( x, y ); }

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

    //! Print spline coefficients
    void
    writeToStream( ostream_type & s ) const
    { return m_pSpline2D->writeToStream( s ); }

    //! Return spline typename
    char const * type_name() const { return m_pSpline2D->type_name(); }

    string
    info() const
    { return m_pSpline2D->info(); }

    void
    info( ostream_type & stream ) const
    { m_pSpline2D->info( stream ); }
  };

}

// EOF Splines2D.hxx
