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
  class Spline2Dblend {
  protected:

    Spline2D m_surf0;
    Spline2D m_surf1;

    void check_compatibility() const;

  public:

    //! \name Constructors
    ///@{

    //!
    //! Build an empty spline of `Spline2D` type
    //!
    //! \param name the name of the spline
    //!
    explicit
    Spline2Dblend( string_view n )
    : m_surf0(fmt::format("{}_0",n))
    , m_surf1(fmt::format("{}_1",n))
    {}

    //!
    //! Spline destructor.
    //!
    virtual
    ~Spline2Dblend() {}

    ///@}

    Spline2D const & get_surf0() const { return m_surf0; }
    Spline2D const & get_surf1() const { return m_surf1; }


    //!
    //! \name Info
    //!
    ///@{

    //!
    //! Return the number of support points of the spline along x direction.
    //!
    integer num_point_x0() const { return m_surf0.num_point_x(); }
    integer num_point_x1() const { return m_surf1.num_point_x(); }

    //!
    //! Return the number of support points of the spline along y direction.
    //!
    integer num_point_y0() const { return m_surf0.num_point_y(); }
    integer num_point_y1() const { return m_surf1.num_point_y(); }

    //!
    //! \name Get bounds
    //!
    ///@{

    //!
    //! Return x-minumum spline value.
    //!
    real_type x_min() const { return m_surf0.x_min(); }

    //!
    //! Return x-maximum spline value.
    //!
    real_type x_max() const { return m_surf0.x_max(); }

    //!
    //! Return y-minumum spline value.
    //!
    real_type y_min() const { return m_surf0.y_min(); }

    //!
    //! Return y-maximum spline value
    //!
    real_type y_max() const { return m_surf0.y_max(); }

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
      SplineType2D    tp0,
      real_type const x0[], integer incx0,
      real_type const y0[], integer incy0,
      real_type const z0[], integer ldZ0,
      integer   const nx0,
      integer   const ny0,
      bool            fortran_storage0,
      bool            transposed0,
      
      SplineType2D    tp1,
      real_type const x1[], integer incx1,
      real_type const y1[], integer incy1,
      real_type const z1[], integer ldZ1,
      integer   const nx1,
      integer   const ny1,
      bool            fortran_storage1,
      bool            transposed1
    ) {
      m_surf0.build( tp0, x0, incx0, y0, incy0, z0, ldZ0, nx0, ny0, fortran_storage0, transposed0 );
      m_surf1.build( tp1, x1, incx1, y1, incy1, z1, ldZ1, nx1, ny1, fortran_storage1, transposed1 );
      check_compatibility();
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
      SplineType2D              tp0,
      vector<real_type> const & x0,
      vector<real_type> const & y0,
      vector<real_type> const & z0,
      bool                      fortran_storage0,
      bool                      transposed0,
      SplineType2D              tp1,
      vector<real_type> const & x1,
      vector<real_type> const & y1,
      vector<real_type> const & z1,
      bool                      fortran_storage1,
      bool                      transposed1
    ) {
      m_surf0.build( tp0, x0, y0, z0, fortran_storage0, transposed0 );
      m_surf1.build( tp1, x1, y1, z1, fortran_storage1, transposed1 );
      check_compatibility();
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
      SplineType2D    tp0,
      real_type const z0[],
      integer   const ldZ0,
      integer   const nx0,
      integer   const ny0,
      bool            fortran_storage0,
      bool            transposed0,
      SplineType2D    tp1,
      real_type const z1[],
      integer   const ldZ1,
      integer   const nx1,
      integer   const ny1,
      bool            fortran_storage1,
      bool            transposed1
    ) {
      m_surf0.build( tp0, z0, ldZ0, nx0, ny0, fortran_storage0, transposed0 );
      m_surf1.build( tp1, z1, ldZ1, nx1, ny1, fortran_storage1, transposed1 );
      check_compatibility();
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
      SplineType2D              tp0,
      vector<real_type> const & z0,
      integer           const   nx0,
      integer           const   ny0,
      bool                      fortran_storage0,
      bool                      transposed0,
      SplineType2D              tp1,
      vector<real_type> const & z1,
      integer           const   nx1,
      integer           const   ny1,
      bool                      fortran_storage1,
      bool                      transposed1
    ) {
      m_surf0.build( tp0, z0, nx0, ny0, fortran_storage0, transposed0 );
      m_surf1.build( tp1, z1, nx1, ny1, fortran_storage1, transposed1 );
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
    operator () ( real_type const x, real_type const y, real_type const s ) const
    { return (1-s) * m_surf0.eval( x, y ) + s * m_surf1.eval( x, y ); }

    //!
    //! Evaluate spline value at `(x,y)`.
    //!
    real_type
    eval( real_type const x, real_type const y, real_type const s ) const
    { return (1-s) * m_surf0.eval( x, y ) + s * m_surf1.eval( x, y ); }

    #ifdef AUTODIFF_SUPPORT
    autodiff::dual1st eval( autodiff::dual1st const & x, autodiff::dual1st const & y, real_type const s ) const
    { return (1-s) * m_surf0.eval( x, y ) + s * m_surf1.eval( x, y ); }
    autodiff::dual2nd eval( autodiff::dual2nd const & x, autodiff::dual2nd const & y, real_type const s ) const
    { return (1-s) * m_surf0.eval( x, y ) + s * m_surf1.eval( x, y ); }

    template <typename T1, typename T2>
    autodiff::HigherOrderDual<autodiff::detail::DualOrder<T1,T2>::value,real_type>
    eval( T1 const & x, T2 const & y, real_type const s ) const {
      autodiff::HigherOrderDual<autodiff::detail::DualOrder<T1,T2>::value,real_type> X{x}, Y{y};
      return (1-s) * m_surf0.eval( X, Y, s ) + s * m_surf1.eval( X, Y, s );
    }
    #endif

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
    D( real_type const x, real_type const y, real_type const s, real_type d[3] ) const {
      real_type d0[3], d1[3];
      m_surf0.D( x, y, d0 );
      m_surf1.D( x, y, d1 );
      d[0] = (1-s) * d0[0] + s * d1[0];
      d[1] = (1-s) * d0[1] + s * d1[1];
      d[2] = (1-s) * d0[2] + s * d1[2];
    }

    //!
    //! First derivatives respect to \f$ x \f$ at point \f$ (x,y) \f$
    //! of the spline: \f$ S_x(x,y) \f$.
    //!
    real_type
    Dx( real_type const x, real_type const y, real_type const s ) const
    { return (1-s) * m_surf0.Dx( x, y ) + s * m_surf1.Dx( x, y ); }

    //!
    //! First derivatives respect to \f$ y \f$ at point \f$ (x,y) \f$
    //! of the spline: \f$ S_x(x,y) \f$.
    //!
    real_type
    Dy( real_type const x, real_type const y, real_type const s ) const
    { return (1-s) * m_surf0.Dy( x, y ) + s * m_surf1.Dy( x, y ); }

    //!
    //! First derivatives respect to \f$ x \f$ at point \f$ (x,y) \f$
    //! of the spline: \f$ S_x(x,y) \f$.
    //!
    real_type
    eval_D_1( real_type const x, real_type const y, real_type const s ) const
    { return this->Dx(x,y,s); }

    //!
    //! First derivatives respect to \f$ y \f$ at point \f$ (x,y) \f$
    //! of the spline: \f$ S_x(x,y) \f$.
    //!
    real_type
    eval_D_2( real_type const x, real_type const y, real_type const s ) const
    { return this->Dy(x,y,s); }

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
    DD( real_type const x, real_type const y, real_type const s, real_type dd[6] ) const {
      real_type dd0[6], dd1[6];
      m_surf0.DD( x, y, dd0 );
      m_surf1.DD( x, y, dd1 );
      dd[0] = (1-s) * dd0[0] + s * dd1[0];
      dd[1] = (1-s) * dd0[1] + s * dd1[1];
      dd[2] = (1-s) * dd0[2] + s * dd1[2];
      dd[3] = (1-s) * dd0[3] + s * dd1[3];
      dd[4] = (1-s) * dd0[4] + s * dd1[4];
      dd[5] = (1-s) * dd0[5] + s * dd1[5];
    }

    //!
    //! Second derivatives respect to \f$ x \f$ at point \f$ (x,y) \f$
    //! of the spline: \f$ S_{xx}(x,y) \f$.
    //!
    real_type
    Dxx( real_type const x, real_type const y, real_type const s ) const
    { return (1-s) * m_surf0.Dxx( x, y ) + s * m_surf1.Dxx( x, y ); }

    //!
    //! Mixed second derivatives: \f$ S_{xy}(x,y) \f$.
    //!
    real_type
    Dxy( real_type const x, real_type const y, real_type const s ) const
    { return (1-s) * m_surf0.Dxy( x, y ) + s * m_surf1.Dxy( x, y ); }

    //!
    //! Second derivatives respect to \f$ y \f$ at point \f$ (x,y) \f$
    //! of the spline: \f$ S_{xx}(x,y) \f$.
    //!
    real_type
    Dyy( real_type const x, real_type const y, real_type const s ) const
    { return (1-s) * m_surf0.Dyy( x, y ) + s * m_surf1.Dyy( x, y ); }

    //!
    //! Second derivatives respect to \f$ x \f$ at point \f$ (x,y) \f$
    //! of the spline: \f$ S_{xx}(x,y) \f$.
    //!
    real_type
    eval_D_1_1( real_type const x, real_type const y, real_type const s ) const
    { return this->Dxx(x,y,s); }

    //!
    //! Mixed second derivatives: \f$ S_{xy}(x,y) \f$.
    //!
    real_type
    eval_D_1_2( real_type const x, real_type const y, real_type const s ) const
    { return this->Dxy(x,y,s); }

    //!
    //! Second derivatives respect to \f$ y \f$ at point \f$ (x,y) \f$
    //! of the spline: \f$ S_{xx}(x,y) \f$.
    //!
    real_type
    eval_D_2_2( real_type const x, real_type const y, real_type const s ) const
    { return this->Dyy(x,y,s); }
    ///@}

  };

}

// EOF Splines2D.hxx
