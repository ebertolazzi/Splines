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

#include "Splines/Splines.hh"

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

  void SplineSurf::load_Z( real_type const z[], integer const ldZ, bool fortran_storage, bool transposed )
  {
    //
    //  +--------------+
    //  |              | ny = nr
    //  |              |
    //  +--------------+
    //      nx = nc
    //
    using StrideType = Eigen::OuterStride<Eigen::Dynamic>;
    integer const tf = ( transposed ? 1 : 0 ) + ( fortran_storage ? 2 : 0 );
    switch ( tf )
    {
      case 0:  // NO transpose NO fortran
      {
        UTILS_ASSERT(
          ldZ >= m_ny,
          "SplineSurf::load_Z( z, ldZ={}, fortran_storage={}, transposed={}) with nx={} and ny={} bad leading dimension",
          ldZ,
          fortran_storage,
          transposed,
          m_nx,
          m_ny );
        Eigen::Map<const MatC, 0, StrideType> ZZ( z, m_nx, m_ny, StrideType( ldZ ) );
        load_Z( ZZ, false );
      }
      break;

      case 1:  // YES transpose NO fortran
      {
        UTILS_ASSERT(
          ldZ >= m_nx,
          "SplineSurf::load_Z( z, ldZ={}, fortran_storage={}, transposed={}) with nx={} and ny={} bad leading dimension",
          ldZ,
          fortran_storage,
          transposed,
          m_nx,
          m_ny );
        Eigen::Map<const MatC, 0, StrideType> ZZ( z, m_ny, m_nx, StrideType( ldZ ) );
        load_Z( ZZ, true );
      }
      break;

      case 2:  // NO transpose YES fortran
      {
        UTILS_ASSERT(
          ldZ >= m_nx,
          "SplineSurf::load_Z( z, ldZ={}, fortran_storage={}, transposed={}) with nx={} and ny={} bad leading dimension",
          ldZ,
          fortran_storage,
          transposed,
          m_nx,
          m_ny );
        Eigen::Map<const Mat, 0, StrideType> ZZ( z, m_nx, m_ny, StrideType( ldZ ) );
        load_Z( ZZ, false );
      }
      break;

      case 3:  // YES transpose YES fortran
      {
        UTILS_ASSERT(
          ldZ >= m_ny,
          "SplineSurf::load_Z( z, ldZ={}, fortran_storage={}, transposed={}) with nx={} and ny={} bad leading dimension",
          ldZ,
          fortran_storage,
          transposed,
          m_nx,
          m_ny );
        Eigen::Map<const Mat, 0, StrideType> ZZ( z, m_ny, m_nx, StrideType( ldZ ) );
        load_Z( ZZ, true );
      }
      break;
    }
  }

  void SplineSurf::resize( integer const nx, integer const ny )
  {
    m_nx = nx;
    m_ny = ny;

    mX.resize( nx );
    mY.resize( ny );

    m_X_ptr = mX.data();
    m_Y_ptr = mY.data();

    mZ.resize( nx, ny );
  }

  void SplineSurf::clear()
  {
    m_nx = 0;
    m_ny = 0;

    m_X_ptr = nullptr;
    m_Y_ptr = nullptr;

    mX.resize( 0 );
    mY.resize( 0 );
    mZ.resize( 0, 0 );

    m_Z_min = 0;
    m_Z_max = 0;
  }

  void SplineSurf::build(
    real_type const x[],
    integer const   incx,
    real_type const y[],
    integer const   incy,
    real_type const z[],
    integer const   ldZ,
    integer const   nx,
    integer const   ny,
    bool            fortran_storage,
    bool            transposed )
  {
    resize( nx, ny );

    // -----------------------------------------------------------
    // OTTIMIZZAZIONE 1: Copia vettoriale con Stride (No cicli for)
    // -----------------------------------------------------------
    // Definiamo una "vista" sui dati di input che salta 'incx' elementi
    // Eigen userà istruzioni AVX/SSE per copiare se possibile.

    using Stride = Eigen::InnerStride<Eigen::Dynamic>;

    if ( incx == 1 )
      mX = Eigen::Map<const Vec>( x, nx );
    else
      mX = Eigen::Map<const Vec, 0, Stride>( x, nx, Stride( incx ) );

    if ( incy == 1 )
      mY = Eigen::Map<const Vec>( y, ny );
    else
      mY = Eigen::Map<const Vec, 0, Stride>( y, ny, Stride( incy ) );

    // -----------------------------------------------------------
    // OTTIMIZZAZIONE 2: Mapping intelligente di Z
    // -----------------------------------------------------------
    // Invece di passare pointer grezzi, creiamo una Map che gestisce
    // il "Leading Dimension" (ldZ) tramite OuterStride.

    integer nr = nx;
    integer nc = ny;
    if ( transposed ) std::swap( nr, nc );

    if ( fortran_storage )
    {
      // Fortran = Column Major.
      // ldZ è la distanza tra l'inizio di due colonne consecutive.
      using StrideType = Eigen::OuterStride<Eigen::Dynamic>;
      using MapType    = Eigen::Map<const Mat, 0, StrideType>;

      MapType mapZ( z, nr, nc, StrideType( ldZ ) );  // Assumo nx righe, ny colonne logiche

      // Chiama il tuo load_Z generico (che accetta Eigen::Ref)
      load_Z( mapZ, transposed );
    }
    else
    {
      // C = Row Major.
      // ldZ è la distanza tra l'inizio di due righe consecutive.
      using StrideType = Eigen::OuterStride<Eigen::Dynamic>;
      using MapType    = Eigen::Map<const MatC, 0, StrideType>;

      MapType mapZ( z, nr, nc, StrideType( ldZ ) );  // Assumo nx righe, ny colonne logiche

      // Chiama il tuo load_Z generico (che accetta Eigen::Ref)
      load_Z( mapZ, transposed );
    }

    make_spline();
  }

  void SplineSurf::build(
    real_type const z[],
    integer         ldZ,
    integer         nx,
    integer         ny,
    bool            fortran_storage,
    bool            transposed )
  {
    resize( nx, ny );
    for ( integer i = 0; i < nx; ++i ) mX.coeffRef( i ) = static_cast<real_type>( i );
    for ( integer j = 0; j < ny; ++j ) mY.coeffRef( j ) = static_cast<real_type>( j );
    load_Z( z, ldZ, fortran_storage, transposed );
    make_spline();
  }

  void SplineSurf::setup( GenericContainer const & gc )
  {
    /*
    // gc["xdata"]
    // gc["ydata"]
    // gc["zdata"]
    //
    */
    string const where = fmt::format( "SplineSurf[{}]::setup( gc ):", m_name );

    std::set<std::string> keywords;
    for ( auto const & pair : gc.get_map( where ) ) keywords.insert( pair.first );

    GenericContainer const & gc_x = gc( "xdata", where );
    keywords.erase( "xdata" );
    GenericContainer const & gc_y = gc( "ydata", where );
    keywords.erase( "ydata" );
    GenericContainer const & gc_z = gc( "zdata", where );
    keywords.erase( "zdata" );
    keywords.erase( "spline_type" );

    integer nx = static_cast<integer>( gc_x.get_num_elements() );
    integer ny = static_cast<integer>( gc_y.get_num_elements() );

    resize( nx, ny );

    for ( integer i = 0; i < m_nx; ++i ) mX.coeffRef( i ) = gc_x.get_number_at( i );
    for ( integer j = 0; j < m_ny; ++j ) mY.coeffRef( j ) = gc_y.get_number_at( j );

    bool fortran_storage = false;
    gc.get_if_exists( "fortran_storage", fortran_storage );
    keywords.erase( "fortran_storage" );
    bool transposed = false;
    gc.get_if_exists( "transposed", transposed );
    keywords.erase( "transposed" );

    /*
    //     +------+
    //  j ny      | (xi,yj)
    //     +  nx  +
    //         i
    */

    // integer const LD = fortran_storage ? NR : NC;

    auto read_mat = [this, &transposed, &fortran_storage, &where]( GenericContainer const & M ) -> void
    {
      bool          trans = transposed == fortran_storage;
      integer const NR    = trans ? m_ny : m_nx;
      integer const NC    = trans ? m_nx : m_ny;
      integer const nr    = static_cast<integer>( M.num_rows() );
      integer const nc    = static_cast<integer>( M.num_cols() );
      UTILS_ASSERT(
        NR == nr && NC == nc,
        "{}, field `zdata` is a matrix expected to be of size {} x {}, found: {} x {}\n",
        where,
        NR,
        NC,
        nr,
        nc );

      if ( GC_type::MAT_REAL == M.get_type() )
      {
        auto &                mat = M.get_mat_real();  // è in fortran storage!
        Eigen::Map<const Mat> Z( mat.data(), NR, NC );
        load_Z( Z, trans );
      }
      else
      {
        GenericContainer::mat_real_type z_tmp;
        M.copyto_mat_real( z_tmp );  // è in fortran storage!
        Eigen::Map<Mat> Z( z_tmp.data(), NR, NC );
        load_Z( Z, trans );
      }
    };

    if (
      GC_type::MAT_REAL == gc_z.get_type() || GC_type::MAT_INTEGER == gc_z.get_type() ||
      GC_type::MAT_LONG == gc_z.get_type() )
    {
      read_mat( gc_z );
    }
    else if (
      GC_type::VEC_REAL == gc_z.get_type() || GC_type::VEC_INTEGER == gc_z.get_type() ||
      GC_type::VEC_LONG == gc_z.get_type() )
    {
      // cosa mi aspetto in lettura
      integer NR  = transposed ? m_ny : m_nx;
      integer NC  = transposed ? m_nx : m_ny;
      integer nz  = static_cast<integer>( gc_z.get_num_elements() );
      integer nxy = m_nx * m_ny;
      UTILS_ASSERT(
        nz == nxy,
        "{}, field `zdata` expected to be of size {} = {}x{}, found: `{}`\n",
        where,
        nxy,
        m_nx,
        m_ny,
        nz );

      bool first_value = true;
      for ( integer k = 0; k < nxy; ++k )
      {
        integer         i, j;
        real_type const v = gc_z.get_number_at( k );
        if ( first_value )
        {
          m_Z_min     = v;
          m_Z_max     = v;
          first_value = false;
        }
        else
        {
          m_Z_min = std::min( m_Z_min, v );
          m_Z_max = std::max( m_Z_max, v );
        }
        if ( fortran_storage )
        {
          i = k % NR;
          j = k / NR;
        }
        else
        {
          j = k % NC;
          i = k / NC;
        }
        if ( transposed )
          z_node_ref( j, i ) = v;
        else
          z_node_ref( i, j ) = v;
      }
    }
    else if ( GC_type::VECTOR == gc_z.get_type() )
    {
      GenericContainer mat;
      mat.load( gc_z );
      mat.collapse();
      UTILS_ASSERT(
        GC_type::MAT_REAL == mat.get_type() || GC_type::MAT_INTEGER == mat.get_type() ||
          GC_type::MAT_LONG == mat.get_type(),
        "{}, field `zdata` cannot be converted to a matrix\n",
        where );
      read_mat( mat );
    }
    else
    {
      UTILS_ERROR(
        "{}, field `zdata` expected to be of type `mat_real_type` or  `vec_real_type` or `vector_type` found: `{}`\n",
        where,
        gc_z.get_type_name() );
    }

    UTILS_WARNING(
      keywords.empty(),
      "{}: unused keys\n{}\n",
      where,
      [&keywords]() -> string
      {
        string res;
        for ( auto const & it : keywords )
        {
          res += it;
          res += ' ';
        };
        return res;
      }() );

    make_spline();
  }

}  // namespace Splines

//
// EOF: SplineSurf.cc
//
