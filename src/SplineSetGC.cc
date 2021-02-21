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

#include "Splines.hh"
#include <limits>
#include <cmath>

#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
#pragma clang diagnostic ignored "-Wpoison-system-directories"
#endif

/**
 *
 */


namespace Splines {

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  using GenericContainerNamespace::GC_INTEGER;
  using GenericContainerNamespace::GC_VEC_BOOL;
  using GenericContainerNamespace::GC_VEC_INTEGER;
  using GenericContainerNamespace::GC_VEC_REAL;
  using GenericContainerNamespace::GC_VEC_STRING;
  using GenericContainerNamespace::GC_VECTOR;
  using GenericContainerNamespace::GC_MAP;
  using GenericContainerNamespace::GC_MAT_REAL;
  using GenericContainerNamespace::mat_real_type;
  using GenericContainerNamespace::vec_int_type;
  using GenericContainerNamespace::vec_real_type;
  using GenericContainerNamespace::vec_string_type;
  using GenericContainerNamespace::vector_type;
  using GenericContainerNamespace::map_type;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineSet::setup( GenericContainer const & gc ) {
    /*
    // gc["spline_type"]
    // gc["xdata"]
    // gc["ydata"]
    // gc["headers"] (opzionale)
    */
    vec_string_type       spline_type_vec;
    vec_real_type         X;
    vector<SplineType1D>  stype;
    vec_string_type       headers;
    vector<vec_real_type> Y, Yp;

    string msg = fmt::format( "SplineSet[{}]::build(): ", m_name );

    UTILS_ASSERT(
      gc.exists("spline_type"), "{}, missing `spline_type` field!\n", msg
    );
    gc("spline_type").copyto_vec_string(
      spline_type_vec, (msg+"in reading `spline_type'\n").c_str()
    );
    m_nspl = integer(spline_type_vec.size());
    stype.resize( size_t(m_nspl) );
    for ( size_t spl = 0; spl < size_t(m_nspl); ++spl )
      stype[spl] = string_to_splineType( spline_type_vec[spl] );

    UTILS_ASSERT(
      gc.exists("xdata"), "{}, missing `xdata` field!\n", msg
    );
    gc("xdata").copyto_vec_real( X, (msg+"reading `xdata'").c_str() );
    m_npts = integer( X.size() );

    UTILS_ASSERT(
      gc.exists("ydata"), "{}, missing `ydata` field!\n", msg
    );
    GenericContainer const & gc_ydata = gc("ydata");

    // allocate for _nspl splines
    Y  . resize( size_t(m_nspl) );
    Yp . resize( size_t(m_nspl) );

    // se tipo vettore o matrice deve esserci headers
    if ( GC_MAT_REAL == gc_ydata.get_type() ||
         GC_VECTOR   == gc_ydata.get_type() ) {
      UTILS_ASSERT(
        gc.exists("headers"),
        "{}, missing `headers` field!\n", msg
      );
      GenericContainer const & gc_headers = gc("headers");
      gc_headers.copyto_vec_string(
        headers, (msg+", reading `headers'").c_str()
      );
      UTILS_ASSERT(
        headers.size() == size_t(m_nspl),
        "{}, field `headers` expected to be of size {} found of size {}\n",
        msg, m_nspl, headers.size()
      );
    }

    if ( GC_MAT_REAL == gc_ydata.get_type() ) {
      // leggo matrice
      mat_real_type const & data = gc_ydata.get_mat_real();
      UTILS_ASSERT(
        size_t(m_nspl) == data.numCols(),
        "{}, number of splines [{}]\n"
        "differs from the number of `ydata` columns [{}] in data\n",
        msg, m_nspl, data.numCols()
      );
      UTILS_ASSERT(
        size_t(m_npts) == data.numRows(),
        "{}, number of points [{}]\n"
        "differs from the number of `ydata` rows [{}] in data\n",
        msg, m_npts, data.numRows()
      );
      for ( size_t i = 0; i < size_t(m_nspl); ++i )
        data.getColumn(unsigned(i),Y[i]);
    } else if ( GC_VECTOR == gc_ydata.get_type() ) {
      vector_type const & data = gc_ydata.get_vector();
      UTILS_ASSERT(
        size_t(m_nspl) == data.size(),
        "{}, field `ydata` expected of size {} found of size {}\n",
        msg, m_nspl, data.size()
      );
      string msg1 = msg+" reading `ydata` columns";
      for ( size_t spl = 0; spl < size_t(m_nspl); ++spl ) {
        GenericContainer const & datai = data[spl];
        integer nrow = m_npts;
        if ( stype[spl] == CONSTANT_TYPE ) --nrow; // constant spline uses n-1 points
        datai.copyto_vec_real( Y[spl], msg1.c_str() );
        UTILS_ASSERT(
          size_t(m_npts) == Y[spl].size(),
          "{}, column {} of `ydata` expected of size {} found of size {}\n",
          msg, spl, m_npts, Y[spl].size()
        );
      }
    } else if ( GC_MAP == gc_ydata.get_type() ) {
      map_type const & data = gc_ydata.get_map();
      UTILS_ASSERT(
        data.size() == size_t(m_nspl),
        "{}, field `ydata` expected of size {} found of size {}\n",
        msg, m_nspl, data.size()
      );
      headers.clear(); headers.reserve(data.size());
      map_type::const_iterator im = data.begin();
      string msg1 = msg+" reading `ydata` columns";
      for ( size_t spl = 0; im != data.end(); ++im, ++spl ) {
        headers.push_back(im->first);
        GenericContainer const & datai = im->second;
        integer nrow = m_npts;
        if ( stype[spl] == CONSTANT_TYPE ) --nrow; // constant spline uses n-1 points
        datai.copyto_vec_real( Y[spl], msg1.c_str() );
        UTILS_ASSERT(
          size_t(m_npts) == Y[spl].size(),
          "{}, column `{}` of `ydata` expected of size {} found of size {}\n",
          msg, im->first, m_npts, Y[spl].size()
        );
      }
    } else {
      UTILS_ERROR(
        "{}, field `data` expected\n"
        "to be of type `mat_real_type`, `vector_type` or `map_type'\n"
        "found: `{}`\n",
        msg, gc_ydata.get_type_name()
      );
    }

    if ( gc.exists("ypdata") ) { // yp puo' essere solo tipo map
      GenericContainer const & gc_ypdata = gc("ypdata");
      UTILS_ASSERT(
        GC_MAP == gc_ypdata.get_type(),
        "{}, field `ypdata` expected to be of type `map_type` found: `{}`\n",
        msg, gc_ypdata.get_type_name()
      );

      std::map<string,integer> h_to_pos;
      vec_string_type::const_iterator is = headers.begin();
      for ( integer idx = 0; is != headers.end(); ++is )
        h_to_pos[*is] = idx++;

      string msg1 = msg+" reading `ypdata` columns";
      map_type const & data = gc_ypdata.get_map();
      map_type::const_iterator im = data.begin();
      for (; im != data.end(); ++im ) {
        // cerca posizione
        std::map<string,integer>::iterator is_pos = h_to_pos.find(im->first);
        UTILS_ASSERT(
          is_pos != h_to_pos.end(),
          "{}, column `{}` of `ypdata` not found\n",
          msg, im->first
        );
        integer spl = is_pos->second;

        GenericContainer const & datai = im->second;
        integer nrow = m_npts;
        if ( stype[size_t(spl)] == CONSTANT_TYPE ) --nrow; // constant spline uses n-1 points
        datai.copyto_vec_real( Yp[size_t(spl)], msg1.c_str() );
        UTILS_ASSERT(
          size_t(m_npts) == Y[spl].size(),
          "{}, column `{}` of `ypdata` expected of size {} found of size {}\n",
          msg, im->first, m_npts, Y[spl].size()
        );
      }
    }

    Utils::Malloc<void*> mem( msg );
    mem.allocate( size_t( 3*m_nspl ) );

    void ** __headers = mem( size_t( m_nspl ) );
    void ** __Y       = mem( size_t( m_nspl ) );
    void ** __Yp      = mem( size_t( m_nspl ) );

    for ( size_t spl = 0; spl < size_t(m_nspl); ++spl ) {
      __headers[spl] = const_cast<void *>( reinterpret_cast<void const *>(
        headers[spl].c_str()
      ) );
      __Y[spl] = reinterpret_cast<void *>(
        &Y[spl].front()
      );
      __Yp[spl] = reinterpret_cast<void *>(
        Yp[spl].size() > 0 ? &Yp[spl].front() : nullptr
      );
    }

    this->build(
      m_nspl, m_npts,
      reinterpret_cast<char const **>(const_cast<void const **>(__headers)),
      &stype.front(), &X.front(),
      reinterpret_cast<real_type const **>(const_cast<void const **>(__Y)),
      reinterpret_cast<real_type const **>(const_cast<void const **>(__Yp))
    );

    if ( gc.exists("boundary") ) {
      GenericContainer const & gc_boundary = gc("boundary");
      integer ne = integer(gc_boundary.get_num_elements());
      UTILS_ASSERT(
        ne == m_nspl,
        "{}, field `boundary` expected a"
        " generic vector of size: {} but is of size: {}\n",
        msg, ne, m_nspl
      );

      for ( integer ispl = 0; ispl < ne; ++ispl ) {
        Spline * S = m_splines[size_t(ispl)];
        GenericContainer const & item = gc_boundary(ispl,"SplineSet boundary data");

        if ( item.exists("closed") && item("closed").get_bool() ) {
          S->make_closed();
        } else {
          S->make_opened();
          if ( item.exists("extend") && item("extend").get_bool() ) {
            S->make_unbounded();
            if ( item("extend_constant").get_bool() )
              S->make_extended_constant();
            else
              S->make_extended_not_constant();
          } else {
            S->make_bounded();
          }
        }
      }
    }
  }

  /*!
   | Evaluate all the splines at `x` and fill
   | a map of values in a GenericContainer
  \*/
  void
  SplineSet::eval( real_type x, GenericContainer & gc ) const {
    map_type & vals = gc.set_map();
    for ( integer pos = 0; pos < m_header_to_position.n_elem(); ++pos ) {
      BinarySearch::DATA_TYPE const & D = m_header_to_position.get_elem( pos );
      vals[D.first] = m_splines[size_t(D.second)]->eval(x);
    }
  }

  /*!
   | Evaluate all the splines at `x` values contained in vec
   | and fill a map of vector in a GenericContainer
  \*/
  void
  SplineSet::eval( vec_real_type const & vec, GenericContainer & gc ) const {
    integer npts = integer(vec.size());
    map_type & vals = gc.set_map();
    for ( integer pos = 0; pos < m_header_to_position.n_elem(); ++pos ) {
      BinarySearch::DATA_TYPE const & D = m_header_to_position.get_elem( pos );
      vec_real_type & v = vals[D.first].set_vec_real(unsigned(npts));
      Spline const * p_spl = m_splines[size_t(D.second)];
      for ( size_t i = 0; i < size_t(npts); ++i ) v[i] = p_spl->eval(vec[i]);
    }
  }

  /*!
   | Evaluate all the splines at `x` and fill a map of values
   | in a GenericContainer with keys in `columns`
  \*/
  void
  SplineSet::eval(
    real_type               x,
    vec_string_type const & columns,
    GenericContainer      & gc
  ) const {
    map_type & vals = gc.set_map();
    for ( vec_string_type::const_iterator is = columns.begin();
          is != columns.end(); ++is ) {
      Spline const * p_spl = getSpline( is->c_str() );
      vals[*is] = p_spl->eval(x);
    }
  }

  /*!
   | Evaluate all the splines at `x` values contained in vec
   | and fill a map of vector in a GenericContainer with keys in `columns`
  \*/
  void
  SplineSet::eval(
    vec_real_type   const & vec,
    vec_string_type const & columns,
    GenericContainer      & gc
  ) const {
    integer npts = integer(vec.size());
    map_type & vals = gc.set_map();
    for ( vec_string_type::const_iterator is = columns.begin();
          is != columns.end(); ++is ) {
      Spline const * p_spl = getSpline( is->c_str() );
      vec_real_type & v = vals[*is].set_vec_real(unsigned(npts));
      for ( size_t i = 0; i < size_t(npts); ++i ) v[i] = p_spl->eval(vec[i]);
    }
  }

  /*!
   | Evaluate all the splines at `zeta` and fill a map of values
   | in a GenericContainer and `indep` as independent spline
  \*/
  void
  SplineSet::eval2(
    real_type          zeta,
    integer            indep,
    GenericContainer & gc
  ) const {
    map_type & vals = gc.set_map();
    real_type x;
    intersect( indep, zeta, x );
    for ( integer pos = 0; pos < m_header_to_position.n_elem(); ++pos ) {
      BinarySearch::DATA_TYPE const & D = m_header_to_position.get_elem( pos );
      vals[D.first] = m_splines[size_t(D.second)]->eval(x);
    }
  }

  /*!
   | Evaluate all the splines at `zeta` values contained in vec
   | and fill a map of vector in a GenericContainer and `indep`
   | as independent spline
  \*/
  void
  SplineSet::eval2(
    vec_real_type const & zetas,
    integer               indep,
    GenericContainer    & gc
  ) const {
    integer npts = integer(zetas.size());
    map_type & vals = gc.set_map();

    // preallocation
    for ( integer pos = 0; pos < m_header_to_position.n_elem(); ++pos ) {
      BinarySearch::DATA_TYPE const & D = m_header_to_position.get_elem( pos );
      vals[D.first].set_vec_real(unsigned(npts));
    }

    for ( size_t i = 0; i < size_t(npts); ++i ) {
      real_type x;
      intersect( indep, zetas[i], x );
      for ( integer pos = 0; pos < m_header_to_position.n_elem(); ++pos ) {
        BinarySearch::DATA_TYPE const & D = m_header_to_position.get_elem( pos );
        vec_real_type & v = vals[D.first].get_vec_real();
        v[i] = m_splines[size_t(D.second)]->eval(x);
      }
    }
  }

  /*!
   | Evaluate all the splines at `zeta` and fill a map of values
   | in a GenericContainer with keys in `columns` and `indep`
   | as independent spline
  \*/
  void
  SplineSet::eval2(
    real_type               zeta,
    integer                 indep,
    vec_string_type const & columns,
    GenericContainer      & gc
  ) const {
    map_type & vals = gc.set_map();
    real_type x;
    intersect( indep, zeta, x );
    for ( vec_string_type::const_iterator is = columns.begin();
          is != columns.end(); ++is ) {
      Spline const * p_spl = getSpline( is->c_str() );
      vals[*is] = p_spl->eval(x);
    }
  }

  /*!
   | Evaluate all the splines at `zeta` values contained
   | in vec and fill a map of vector in a GenericContainer
   | with keys in `columns` and `indep` as independent spline
  \*/
  void
  SplineSet::eval2(
    vec_real_type   const & zetas,
    integer                 indep,
    vec_string_type const & columns,
    GenericContainer      & gc
  ) const {
    integer npts = integer(zetas.size());
    map_type & vals = gc.set_map();

    // preallocation
    for ( vec_string_type::const_iterator is = columns.begin();
          is != columns.end(); ++is )
      vals[*is].set_vec_real(unsigned(npts));

    for ( size_t i = 0; i < size_t(npts); ++i ) {
      real_type x;
      intersect( indep, zetas[i], x );
      for ( vec_string_type::const_iterator is = columns.begin();
            is != columns.end(); ++is ) {
        vec_real_type & v = vals[*is].get_vec_real();
        Spline const * p_spl = getSpline( is->c_str() );
        v[i] = p_spl->eval(x);
      }
    }
  }

  /*
  //    __ _          _         _           _            _   _
  //   / _(_)_ __ ___| |_    __| | ___ _ __(_)_   ____ _| |_(_)_   _____
  //  | |_| | '__/ __| __|  / _` |/ _ \ '__| \ \ / / _` | __| \ \ / / _ \
  //  |  _| | |  \__ \ |_  | (_| |  __/ |  | |\ V / (_| | |_| |\ V /  __/
  //  |_| |_|_|  |___/\__|  \__,_|\___|_|  |_| \_/ \__,_|\__|_| \_/ \___|
  */
  void
  SplineSet::eval_D( real_type x, GenericContainer & gc ) const {
    map_type & vals = gc.set_map();
    for ( integer pos = 0; pos < m_header_to_position.n_elem(); ++pos ) {
      BinarySearch::DATA_TYPE const & D = m_header_to_position.get_elem( pos );
      vals[D.first] = m_splines[size_t(D.second)]->eval_D(x);
    }
  }

  /*!
   | Evaluate all the splines at `x` values contained in vec and
   | fill a map of vector in a GenericContainer
  \*/
  void
  SplineSet::eval_D( vec_real_type const & vec, GenericContainer & gc ) const {
    integer npts = integer(vec.size());
    map_type & vals = gc.set_map();
    for ( integer pos = 0; pos < m_header_to_position.n_elem(); ++pos ) {
      BinarySearch::DATA_TYPE const & D = m_header_to_position.get_elem( pos );
      vec_real_type & v = vals[D.first].set_vec_real(unsigned(npts));
      Spline const * p_spl = m_splines[size_t(D.second)];
      for ( size_t i = 0; i < size_t(npts); ++i ) v[i] = p_spl->eval_D(vec[i]);
    }
  }

  /*!
   | Evaluate all the splines at `x` and fill a map of values in
   | a GenericContainer with keys in `columns`
  \*/
  void
  SplineSet::eval_D(
    real_type               x,
    vec_string_type const & columns,
    GenericContainer      & gc
  ) const {
    map_type & vals = gc.set_map();
    for ( vec_string_type::const_iterator is = columns.begin();
          is != columns.end(); ++is ) {
      Spline const * p_spl = getSpline( is->c_str() );
      vals[*is] = p_spl->eval_D(x);
    }
  }

  /*!
   | Evaluate all the splines at `x` values contained in vec
   | and fill a map of vector in a GenericContainer with keys in `columns`
  \*/
  void
  SplineSet::eval_D(
    vec_real_type   const & vec,
    vec_string_type const & columns,
    GenericContainer      & gc
  ) const {
    integer npts = integer(vec.size());
    map_type & vals = gc.set_map();
    for ( vec_string_type::const_iterator is = columns.begin();
          is != columns.end(); ++is ) {
      Spline const * p_spl = getSpline( is->c_str() );
      vec_real_type & v = vals[*is].set_vec_real(unsigned(npts));
      for ( size_t i = 0; i < size_t(npts); ++i ) v[i] = p_spl->eval_D(vec[i]);
    }
  }

  /*!
   | Evaluate all the splines at `zeta` and fill a map of values
   | in a GenericContainer and `indep` as independent spline
  \*/
  void
  SplineSet::eval2_D(
    real_type          zeta,
    integer            indep,
    GenericContainer & gc
  ) const {
    map_type & vals = gc.set_map();
    real_type x;
    intersect( indep, zeta, x );
    for ( integer pos = 0; pos < m_header_to_position.n_elem(); ++pos ) {
      BinarySearch::DATA_TYPE const & D = m_header_to_position.get_elem( pos );
      vals[D.first] = m_splines[size_t(D.second)]->eval_D(x);
    }
  }

  /*!
   | Evaluate all the splines at `zeta` values contained in vec
   | and fill a map of vector in a GenericContainer and `indep`
   | as independent spline
  \*/
  void
  SplineSet::eval2_D(
    vec_real_type const & zetas,
    integer               indep,
    GenericContainer    & gc
  ) const {
    integer npts = integer(zetas.size());
    map_type & vals = gc.set_map();

    // preallocation
    for ( integer pos = 0; pos < m_header_to_position.n_elem(); ++pos ) {
      BinarySearch::DATA_TYPE const & D = m_header_to_position.get_elem( pos );
      vals[D.first].set_vec_real(unsigned(npts));
    }

    for ( size_t i = 0; i < size_t(npts); ++i ) {
      real_type x;
      intersect( indep, zetas[i], x );
      for ( integer pos = 0; pos < m_header_to_position.n_elem(); ++pos ) {
        BinarySearch::DATA_TYPE const & D = m_header_to_position.get_elem( pos );
        vec_real_type & v = vals[D.first].get_vec_real();
        v[i] = m_splines[size_t(D.second)]->eval_D(x);
      }
    }
  }

  /*!
   | Evaluate all the splines at `zeta` and fill a map of values
   | in a GenericContainer with keys in `columns` and `indep`
   | as independent spline
  \*/
  void
  SplineSet::eval2_D(
    real_type               zeta,
    integer                 indep,
    vec_string_type const & columns,
    GenericContainer      & gc
  ) const {
    map_type & vals = gc.set_map();
    real_type x;
    intersect( indep, zeta, x );
    for ( vec_string_type::const_iterator is = columns.begin();
          is != columns.end(); ++is ) {
      Spline const * p_spl = getSpline( is->c_str() );
      vals[*is] = p_spl->eval_D(x);
    }
  }

  /*!
   | Evaluate all the splines at `zeta` values contained in vec
   | and fill a map of vector in a GenericContainer with keys
   | in `columns` and `indep` as independent spline
  \*/
  void
  SplineSet::eval2_D(
    vec_real_type   const & zetas,
    integer                 indep,
    vec_string_type const & columns,
    GenericContainer      & gc
  ) const {
    integer npts = integer(zetas.size());
    map_type & vals = gc.set_map();

    // preallocation
    for ( vec_string_type::const_iterator is = columns.begin();
          is != columns.end(); ++is )
      vals[*is].set_vec_real(unsigned(npts));

    for ( size_t i = 0; i < size_t(npts); ++i ) {
      real_type x;
      intersect( indep, zetas[i], x );
      for ( vec_string_type::const_iterator is = columns.begin();
            is != columns.end(); ++is ) {
        vec_real_type & v = vals[*is].get_vec_real();
        Spline const * p_spl = getSpline( is->c_str() );
        v[i] = p_spl->eval_D(x);
      }
    }
  }

  /*                                _       _           _            _   _
  //   ___  ___  ___ ___  _ __   __| |   __| | ___ _ __(_)_   ____ _| |_(_)_   _____
  //  / __|/ _ \/ __/ _ \| '_ \ / _` |  / _` |/ _ \ '__| \ \ / / _` | __| \ \ / / _ \
  //  \__ \  __/ (_| (_) | | | | (_| | | (_| |  __/ |  | |\ V / (_| | |_| |\ V /  __/
  //  |___/\___|\___\___/|_| |_|\__,_|  \__,_|\___|_|  |_| \_/ \__,_|\__|_| \_/ \___|
  */
  void
  SplineSet::eval_DD( real_type x, GenericContainer & gc ) const {
    map_type & vals = gc.set_map();
    for ( integer pos = 0; pos < m_header_to_position.n_elem(); ++pos ) {
      BinarySearch::DATA_TYPE const & D = m_header_to_position.get_elem( pos );
      vals[D.first] = m_splines[size_t(D.second)]->eval_DD(x);
    }
  }

  /*!
   | Evaluate all the splines at `x` values contained in vec
   | and fill a map of vector in a GenericContainer
  \*/
  void
  SplineSet::eval_DD( vec_real_type const & vec, GenericContainer & gc ) const {
    integer npts = integer(vec.size());
    map_type & vals = gc.set_map();
    for ( integer pos = 0; pos < m_header_to_position.n_elem(); ++pos ) {
      BinarySearch::DATA_TYPE const & D = m_header_to_position.get_elem( pos );
      vec_real_type & v = vals[D.first].set_vec_real(unsigned(npts));
      Spline const * p_spl = m_splines[size_t(D.second)];
      for ( size_t i = 0; i < size_t(npts); ++i ) v[i] = p_spl->eval_DD(vec[i]);
    }
  }

  /*!
   | Evaluate all the splines at `x` and fill a map of values in
   | a GenericContainer with keys in `columns`
  \*/
  void
  SplineSet::eval_DD(
    real_type               x,
    vec_string_type const & columns,
    GenericContainer      & gc
  ) const {
    map_type & vals = gc.set_map();
    for ( vec_string_type::const_iterator is = columns.begin();
          is != columns.end(); ++is ) {
      Spline const * p_spl = getSpline( is->c_str() );
      vals[*is] = p_spl->eval_DD(x);
    }
  }

  /*!
   | Evaluate all the splines at `x` values contained in vec and
   | fill a map of vector in a GenericContainer with keys in `columns`
  \*/
  void
  SplineSet::eval_DD(
    vec_real_type   const & vec,
    vec_string_type const & columns,
    GenericContainer      & gc
  ) const {
    integer npts = integer(vec.size());
    map_type & vals = gc.set_map();
    for ( vec_string_type::const_iterator is = columns.begin();
          is != columns.end(); ++is ) {
      Spline const * p_spl = getSpline( is->c_str() );
      vec_real_type & v = vals[*is].set_vec_real(unsigned(npts));
      for ( size_t i = 0; i < size_t(npts); ++i ) v[i] = p_spl->eval_DD(vec[i]);
    }
  }

  /*!
   | Evaluate all the splines at `zeta` and fill a map of values
   | in a GenericContainer and `indep` as independent spline
  \*/
  void
  SplineSet::eval2_DD(
    real_type          zeta,
    integer            indep,
    GenericContainer & gc
  ) const {
    map_type & vals = gc.set_map();
    real_type x;
    intersect( indep, zeta, x );
    for ( integer pos = 0; pos < m_header_to_position.n_elem(); ++pos ) {
      BinarySearch::DATA_TYPE const & D = m_header_to_position.get_elem( pos );
      vals[D.first] = m_splines[size_t(D.second)]->eval_DD(x);
    }
  }

  /*!
   | Evaluate all the splines at `zeta` values contained in vec
   | and fill a map of vector in a GenericContainer and `indep`
   | as independent spline
  \*/
  void
  SplineSet::eval2_DD(
    vec_real_type const & zetas,
    integer               indep,
    GenericContainer    & gc
  ) const {
    integer npts = integer(zetas.size());
    map_type & vals = gc.set_map();

    // preallocation
    for ( integer pos = 0; pos < m_header_to_position.n_elem(); ++pos ) {
      BinarySearch::DATA_TYPE const & D = m_header_to_position.get_elem( pos );
      vals[D.first].set_vec_real(unsigned(npts));
    }

    for ( size_t i = 0; i < size_t(npts); ++i ) {
      real_type x;
      intersect( indep, zetas[i], x );
      for ( integer pos = 0; pos < m_header_to_position.n_elem(); ++pos ) {
        BinarySearch::DATA_TYPE const & D = m_header_to_position.get_elem( pos );
        vec_real_type & v = vals[D.first].get_vec_real();
        v[i] = m_splines[size_t(D.second)]->eval_DD(x);
      }
    }
  }

  /*!
   | Evaluate all the splines at `zeta` and fill a map of values
   | in a GenericContainer with keys in `columns` and `indep`
   | as independent spline
  \*/
  void
  SplineSet::eval2_DD(
    real_type               zeta,
    integer                 indep,
    vec_string_type const & columns,
    GenericContainer      & gc
  ) const {
    map_type & vals = gc.set_map();
    real_type x;
    intersect( indep, zeta, x );
    for ( vec_string_type::const_iterator is = columns.begin();
          is != columns.end(); ++is ) {
      Spline const * p_spl = getSpline( is->c_str() );
      vals[*is] = p_spl->eval_DD(x);
    }
  }

  /*!
   | Evaluate all the splines at `zeta` values contained in vec
   | and fill a map of vector in a GenericContainer with keys
   | in `columns` and `indep` as independent spline
  \*/
  void
  SplineSet::eval2_DD(
    vec_real_type   const & zetas,
    integer                 indep,
    vec_string_type const & columns,
    GenericContainer      & gc
  ) const {
    integer npts = integer(zetas.size());
    map_type & vals = gc.set_map();

    // preallocation
    for ( vec_string_type::const_iterator is = columns.begin();
          is != columns.end(); ++is )
      vals[*is].set_vec_real(unsigned(npts));

    for ( size_t i = 0; i < size_t(npts); ++i ) {
      real_type x;
      intersect( indep, zetas[i], x );
      for ( vec_string_type::const_iterator is = columns.begin();
            is != columns.end(); ++is ) {
        vec_real_type & v = vals[*is].get_vec_real();
        Spline const * p_spl = getSpline( is->c_str() );
        v[i] = p_spl->eval_DD(x);
      }
    }
  }

  /*
  //   _   _     _         _       _           _            _   _
  //  | |_| |__ (_)_ __ __| |   __| | ___ _ __(_)_   ____ _| |_(_)_   _____
  //  | __| '_ \| | '__/ _` |  / _` |/ _ \ '__| \ \ / / _` | __| \ \ / / _ \
  //  | |_| | | | | | | (_| | | (_| |  __/ |  | |\ V / (_| | |_| |\ V /  __/
  //   \__|_| |_|_|_|  \__,_|  \__,_|\___|_|  |_| \_/ \__,_|\__|_| \_/ \___|
  */
  void
  SplineSet::eval_DDD( real_type x, GenericContainer & gc ) const {
    map_type & vals = gc.set_map();
    for ( integer pos = 0; pos < m_header_to_position.n_elem(); ++pos ) {
      BinarySearch::DATA_TYPE const & D = m_header_to_position.get_elem( pos );
      vals[D.first] = m_splines[size_t(D.second)]->eval_DDD(x);
    }
  }

  /*!
   | Evaluate all the splines at `x` values contained in vec
   | and fill a map of vector in a GenericContainer
  \*/
  void
  SplineSet::eval_DDD(
    vec_real_type const & vec,
    GenericContainer    & gc
  ) const {
    integer npts = integer(vec.size());
    map_type & vals = gc.set_map();
    for ( integer pos = 0; pos < m_header_to_position.n_elem(); ++pos ) {
      BinarySearch::DATA_TYPE const & D = m_header_to_position.get_elem( pos );
      vec_real_type & v = vals[D.first].set_vec_real(unsigned(npts));
      Spline const * p_spl = m_splines[size_t(D.second)];
      for ( size_t i = 0; i < size_t(npts); ++i ) v[i] = p_spl->eval_DDD(vec[i]);
    }
  }

  /*!
   | Evaluate all the splines at `x` and fill a map of values
   | in a GenericContainer with keys in `columns`
  \*/
  void
  SplineSet::eval_DDD(
    real_type               x,
    vec_string_type const & columns,
    GenericContainer      & gc
  ) const {
    map_type & vals = gc.set_map();
    for ( vec_string_type::const_iterator is = columns.begin();
          is != columns.end(); ++is ) {
      Spline const * p_spl = getSpline( is->c_str() );
      vals[*is] = p_spl->eval_DDD(x);
    }
  }

  /*!
   | Evaluate all the splines at `x` values contained in vec and
   | fill a map of vector in a GenericContainer with keys in `columns`
  \*/
  void
  SplineSet::eval_DDD(
    vec_real_type   const & vec,
    vec_string_type const & columns,
    GenericContainer      & gc
  ) const {
    integer npts = integer(vec.size());
    map_type & vals = gc.set_map();
    for ( vec_string_type::const_iterator is = columns.begin();
          is != columns.end(); ++is ) {
      Spline const * p_spl = getSpline( is->c_str() );
      vec_real_type & v = vals[*is].set_vec_real(unsigned(npts));
      for ( size_t i = 0; i < size_t(npts); ++i ) v[i] = p_spl->eval_DDD(vec[i]);
    }
  }

  /*!
   | Evaluate all the splines at `zeta` and fill a map of values
   | in a GenericContainer and `indep` as independent spline
  \*/
  void
  SplineSet::eval2_DDD(
    real_type          zeta,
    integer            indep,
    GenericContainer & gc
  ) const {
    map_type & vals = gc.set_map();
    real_type x;
    intersect( indep, zeta, x );
    for ( integer pos = 0; pos < m_header_to_position.n_elem(); ++pos ) {
      BinarySearch::DATA_TYPE const & D = m_header_to_position.get_elem( pos );
      vals[D.first] = m_splines[size_t(D.second)]->eval_DDD(x);
    }
  }

  /*!
   | Evaluate all the splines at `zeta` values contained in vec
   | and fill a map of vector in a GenericContainer and `indep`
   | as independent spline
  \*/
  void
  SplineSet::eval2_DDD(
    vec_real_type const & zetas,
    integer               indep,
    GenericContainer    & gc
  ) const {
    integer npts = integer(zetas.size());
    map_type & vals = gc.set_map();

    // preallocation
    for ( integer pos = 0; pos < m_header_to_position.n_elem(); ++pos ) {
      BinarySearch::DATA_TYPE const & D = m_header_to_position.get_elem( pos );
      vals[D.first].set_vec_real(unsigned(npts));
    }

    for ( size_t i = 0; i < size_t(npts); ++i ) {
      real_type x;
      intersect( indep, zetas[i], x );
      for ( integer pos = 0; pos < m_header_to_position.n_elem(); ++pos ) {
        BinarySearch::DATA_TYPE const & D = m_header_to_position.get_elem( pos );
        vec_real_type & v = vals[D.first].get_vec_real();
        v[i] = m_splines[size_t(D.second)]->eval_DDD(x);
      }
    }
  }

  /*!
   | Evaluate all the splines at `zeta` and fill a map of values
   | in a GenericContainer with keys in `columns` and `indep`
   | as independent spline
  \*/
  void
  SplineSet::eval2_DDD(
    real_type               zeta,
    integer                 indep,
    vec_string_type const & columns,
    GenericContainer      & gc
  ) const {
    map_type & vals = gc.set_map();
    real_type x;
    intersect( indep, zeta, x );
    for ( vec_string_type::const_iterator is = columns.begin();
          is != columns.end(); ++is ) {
      Spline const * p_spl = getSpline( is->c_str() );
      vals[*is] = p_spl->eval_DDD(x);
    }
  }

  /*!
   | Evaluate all the splines at `zeta` values contained in vec
   | and fill a map of vector in a GenericContainer with keys
   | in `columns` and `indep` as independent spline
  \*/
  void
  SplineSet::eval2_DDD(
    vec_real_type   const & zetas,
    integer                 indep,
    vec_string_type const & columns,
    GenericContainer      & gc
  ) const {
    integer npts = integer(zetas.size());
    map_type & vals = gc.set_map();

    // preallocation
    for ( vec_string_type::const_iterator is = columns.begin();
          is != columns.end(); ++is )
      vals[*is].set_vec_real( unsigned(npts) );

    for ( size_t i = 0; i < size_t(npts); ++i ) {
      real_type x;
      intersect( indep, zetas[i], x );
      for ( vec_string_type::const_iterator is = columns.begin();
            is != columns.end(); ++is ) {
        vec_real_type & v = vals[*is].get_vec_real();
        Spline const * p_spl = getSpline( is->c_str() );
        v[i] = p_spl->eval_DDD(x);
      }
    }
  }
}
