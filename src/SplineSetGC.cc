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

    SPLINE_ASSERT(
      gc.exists("spline_type"),
      "[SplineSet[" << m_name << "]::setup] missing `spline_type` field!"
    )
    gc("spline_type").copyto_vec_string(
      spline_type_vec,
      "SplineSet::setup -- in reading `spline_type'\n"
    );
    m_nspl = integer(spline_type_vec.size());
    stype.resize( size_t(m_nspl) );
    headers.resize( size_t(m_nspl) );
    for ( size_t spl = 0; spl < size_t(m_nspl); ++spl )
      stype[spl] = string_to_splineType( spline_type_vec[spl] );

    SPLINE_ASSERT(
      gc.exists("xdata"),
      "[SplineSet[" << m_name << "]::setup] missing `xdata` field!"
    )
    gc("xdata").copyto_vec_real( X, "SplineSet::setup reading `xdata'" );
    m_npts = integer( X.size() );

    SPLINE_ASSERT(
      gc.exists("ydata"),
      "[SplineSet[" << m_name << "]::setup] missing `ydata` field!"
    )
    GenericContainer const & gc_ydata = gc("ydata");

    // allocate for _nspl splines
    Y  . resize( size_t(m_nspl) );
    Yp . resize( size_t(m_nspl) );

    // se tipo vettore o matrice deve esserci headers
    if ( GC_MAT_REAL == gc_ydata.get_type() ||
         GC_VECTOR   == gc_ydata.get_type() ) {
      SPLINE_ASSERT(
        gc.exists("headers"),
        "[SplineSet[" << m_name << "]::setup] missing `headers` field!"
      )
      GenericContainer const & gc_headers = gc("headers");
      gc_headers.copyto_vec_string(
        headers, "SplineSet::setup reading `headers'\n"
      );
      SPLINE_ASSERT(
        headers.size() == size_t(m_nspl),
        "[SplineSet[" << m_name <<
        "]::setup] field `headers` expected to be of size " << m_nspl <<
        " found of size " << headers.size()
      )
    }

    if ( GC_MAT_REAL == gc_ydata.get_type() ) {
      // leggo matrice
      mat_real_type const & data = gc_ydata.get_mat_real();
      SPLINE_ASSERT(
        size_t(m_nspl) == data.numCols(),
        "[SplineSet[" << m_name <<
        "]::setup] number of splines [" << m_nspl <<
        "] differs from the number of `ydata` columns [" <<
        data.numCols() << "] in data"
      )
      SPLINE_ASSERT(
        size_t(m_npts) == data.numRows(),
        "[SplineSet[" << m_name <<
        "]::setup] number of points [" << m_npts <<
        "] differs from the numeber of `ydata` rows [" <<
        data.numRows() << "] in data"
      )
      for ( size_t i = 0; i < size_t(m_nspl); ++i )
        data.getColumn(unsigned(i),Y[i]);
    } else if ( GC_VECTOR == gc_ydata.get_type() ) {
      vector_type const & data = gc_ydata.get_vector();
      SPLINE_ASSERT(
        size_t(m_nspl) == data.size(),
        "[SplineSet[" << m_name <<
        "]::setup] field `ydata` expected of size " << m_nspl <<
        " found of size " << data.size()
      )
      for ( size_t i = 0; i < size_t(m_nspl); ++i ) {
        GenericContainer const & datai = data[i];
        integer nrow = m_npts;
        if ( stype[i] == CONSTANT_TYPE ) --nrow; // constant spline uses n-1 points
        datai.copyto_vec_real( Y[i], "SplineSet::setup reading `ydata'" );
      }
    } else if ( GC_MAP == gc_ydata.get_type() ) {
      map_type const & data = gc_ydata.get_map();
      SPLINE_ASSERT(
        data.size() == size_t(m_nspl),
        "[SplineSet[" << m_name <<
        "]::setup] field `ydata` expected of size " << m_nspl <<
        " found of size " << data.size()
      )
      map_type::const_iterator im = data.begin();
      for ( size_t spl = 0; im != data.end(); ++im, ++spl ) {
        headers[spl] = im->first;
        GenericContainer const & datai = im->second;
        integer nrow = m_npts;
        if ( stype[spl] == CONSTANT_TYPE ) --nrow; // constant spline uses n-1 points
        datai.copyto_vec_real( Y[spl], "SplineSet::setup reading `ydata'" );
      }
    } else {
      SPLINE_DO_ERROR(
        "[SplineSet[" << m_name <<
        "]::setup] field `data` expected to be of type" <<
        " `mat_real_type`, `vector_type` or `map_type' " <<
        "found: `" << gc_ydata.get_type_name() << "`"
      )
    }

    if ( gc.exists("ypdata") ) { // yp puo essere solo tipo map
      GenericContainer const & gc_ypdata = gc("ypdata");
      SPLINE_ASSERT(
        GC_MAP == gc_ypdata.get_type(),
        "[SplineSet[" << m_name <<
        "]::setup] field `ypdata` expected to be of type "
        "`map_type` found: ` " << gc_ypdata.get_type_name() << "`"
      )
      for ( integer spl = 0; spl < m_nspl; ++spl )
        m_header_to_position.insert( headers[size_t(spl)], spl );

      map_type const & data = gc_ypdata.get_map();
      map_type::const_iterator im = data.begin();
      for (; im != data.end(); ++im ) {
        integer spl = getPosition(im->first.c_str());
        GenericContainer const & datai = im->second;
        integer nrow = m_npts;
        if ( stype[size_t(spl)] == CONSTANT_TYPE ) --nrow; // constant spline uses n-1 points
        datai.copyto_vec_real(
          Yp[size_t(spl)], "SplineSet::setup reading `ypdata'"
        );
      }
    }

    vector<char const*>      __headers;
    vector<real_type const*> __Y, __Yp;

    __headers.resize( size_t( m_nspl ) );
    __Y.resize( size_t( m_nspl ) );
    __Yp.resize( size_t( m_nspl ) );

    for ( size_t spl = 0; spl < size_t(m_nspl); ++spl ) {
      __headers[spl] = headers[spl].c_str();
      __Y[spl]       = &Y[spl].front();
      __Yp[spl]      = Yp[spl].size() > 0 ? &Yp[spl].front() : nullptr;
    }

    this->build(
      m_nspl, m_npts,
      &__headers.front(),
      &stype.front(),
      &X.front(),
      &__Y.front(),
      &__Yp.front()
    );

    if ( gc.exists("boundary") ) {
      GenericContainer const & gc_boundary = gc("boundary");
      unsigned ne = gc_boundary.get_num_elements();
      SPLINE_ASSERT(
        ne == m_splines.size(),
        "[SplineSet[" << m_name <<
        "]::setup] field `boundary` expected a generic vector of size: " << ne <<
        " but is of size: " << m_splines.size()
      )

      for ( unsigned ispl = 0; ispl < ne; ++ispl ) {
        Spline * S = m_splines[ispl];
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
