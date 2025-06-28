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

#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wpoison-system-directories"
#pragma clang diagnostic ignored "-Wundefined-func-template"
#pragma clang diagnostic ignored "-Wsign-conversion"
#endif

#include "Splines.hh"
#include "Utils_fmt.hh"

#include <limits>
#include <cmath>
#include <set>

namespace Splines {

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  using GC_namespace::GC_type;
  using GC_namespace::mat_real_type;
  using GC_namespace::vec_int_type;
  using GC_namespace::vec_real_type;
  using GC_namespace::vec_string_type;
  using GC_namespace::vector_type;
  using GC_namespace::map_type;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineSet::setup( GenericContainer const & gc ) {

    string const where{ fmt::format( "SplineSet[{}]::setup( gc ): ", m_name ) };

    std::set<std::string> keywords;
    for ( auto const & pair : gc.get_map(where) ) { keywords.insert(pair.first); }

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

    GenericContainer const & gc_stype{ gc("spline_type",where) }; keywords.erase("spline_type");
    GenericContainer const & gc_xdata{ gc("xdata",where)       }; keywords.erase("xdata");
    GenericContainer const & gc_ydata{ gc("ydata",where)       }; keywords.erase("ydata");

    //
    // gc["spline_type"]
    //
    gc_stype.copyto_vec_string( spline_type_vec, (where+"in reading `spline_type'\n") );
    m_nspl = static_cast<integer>(spline_type_vec.size());
    stype.resize( m_nspl );
    for ( integer spl{0}; spl < m_nspl; ++spl )
      stype[spl] = string_to_splineType1D( spline_type_vec[spl] );

    // se non tipo MAP deve esserci headers
    //
    // gc["headers"] (opzionale)
    //
    if ( GC_type::MAP != gc_ydata.get_type() ) {
      GenericContainer const & gc_headers{ gc("headers",where) }; keywords.erase("headers");
      gc_headers.copyto_vec_string( headers, (where+", reading `headers'") );
      UTILS_ASSERT(
        headers.size() == static_cast<size_t>(m_nspl),
        "{}, field `headers` expected to be of size {} found of size {}\n",
        where, m_nspl, headers.size()
      );
    }

    //
    // gc["xdata"]
    //
    gc_xdata.copyto_vec_real( X, (where+"reading `xdata'") );
    m_npts = static_cast<integer>(X.size());

    // allocate for _nspl splines
    Y  . resize( m_nspl );
    Yp . resize( m_nspl );

    switch ( gc_ydata.get_type() ) {
    case GC_type::VEC_BOOL:
    case GC_type::VEC_INTEGER:
    case GC_type::VEC_LONG:
    case GC_type::VEC_REAL:
    case GC_type::VEC_COMPLEX:
      {
        // leggo vettore come matrice di una sola colonna
        vec_real_type data;
        gc_ydata.copyto_vec_real( data, (where+"reading `ydata'") );
        UTILS_ASSERT(
          m_nspl == 1,
          "{}, number of splines [{}]\n"
          "is incompatible with the type of `ydata` a vector (or a matrix with 1 column) in data\n",
          where, m_nspl
        );
        UTILS_ASSERT(
          static_cast<size_t>(m_npts) == data.size(),
          "{}, number of points [{}]\n"
          "differs from the number of `ydata` rows [{}] in data\n",
          where, m_npts, data.size()
        );
        Y[0].reserve(data.size());
        std::copy( data.begin(), data.end(), std::back_inserter(Y[0]) );
      }
      break;
    case GC_type::MAT_INTEGER:
    case GC_type::MAT_LONG:
    case GC_type::MAT_REAL:
      {
        mat_real_type data;
        gc_ydata.copyto_mat_real( data, (where+"reading `ydata'") );
        UTILS_ASSERT(
          static_cast<size_t>(m_nspl) == data.num_cols(),
          "{}, number of splines [{}]\n"
          "differs from the number of `ydata` columns [{}] in data\n",
          where, m_nspl, data.num_cols()
        );
        UTILS_ASSERT(
          static_cast<size_t>(m_npts) == data.num_rows(),
          "{}, number of points [{}]\n"
          "differs from the number of `ydata` rows [{}] in data\n",
          where, m_npts, data.num_rows()
        );
        for ( integer i{0}; i < m_nspl; ++i )
          data.get_column(static_cast<unsigned>(i),Y[i]);
      }
      break;
    case GC_type::VECTOR:
      {
        vector_type const & data{ gc_ydata.get_vector() };
        UTILS_ASSERT(
          static_cast<size_t>(m_nspl) == data.size(),
          "{}, field `ydata` expected of size {} found of size {}\n",
          where, m_nspl, data.size()
        );
        string msg1{ where+" reading `ydata` columns" };
        for ( integer spl = 0; spl < m_nspl; ++spl ) {
          GenericContainer const & datai{ data[spl] };
          integer nrow{ m_npts };
          if ( stype[spl] == SplineType1D::CONSTANT ) --nrow; // constant spline uses n-1 points
          datai.copyto_vec_real( Y[spl], msg1 );
          UTILS_ASSERT(
            static_cast<size_t>(nrow) == Y[spl].size(),
            "{}, column {} of `ydata` of type `{}` expected of size {} found of size {}\n",
            where, spl, spline_type_vec[spl], nrow, Y[spl].size()
          );
        }
      }
      break;
    case GC_type::MAP:
      {
        map_type const & data{ gc_ydata.get_map() };
        UTILS_ASSERT(
          data.size() == static_cast<size_t>(m_nspl),
          "{}, field `ydata` expected of size {} found of size {}\n",
          where, m_nspl, data.size()
        );
        headers.clear(); headers.reserve(data.size());
        auto im{ data.begin() };
        string msg1{ where+" reading `ydata` columns" };
        for ( integer spl{0}; im != data.end(); ++im, ++spl ) {
          headers.emplace_back(im->first);
          GenericContainer const & datai{ im->second };
          integer nrow{ m_npts };
          if ( stype[spl] == SplineType1D::CONSTANT ) --nrow; // constant spline uses n-1 points
          datai.copyto_vec_real( Y[spl], msg1 );
          UTILS_ASSERT(
            static_cast<size_t>(nrow) == Y[spl].size(),
            "{}, column `{}` of `ydata` ot type `{}` expected of size {} found of size {}\n",
            where, im->first, spline_type_vec[spl], nrow, Y[spl].size()
          );
        }
      }
      break;
    default:
      UTILS_ERROR(
        "{}, field `data` expected\n"
        "to be of type `vec_[int/long/real]_type`, `mat_[int/long/real]_type`, `vector_type` or `map_type'\n"
        "found: `{}`\n",
        where, gc_ydata.get_type_name()
      );
      break;
    }

    if ( gc.exists("ypdata") ) { // yp puo' essere solo tipo map
      GenericContainer const & gc_ypdata{ gc("ypdata",where) }; keywords.erase("ypdata");
      UTILS_ASSERT(
        GC_type::MAP == gc_ypdata.get_type(),
        "{}, field `ypdata` expected to be of type `map_type` found: `{}`\n",
        where, gc_ypdata.get_type_name()
      );

      std::map<string,integer> h_to_pos;
      vec_string_type::const_iterator is{ headers.begin() };
      for ( integer idx{0}; is != headers.end(); ++is )
        h_to_pos[*is] = idx++;

      string msg1{ where+" reading `ypdata` columns" };
      map_type const & data{gc_ypdata.get_map()};
      for ( auto const & [fst, snd] : data ) {
        // cerca posizione
        auto is_pos{ h_to_pos.find(fst) };
        UTILS_ASSERT(
          is_pos != h_to_pos.end(),
          "{}, column `{}` of `ypdata` not found\n",
          where, fst
        );
        integer spl{ is_pos->second };

        GenericContainer const & datai{snd};
        integer nrow{ m_npts };
        if ( stype[ spl ] == SplineType1D::CONSTANT ) --nrow; // constant spline uses n-1 points
        datai.copyto_vec_real( Yp[ spl ], msg1 );
        UTILS_ASSERT(
          static_cast<size_t>(nrow) == Y[spl].size(),
          "{}, column `{}` of `ypdata` or type `{}` expected of size {} found of size {}\n",
          where, fst, spline_type_vec[spl], nrow, Y[spl].size()
        );
      }
    }

    Utils::Malloc<void*> mem( where );
    mem.allocate( 3*m_nspl );

    void ** pp__headers = mem( m_nspl );
    void ** pp__Y       = mem( m_nspl );
    void ** pp__Yp      = mem( m_nspl );

    for ( integer spl{0}; spl < m_nspl; ++spl ) {
      pp__headers[spl] = const_cast<void *>( reinterpret_cast<void const *>( headers[spl].data() ) );
      pp__Y[spl]       = reinterpret_cast<void *>( Y[spl].data() );
      pp__Yp[spl]      = reinterpret_cast<void *>( Yp[spl].empty() ? nullptr : Yp[spl].data() );
    }

    this->build(
      m_nspl, m_npts,
      reinterpret_cast<char const **>(const_cast<void const **>(pp__headers)),
      stype.data(), X.data(),
      reinterpret_cast<real_type const **>(const_cast<void const **>(pp__Y)),
      reinterpret_cast<real_type const **>(const_cast<void const **>(pp__Yp))
    );

    if ( gc.exists("boundary") ) {
      GenericContainer const & gc_boundary{ gc("boundary") }; keywords.erase("boundary");
      integer ne{ static_cast<integer>(gc_boundary.get_num_elements()) };
      UTILS_ASSERT(
        ne == m_nspl,
        "{}, field `boundary` expected a"
        " generic vector of size: {} but is of size: {}\n",
        where, ne, m_nspl
      );

      for ( integer ispl{0}; ispl < ne; ++ispl ) {
        Spline * S{ m_splines[ispl].get() };
        GenericContainer const & item{ gc_boundary(ispl,"SplineSet boundary data") };

        std::set<std::string> keywords2;
        for ( auto const & pair : item.get_map(where) ) { keywords2.insert(pair.first); }

        keywords2.erase("extend");
        keywords2.erase("can_extend");
        keywords2.erase("extend_constant");
        keywords2.insert("closed");

        bool is_closed{false};
        item.get_if_exists("closed",is_closed);
        if ( is_closed ) {
          S->make_closed();
        } else {
          S->make_opened();
          bool can_extend{false};
          if ( !item.get_if_exists("extend",can_extend) ) item.get_if_exists("can_extend",can_extend) ;
          if ( can_extend ) {
            S->make_unbounded();
            bool extend_constant{false};
            item.get_if_exists("extend_constant",extend_constant);
            if ( extend_constant ) S->make_extended_constant();
            else                   S->make_extended_not_constant();
          } else {
            S->make_bounded();
          }
        }
        UTILS_WARNING(
          keywords2.empty(), "{}, spline N.{} of {} unused keys\n{}\n", where, ispl, ne, 
          [&keywords2]()->string {
            string res;
            for ( auto const & it : keywords2 ) { res += it; res += ' '; };
            return res;
          }()
        );
      }
    }
    
    UTILS_WARNING(
      keywords.empty(), "{}: unused keys\n{}\n", where,
      [&keywords]()->string {
        string res;
        for ( auto const & it : keywords ) { res += it; res += ' '; };
        return res;
      }()
    );

  }

  //!
  //! Evaluate all the splines at `x` and
  //! fill a `map of values` in a GenericContainer
  //!
  void
  SplineSet::eval( real_type const x, GenericContainer & gc ) const {
    map_type & vals{ gc.set_map() };
    for ( auto const & [fst, snd] : m_header_to_position )
      vals[fst] = m_splines[snd]->eval(x);
  }

  //!
  //! Evaluate all the splines at `x` values contained in vec
  //! and fill a `map of vector` in a GenericContainer
  //!
  void
  SplineSet::eval( vec_real_type const & vec, GenericContainer & gc ) const {
    integer const npts{ static_cast<integer>(vec.size()) };
    map_type & vals{ gc.set_map() };
    for ( auto const & [fst, snd] : m_header_to_position ) {
      vec_real_type & v{ vals[fst].set_vec_real(snd) };
      Spline const * p_spl{ m_splines[snd].get() };
      for ( integer i{0}; i < npts; ++i ) v[i] = p_spl->eval(vec[i]);
    }
  }

  //!
  //! Evaluate all the splines at `x` and fill a map of values
  //! in a GenericContainer with keys in `columns`
  //!
  void
  SplineSet::eval(
    real_type       const   x,
    vec_string_type const & columns,
    GenericContainer      & gc
  ) const {
    map_type & vals{ gc.set_map() };
    for ( auto const & S : columns ) {
      Spline const * p_spl{ get_spline(S) };
      vals[S] = p_spl->eval(x);
    }
  }

  //!
  //! Evaluate all the splines at `x` values contained in vec
  //! and fill a map of vector in a GenericContainer with keys in `columns`
  //!
  void
  SplineSet::eval(
    vec_real_type   const & vec,
    vec_string_type const & columns,
    GenericContainer      & gc
  ) const {
    integer const npts{ static_cast<integer>(vec.size()) };
    map_type & vals{ gc.set_map() };
    for ( auto const & S : columns ) {
      Spline const * p_spl{ get_spline( S ) };
      vec_real_type & v{ vals[S].set_vec_real(npts) };
      for ( integer i{0}; i < npts; ++i ) v[i] = p_spl->eval(vec[i]);
    }
  }

  //!
  //! Evaluate all the splines at `zeta` and fill a `map of values`
  //! in a GenericContainer and `indep` as independent spline
  //!
  void
  SplineSet::eval2(
    real_type const    zeta,
    integer   const    indep,
    GenericContainer & gc
  ) const {
    map_type & vals{ gc.set_map() };
    real_type x;
    intersect( indep, zeta, x );
    for ( auto const & [fst, snd] : m_header_to_position )
      vals[fst] = m_splines[snd]->eval(x);
  }

  //!
  //! Evaluate all the splines at `zeta` values contained in vec
  //! and fill a `map of vector` in a GenericContainer and `indep`
  //! as independent spline
  //!
  void
  SplineSet::eval2(
    vec_real_type const & zetas,
    integer       const   indep,
    GenericContainer    & gc
  ) const {
    integer const npts{ static_cast<integer>(zetas.size()) };
    map_type & vals{ gc.set_map() };

    // pre-allocation
    for ( auto const & [fst, snd] : m_header_to_position )
      vals[fst].set_vec_real(npts);

    for ( integer i{0}; i < npts; ++i ) {
      real_type x;
      intersect( indep, zetas[i], x );
     for ( auto const & [fst, snd] : m_header_to_position ) {
        vec_real_type & v{ vals[fst].get_vec_real() };
        v[i] = m_splines[snd]->eval(x);
      }
    }
  }

  //!
  //! Evaluate all the splines at `zeta` and fill a map of values
  //! in a GenericContainer with keys in `columns` and `indep`
  //! as independent spline
  //!
  void
  SplineSet::eval2(
    real_type       const   zeta,
    integer         const   indep,
    vec_string_type const & columns,
    GenericContainer      & gc
  ) const {
    map_type & vals = gc.set_map();
    real_type x;
    intersect( indep, zeta, x );
    for ( auto const & S : columns ) {
      Spline const * p_spl{ get_spline( S ) };
      vals[S] = p_spl->eval(x);
    }
  }

  //!
  //! Evaluate all the splines at `zeta` values contained
  //! in vec and fill a map of vector in a GenericContainer
  //! with keys in `columns` and `indep` as independent spline
  //!
  void
  SplineSet::eval2(
    vec_real_type   const & zetas,
    integer         const   indep,
    vec_string_type const & columns,
    GenericContainer      & gc
  ) const {
    integer const npts{ static_cast<integer>(zetas.size()) };
    map_type & vals{ gc.set_map() };

    // pre-allocation
    for ( auto const & S : columns ) vals[S].set_vec_real(npts);

    for ( integer i{0}; i < npts; ++i ) {
      real_type x;
      intersect( indep, zetas[i], x );
      for ( auto const & S : columns ) {
        vec_real_type & v{ vals[S].get_vec_real() };
        Spline const * p_spl{ get_spline( S ) };
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
  SplineSet::eval_D( real_type const x, GenericContainer & gc ) const {
    map_type & vals{ gc.set_map() };
    for ( auto const & [fst, snd] : m_header_to_position )
      vals[fst] = m_splines[snd]->eval_D(x);
  }

  //!
  //! Evaluate all the splines at `x` values contained in vec and
  //! fill a `map of vector` in a GenericContainer
  //!
  void
  SplineSet::eval_D( vec_real_type const & vec, GenericContainer & gc ) const {
    integer const npts{ static_cast<integer>(vec.size()) };
    map_type & vals{ gc.set_map() };
    for ( auto const & [fst, snd] : m_header_to_position ) {
      vec_real_type & v{ vals[fst].set_vec_real(npts) };
      Spline const * p_spl{ m_splines[snd].get() };
      for ( integer i{0}; i < npts; ++i ) v[i] = p_spl->eval_D(vec[i]);
    }
  }

  //!
  //! Evaluate all the splines at `x` and fill a map of values in
  //! a GenericContainer with keys in `columns`
  //!
  void
  SplineSet::eval_D(
    real_type       const   x,
    vec_string_type const & columns,
    GenericContainer      & gc
  ) const {
    map_type & vals{ gc.set_map() };
    for ( auto const & S : columns ) {
      Spline const * p_spl{ get_spline( S ) };
      vals[S] = p_spl->eval_D(x);
    }
  }

  //!
  //! Evaluate all the splines at `x` values contained in vec
  //! and fill a map of vector in a GenericContainer with keys in `columns`
  //!
  void
  SplineSet::eval_D(
    vec_real_type   const & vec,
    vec_string_type const & columns,
    GenericContainer      & gc
  ) const {
    integer const npts{ static_cast<integer>(vec.size()) };
    map_type & vals{ gc.set_map() };
    for ( auto const & S : columns ) {
      Spline const * p_spl{ get_spline( S ) };
      vec_real_type & v{ vals[S].set_vec_real(npts) };
      for ( integer i{0}; i < npts; ++i ) v[i] = p_spl->eval_D(vec[i]);
    }
  }

  //!
  //! Evaluate all the splines at `zeta` and fill a `map of values`
  //! in a GenericContainer and `indep` as independent spline
  //!
  void
  SplineSet::eval2_D(
    real_type const    zeta,
    integer   const    indep,
    GenericContainer & gc
  ) const {
    map_type & vals{ gc.set_map() };
    real_type x;
    intersect( indep, zeta, x );
    for ( auto const & [fst, snd] : m_header_to_position )
      vals[fst] = m_splines[snd]->eval_D(x);
  }

  //!
  //! Evaluate all the splines at `zeta` values contained in vec
  //! and fill a `map of vector` in a GenericContainer and `indep`
  //! as independent spline
  //!
  void
  SplineSet::eval2_D(
    vec_real_type const & zetas,
    integer       const   indep,
    GenericContainer    & gc
  ) const {
    integer const npts{ static_cast<integer>(zetas.size()) };
    map_type & vals{ gc.set_map() };

    // pre-allocation
    for ( auto const & [fst, snd] : m_header_to_position )
      vals[fst].set_vec_real(npts);

    for ( integer i{0}; i < npts; ++i ) {
      real_type x;
      intersect( indep, zetas[i], x );
      for ( auto const & [fst, snd] : m_header_to_position ) {
        vec_real_type & v{ vals[fst].get_vec_real() };
        v[i] = m_splines[snd]->eval_D(x);
      }
    }
  }

  //!
  //! Evaluate all the splines at `zeta` and fill a map of values
  //! in a GenericContainer with keys in `columns` and `indep`
  //! as independent spline
  //!
  void
  SplineSet::eval2_D(
    real_type       const   zeta,
    integer         const   indep,
    vec_string_type const & columns,
    GenericContainer      & gc
  ) const {
    map_type & vals{ gc.set_map() };
    real_type x;
    intersect( indep, zeta, x );
    for ( auto const & S : columns ) {
      Spline const * p_spl{ get_spline( S ) };
      vals[S] = p_spl->eval_D(x);
    }
  }

  //!
  //! Evaluate all the splines at `zeta` values contained in vec
  //! and fill a map of vector in a GenericContainer with keys
  //! in `columns` and `indep` as independent spline
  //!
  void
  SplineSet::eval2_D(
    vec_real_type   const & zetas,
    integer         const   indep,
    vec_string_type const & columns,
    GenericContainer      & gc
  ) const {
    integer const npts{ static_cast<integer>(zetas.size()) };
    map_type & vals{ gc.set_map() };

    // pre-allocation
    for ( auto const & S : columns ) vals[S].set_vec_real(npts);

    for ( integer i{0}; i < npts; ++i ) {
      real_type x;
      intersect( indep, zetas[i], x );
      for ( auto const & S : columns ) {
        vec_real_type & v{ vals[S].get_vec_real() };
        Spline const * p_spl{ get_spline( S ) };
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
  SplineSet::eval_DD( real_type const x, GenericContainer & gc ) const {
    map_type & vals{ gc.set_map() };
    for ( auto const & [fst, snd] : m_header_to_position )
      vals[fst] = m_splines[snd]->eval_DD(x);
  }

  //!
  //! Evaluate all the splines at `x` values contained in vec
  //! and fill a `map of vector` in a GenericContainer
  //!
  void
  SplineSet::eval_DD( vec_real_type const & vec, GenericContainer & gc ) const {
    integer const npts{ static_cast<integer>(vec.size()) };
    map_type & vals{ gc.set_map() };
    for ( auto const & [fst, snd] : m_header_to_position ) {
      vec_real_type & v = vals[fst].set_vec_real(npts);
      auto & p_spl{ m_splines[snd] };
      for ( integer i{0}; i < npts; ++i ) v[i] = p_spl->eval_DD(vec[i]);
    }
  }

  //!
  //! Evaluate all the splines at `x` and fill a map of values in
  //! a GenericContainer with keys in `columns`
  //!
  void
  SplineSet::eval_DD(
    real_type       const   x,
    vec_string_type const & columns,
    GenericContainer      & gc
  ) const {
    map_type & vals{ gc.set_map() };
    for ( auto const & S : columns ) {
      Spline const * p_spl{ get_spline( S ) };
      vals[S] = p_spl->eval_DD(x);
    }
  }

  //!
  //! Evaluate all the splines at `x` values contained in vec and
  //! fill a map of vector in a GenericContainer with keys in `columns`
  //!
  void
  SplineSet::eval_DD(
    vec_real_type   const & vec,
    vec_string_type const & columns,
    GenericContainer      & gc
  ) const {
    integer const npts{ static_cast<integer>(vec.size()) };
    map_type & vals{ gc.set_map() };
    for ( auto const & S : columns ) {
      Spline const * p_spl{ get_spline( S ) };
      vec_real_type & v{ vals[S].set_vec_real(npts) };
      for ( integer i{0}; i < npts; ++i ) v[i] = p_spl->eval_DD(vec[i]);
    }
  }

  //!
  //! Evaluate all the splines at `zeta` and fill a `map of values`
  //! in a GenericContainer and `indep` as independent spline
  //!
  void
  SplineSet::eval2_DD(
    real_type const    zeta,
    integer   const    indep,
    GenericContainer & gc
  ) const {
    map_type & vals = gc.set_map();
    real_type x;
    intersect( indep, zeta, x );
    for ( auto const & [fst, snd] : m_header_to_position )
      vals[fst] = m_splines[snd]->eval_DD(x);
  }

  //!
  //! Evaluate all the splines at `zeta` values contained in vec
  //! and fill a `map of vector` in a GenericContainer and `indep`
  //! as independent spline
  //!
  void
  SplineSet::eval2_DD(
    vec_real_type const & zetas,
    integer       const   indep,
    GenericContainer    & gc
  ) const {
    integer const npts{ static_cast<integer>(zetas.size()) };
    map_type & vals{ gc.set_map() };

    // pre-allocation
    for ( auto const & [fst, snd] : m_header_to_position )
      vals[fst].set_vec_real(npts);

    for ( integer i{0}; i < npts; ++i ) {
      real_type x;
      intersect( indep, zetas[i], x );
      for ( auto const & [fst, snd] : m_header_to_position ) {
        vec_real_type & v = vals[fst].get_vec_real();
        v[i] = m_splines[snd]->eval_DD(x);
      }
    }
  }

  //!
  //! Evaluate all the splines at `zeta` and fill a map of values
  //! in a GenericContainer with keys in `columns` and `indep`
  //! as independent spline
  //!
  void
  SplineSet::eval2_DD(
    real_type       const   zeta,
    integer         const   indep,
    vec_string_type const & columns,
    GenericContainer      & gc
  ) const {
    map_type & vals{ gc.set_map() };
    real_type x;
    intersect( indep, zeta, x );
    for ( auto const & S : columns ) {
      Spline const * p_spl{ get_spline( S ) };
      vals[S] = p_spl->eval_DD(x);
    }
  }

  //!
  //! Evaluate all the splines at `zeta` values contained in vec
  //! and fill a map of vector in a GenericContainer with keys
  //! in `columns` and `indep` as independent spline
  //!
  void
  SplineSet::eval2_DD(
    vec_real_type   const & zetas,
    integer         const   indep,
    vec_string_type const & columns,
    GenericContainer      & gc
  ) const {
    integer const npts{ static_cast<integer>(zetas.size()) };
    map_type & vals{ gc.set_map() };

    // pre-allocation
    for ( auto const & S : columns ) vals[S].set_vec_real(npts);

    for ( integer i{0}; i < npts; ++i ) {
      real_type x;
      intersect( indep, zetas[i], x );
      for ( auto const & S : columns ) {
        vec_real_type & v{ vals[S].get_vec_real() };
        Spline const * p_spl{ get_spline( S ) };
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
  SplineSet::eval_DDD( real_type const x, GenericContainer & gc ) const {
    map_type & vals{ gc.set_map() };
    for ( auto const & [fst, snd] : m_header_to_position )
      vals[fst] = m_splines[snd]->eval_DDD(x);
  }

  //!
  //! Evaluate all the splines at `x` values contained in vec
  //! and fill a `map of vector` in a GenericContainer
  //!
  void
  SplineSet::eval_DDD(
    vec_real_type const & vec,
    GenericContainer    & gc
  ) const {
    integer const npts{ static_cast<integer>(vec.size()) };
    map_type & vals{ gc.set_map() };
    for ( auto const & [fst, snd] : m_header_to_position ) {
      vec_real_type & v{ vals[fst].set_vec_real(npts) };
      Spline const * p_spl{ m_splines[snd].get() };
      for ( integer i{0}; i < npts; ++i ) v[i] = p_spl->eval_DDD(vec[i]);
    }
  }

  //!
  //! Evaluate all the splines at `x` and fill a map of values
  //! in a GenericContainer with keys in `columns`
  //!
  void
  SplineSet::eval_DDD(
    real_type       const   x,
    vec_string_type const & columns,
    GenericContainer      & gc
  ) const {
    map_type & vals{ gc.set_map() };
    for ( auto const & S : columns ) {
      Spline const * p_spl{ get_spline( S ) };
      vals[S] = p_spl->eval_DDD(x);
    }
  }

  //!
  //! Evaluate all the splines at `x` values contained in vec and
  //! fill a map of vector in a GenericContainer with keys in `columns`
  //!
  void
  SplineSet::eval_DDD(
    vec_real_type   const & vec,
    vec_string_type const & columns,
    GenericContainer      & gc
  ) const {
    integer const npts{ static_cast<integer>(vec.size()) };
    map_type & vals{ gc.set_map() };
    for ( auto const & S : columns ) {
      Spline const * p_spl{ get_spline( S ) };
      vec_real_type & v{ vals[S].set_vec_real(npts) };
      for ( integer i{0}; i < npts; ++i ) v[i] = p_spl->eval_DDD(vec[i]);
    }
  }

  //!
  //! Evaluate all the splines at `zeta` and fill a `map of values`
  //! in a GenericContainer and `indep` as independent spline
  //!
  void
  SplineSet::eval2_DDD(
    real_type const    zeta,
    integer   const    indep,
    GenericContainer & gc
  ) const {
    map_type & vals{ gc.set_map() };
    real_type x;
    intersect( indep, zeta, x );
    for ( auto const & [fst, snd] : m_header_to_position )
      vals[fst] = m_splines[snd]->eval_DDD(x);
  }

  //!
  //! Evaluate all the splines at `zeta` values contained in vec
  //! and fill a `map of vector` in a GenericContainer and `indep`
  //! as independent spline
  //!
  void
  SplineSet::eval2_DDD(
    vec_real_type const & zetas,
    integer       const   indep,
    GenericContainer    & gc
  ) const {
    integer const npts{ static_cast<integer>(zetas.size()) };
    map_type & vals{ gc.set_map() };

    // pre-allocation
    for ( auto const & [fst, snd] : m_header_to_position )
      vals[fst].set_vec_real(npts);

    for ( integer i{0}; i < npts; ++i ) {
      real_type x;
      intersect( indep, zetas[i], x );
      for ( auto const & [fst, snd] : m_header_to_position ){
        vec_real_type & v{ vals[fst].get_vec_real() };
        v[i] = m_splines[snd]->eval_DDD(x);
      }
    }
  }

  //!
  //! Evaluate all the splines at `zeta` and fill a map of values
  //! in a GenericContainer with keys in `columns` and `indep`
  //! as independent spline
  //!
  void
  SplineSet::eval2_DDD(
    real_type const         zeta,
    integer   const         indep,
    vec_string_type const & columns,
    GenericContainer      & gc
  ) const {
    map_type & vals{ gc.set_map() };
    real_type x;
    intersect( indep, zeta, x );
    for ( auto const & S : columns ) {
      Spline const * p_spl{ get_spline( S ) };
      vals[S] = p_spl->eval_DDD(x);
    }
  }

  //!
  //! Evaluate all the splines at `zeta` values contained in vec
  //! and fill a map of vector in a GenericContainer with keys
  //! in `columns` and `indep` as independent spline
  //!
  void
  SplineSet::eval2_DDD(
    vec_real_type   const & zetas,
    integer         const   indep,
    vec_string_type const & columns,
    GenericContainer      & gc
  ) const {
    integer const npts{ static_cast<integer>(zetas.size()) };
    map_type & vals{ gc.set_map() };

    // pre-allocation
    for ( auto const & S : columns )
      vals[S].set_vec_real( npts );

    for ( integer i{0}; i < npts; ++i ) {
      real_type x;
      intersect( indep, zetas[i], x );
      for ( auto const & S : columns ) {
        vec_real_type & v{ vals[S].get_vec_real() };
        Spline const * p_spl{ get_spline( S ) };
        v[i] = p_spl->eval_DDD(x);
      }
    }
  }
}
