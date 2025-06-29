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
#include "PolynomialRoots.hh"
#include "Utils_fmt.hh"

#include <limits>
#include <cmath>
#include <set>

namespace Splines {

  using std::abs;
  using std::sqrt;
  using std::copy_n;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer
  SplineSet::BinarySearch::search( string_view const id ) const {
    size_t U{ data.size() };
    size_t L{ 0 };
    while ( U-L > 1 ) {
      size_t const pos{ (L+U)>>1 }; // se L=U+1 --> (L+U)/2 ==> L
      string_view id_pos{ data[pos].first };
      if ( id_pos < id ) L = pos; else U = pos;
    }
    if ( data[L].first == id ) return data[L].second;
    if ( U < data.size() )
      if ( data[U].first == id )
        return data[U].second;
    return -1; // non trovato
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineSet::BinarySearch::insert( string_view const id, integer const position ) {
    size_t pos{ data.size() };
    data.emplace_back(id,position);
    while ( pos > 0 ) {
      size_t const pos1{pos-1};
      data[pos].first  = data[pos1].first;
      data[pos].second = data[pos1].second;
      if ( data[pos1].first < id ) break;
      pos = pos1;
    }
    data[pos] = DATA_TYPE(id,position);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  //! spline constructor
  SplineSet::SplineSet( string_view name )
  : m_name(name)
  , m_mem( fmt::format( "SplineSet[{}]::m_mem", name ) )
  , m_mem_p( fmt::format( "SplineSet[{}]::m_mem_p", name ) )
  , m_mem_int( fmt::format( "SplineSet[{}]::m_mem_int", name ) )
  {
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  //! spline destructor
  SplineSet::~SplineSet() {
    m_mem.free();
    m_mem_p.free();
    m_mem_int.free();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  string
  SplineSet::info() const {
    string res{ fmt::format( "SplineSet[{}] n.points={} n.splines={}", name(), m_npts, m_nspl ) };

    for ( integer i{0}; i < m_nspl; ++i ) {
      res += fmt::format("\nSpline n.{} ", i);
      switch ( m_is_monotone[i] ) {
        case -2: res += " with NON monotone data\n"; break;
        case -1: res += " is NOT monotone\n";        break;
        case  0: res += " is monotone\n";            break;
        case  1: res += " is strictly monotone\n";   break;
        default: UTILS_ERROR(
          "SplineSet::info classification: {} not in range {-2,-1,0,1}\n",
          m_is_monotone[i]
        );
      }
      res += m_splines[i]->info();
    }
    return res;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineSet::get_headers( std::vector<std::string> & names ) const {
    names.clear();
    names.reserve( m_nspl );
    for ( integer i{0}; i < m_nspl; ++i )
      names.emplace_back( m_splines[i]->name() );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  string
  SplineSet::name_list() const {
    string tmp{"[ "};
    for ( integer i{0}; i < m_nspl; ++i )
      tmp += fmt::format( "'{}' ", m_splines[i]->name() );
    tmp += "]";
    return tmp;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type const *
  SplineSet::y_nodes( integer i ) const {
    UTILS_ASSERT(
      i >= 0 && i < m_nspl,
      "SplineSet[{}]::y_nodes({}) argument out of range [0,{}]\n",
      m_name, i, m_nspl-1
    );
    return m_Y[i];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Spline *
  SplineSet::get_spline( integer i ) const {
    UTILS_ASSERT(
      i >= 0 && i < m_nspl,
      "SplineSet[{}]::get_spline({}) argument out of range [0,{}]\n",
      m_name, i, m_nspl-1
    );
    return m_splines[i].get();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineSet::dump_table( ostream_type & stream, integer const num_points ) const {
    Malloc_real mem("SplineSet::dump_table");
    mem.allocate( m_nspl );
    real_type * vals{ mem( m_nspl ) };
    stream << 's';
    for ( integer i{0}; i < m_nspl; ++i ) stream << '\t' << header(i);
    stream << '\n';

    for ( integer j{0}; j < num_points; ++j ) {
      real_type const s{ x_min() + ((x_max()-x_min())*j)/(num_points-1) };
      this->eval( s, vals, 1 );
      stream << s;
      for ( integer i{0}; i < m_nspl; ++i ) stream << '\t' << vals[i];
      stream << '\n';
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer
  SplineSet::get_position( string_view hdr ) const {
    integer const pos{ m_header_to_position.at(hdr.data()) };
    UTILS_ASSERT(
      pos >= 0 && pos < m_nspl,
      "SplineSet[{}]::get_position(\"{}\") not found!\n"
      "available keys: {}\n",
      m_name, hdr, name_list()
    );
    return pos;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineSet::deep_copy_to( SplineSet & S ) const {

    std::vector<char const *> headers(m_nspl);
    std::vector<SplineType1D> stype(m_nspl);

    for ( integer i{0}; i < m_nspl; ++i ) {
      headers[i] = m_splines[i]->name().data();
      stype[i]   = m_splines[i]->type();
    }

    S.build( m_nspl, m_npts, headers.data(), stype.data(), m_X, m_Y, m_Yp );

    // propagate flags
    for ( integer i{0}; i < m_nspl; ++i ) S.m_splines[i]->copy_flags(*m_splines[i]);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineSet::build(
    integer      const         nspl,
    integer      const         npts,
    char         const * const headers[],
    SplineType1D const         stype[],
    real_type    const         data_X[],
    real_type    const * const data_Y[],
    real_type    const * const data_Yp[]
  ) {
    string msg{ fmt::format("SplineSet[{}]::build(...):", m_name ) };
    UTILS_ASSERT( nspl > 0, "{} expected positive nspl = {}\n", msg, nspl );
    UTILS_ASSERT( npts > 1, "{} expected npts = {} greater than 1", msg, npts );
    m_nspl = nspl;
    m_npts = npts;
    // allocate memory
    m_splines.resize( m_nspl );
    m_is_monotone = m_mem_int.realloc( m_nspl );

    m_header_to_position.clear();

    integer mem{npts};
    for ( integer spl{0}; spl < nspl; ++spl ) {
      switch (stype[spl]) {
      case SplineType1D::QUINTIC:
        mem += npts; // Y, Yp, Ypp
        [[fallthrough]];
      case SplineType1D::CUBIC:
      case SplineType1D::AKIMA:
      case SplineType1D::BESSEL:
      case SplineType1D::PCHIP:
      case SplineType1D::HERMITE:
        mem += npts; // Y, Yp
        [[fallthrough]];
      case SplineType1D::CONSTANT:
      case SplineType1D::LINEAR:
        mem += npts;
        break;
      case SplineType1D::SPLINE_SET:
      case SplineType1D::SPLINE_VEC:
      //default:
        UTILS_ERROR(
          "{} At spline n.{} named {} cannot be done for type = {}\n",
          msg, spl, headers[spl], to_string(stype[spl])
        );
      }
    }

    m_mem.reallocate( mem + 2*nspl );
    m_mem_p.reallocate( 3*nspl );

    m_Y    = m_mem_p ( m_nspl );
    m_Yp   = m_mem_p ( m_nspl );
    m_Ypp  = m_mem_p ( m_nspl );
    m_X    = m_mem   ( m_npts );
    m_Ymin = m_mem   ( m_nspl );
    m_Ymax = m_mem   ( m_nspl );

    copy_n( data_X, npts, m_X );
    for ( integer spl{0}; spl < nspl; ++spl ) {
      real_type * & pY{ m_Y[spl] };
      real_type * & pYp{ m_Yp[spl] };
      real_type * & pYpp{ m_Ypp[spl] };
      pY = m_mem( m_npts );
      copy_n( data_Y[spl], npts, pY );
      if ( stype[spl] == SplineType1D::CONSTANT ) {
        m_Ymin[spl] = *std::min_element( pY, pY+npts-1 );
        m_Ymax[spl] = *std::max_element( pY, pY+npts-1 );
      } else {
        m_Ymin[spl] = *std::min_element( pY, pY+npts );
        m_Ymax[spl] = *std::max_element( pY, pY+npts );
      }
      pYpp = pYp = nullptr;
      switch ( stype[spl] ) {
      case SplineType1D::QUINTIC:
        pYpp = m_mem( m_npts );
        [[fallthrough]];
      case SplineType1D::CUBIC:
      case SplineType1D::AKIMA:
      case SplineType1D::BESSEL:
      case SplineType1D::PCHIP:
      case SplineType1D::HERMITE:
        pYp = m_mem( m_npts );
        if ( stype[spl] == SplineType1D::HERMITE ) {
          UTILS_ASSERT(
            data_Yp != nullptr && data_Yp[spl] != nullptr,
            "{} At spline n.{} named {}\n"
            "expect to find derivative values",
            msg, spl, headers[spl]
          );
          copy_n( data_Yp[spl], npts, pYp );
        }
        [[fallthrough]];
      case SplineType1D::CONSTANT:
      case SplineType1D::LINEAR:
      case SplineType1D::SPLINE_SET:
      case SplineType1D::SPLINE_VEC:
      //default:
        break;
      }
      string h{ headers[spl] };
      std::unique_ptr<Spline> & s{ m_splines[spl] };

      m_is_monotone[spl] = -1;
      switch ( stype[ spl ] ) {
      case SplineType1D::CONSTANT:
        { auto S = std::make_unique<ConstantSpline>(h);
          S->reserve_external( m_npts, m_X, pY );
          S->m_npts = m_npts;
          S->build();
          s = std::move(S);
        }
        break;

      case SplineType1D::LINEAR:
        { auto S = std::make_unique<LinearSpline>(h);
          S->reserve_external( m_npts, m_X, pY );
          S->m_npts = m_npts;
          S->build();
          // check monotonicity of data
          integer flag{1};
          for ( integer j{1}; j < m_npts; ++j ) {
            if ( pY[j-1] > pY[j] ) { flag = -1; break; } // non monotone data
            if ( Utils::is_zero(pY[j-1]-pY[j]) && m_X[j-1] < m_X[j] ) flag = 0; // non strict monotone
          }
          m_is_monotone[spl] = flag;
          s = std::move(S);
        }
        break;

      case SplineType1D::CUBIC:
        { auto S = std::make_unique<CubicSpline>(h);
          S->reserve_external( m_npts, m_X, pY, pYp );
          S->m_npts = m_npts;
          S->build();
          m_is_monotone[spl] = check_cubic_spline_monotonicity( m_X, pY, pYp, m_npts );
          s = std::move(S);
        }
        break;

      case SplineType1D::AKIMA:
        { auto S = std::make_unique<AkimaSpline>(h);
          S->reserve_external( m_npts, m_X, pY, pYp );
          S->m_npts = m_npts;
          S->build();
          m_is_monotone[spl] = check_cubic_spline_monotonicity( m_X, pY, pYp, m_npts );
          s = std::move(S);
        }
        break;

      case SplineType1D::BESSEL:
        { auto S = std::make_unique<BesselSpline>(h);
          S->reserve_external( m_npts, m_X, pY, pYp );
          S->m_npts = m_npts;
          S->build();
          m_is_monotone[spl] = check_cubic_spline_monotonicity( m_X, pY, pYp, m_npts );
          s = std::move(S);
        }
        break;

      case SplineType1D::PCHIP:
        { auto S = std::make_unique<PchipSpline>(h);
          S->reserve_external( m_npts, m_X, pY, pYp );
          S->m_npts = m_npts;
          S->build();
          m_is_monotone[spl] = check_cubic_spline_monotonicity( m_X, pY, pYp, m_npts );
          s = std::move(S);
        }
        break;

      case SplineType1D::HERMITE:
        { auto S = std::make_unique<HermiteSpline>(h);
          S->reserve_external( m_npts, m_X, pY, pYp );
          S->m_npts = m_npts;
          S->build();
          m_is_monotone[spl] = check_cubic_spline_monotonicity( m_X, pY, pYp, m_npts );
          s = std::move(S);
        }
        break;

      case SplineType1D::QUINTIC:
        { auto S = std::make_unique<QuinticSpline>(h);
          S->reserve_external( m_npts, m_X, pY, pYp, pYpp );
          S->m_npts = m_npts;
          S->build();
          s = std::move(S);
        }
        break;

      case SplineType1D::SPLINE_SET:
      case SplineType1D::SPLINE_VEC:
        //default:
        UTILS_ERROR(
          "{} At spline n.{} named {}\n"
          "{} not allowed as spline type\n"
          "in SplineSet::build for {}-th spline\n",
          msg, spl, headers[spl], to_string(stype[ spl ]), spl
        );
      }
      m_header_to_position.insert( {s->name().data(), static_cast<integer>(spl)} );
    }

    m_mem.must_be_empty( "SplineSet::build, baseValue" );
    m_mem_p.must_be_empty( "SplineSet::build, basePointer" );
  }
  
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineSet::eval( real_type const x, vector<real_type> & vals ) const {
    vals.resize( m_nspl );
    for ( integer i{0}; i < m_nspl; ++i )
      vals[i] = m_splines[i]->eval(x);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineSet::eval(
    real_type const x,
    real_type       vals[],
    integer   const incy
  ) const {
    integer ii{0};
    for ( integer i{0}; i < m_nspl; ++i, ii += incy )
      vals[ii] = m_splines[i]->eval(x);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineSet::eval_D( real_type const x, vector<real_type> & vals ) const {
    vals.resize( m_nspl );
    for ( integer i{0}; i < m_nspl; ++i )
      vals[i] = m_splines[i]->D(x);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineSet::eval_D(
    real_type const x,
    real_type       vals[],
    integer   const incy
  ) const {
    size_t ii{0};
    for ( integer i{0}; i < m_nspl; ++i, ii += incy )
      vals[ii] = m_splines[i]->D(x);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineSet::eval_DD( real_type const x, vector<real_type> & vals ) const {
    vals.resize( m_nspl );
    for ( integer i{0}; i < m_nspl; ++i )
      vals[i] = m_splines[i]->DD(x);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineSet::eval_DD(
    real_type const x,
    real_type       vals[],
    integer   const incy
  ) const {
    size_t ii{0};
    for ( integer i{0}; i < m_nspl; ++i, ii += incy )
      vals[ii] = m_splines[i]->DD(x);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineSet::eval_DDD( real_type const x, vector<real_type> & vals ) const {
    vals.resize( m_nspl );
    for ( integer i{0}; i < m_nspl; ++i )
      vals[i] = m_splines[i]->DDD(x);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineSet::eval_DDD(
    real_type const x,
    real_type       vals[],
    integer   const incy
  ) const {
    size_t ii{0};
    for ( integer i{0}; i < m_nspl; ++i, ii += incy )
      vals[ii] = m_splines[i]->DDD(x);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // vectorial values

  Spline const *
  SplineSet::intersect(
    integer   const spl,
    real_type const zeta,
    real_type &     x
  ) const {
    string msg = fmt::format("SplineSet[{}]::intersect(...):", m_name );
    UTILS_ASSERT(
      spl >= 0 && spl < m_nspl,
      "{}\nSpline n.{} is not in SplineSet", msg, spl
    );
    UTILS_ASSERT(
      m_is_monotone[ spl ] > 0,
      "{}\nSpline n.{} is not monotone and can't be used as independent",
      msg, spl
    );
    auto & S{ m_splines[spl] };
    // cerco intervallo intersezione
    real_type const * X{ m_Y[spl] };
    UTILS_ASSERT(
      zeta >= X[0] && zeta <= X[m_npts-1],
      "{} evaluation at zeta = {} is out of range: [{},{}]\n",
      msg, zeta, X[0], X[m_npts-1]
    );

    integer interval{ static_cast<integer>(lower_bound(X, X + m_npts, zeta) - X) };
    if ( interval > 0 ) --interval;
    if ( Utils::is_zero(X[interval]-X[interval+1]) ) ++interval; // degenerate interval for duplicated nodes
    if ( interval >= m_npts-1 ) interval = m_npts-2;

    // compute intersection
    real_type const a{ m_X[interval] };
    real_type const b{ m_X[interval+1] };
    real_type const ya{ X[interval] };
    real_type const yb{ X[interval+1] };
    real_type const DX{ b-a };
    real_type const DY{ yb-ya };
    UTILS_ASSERT(
      zeta >= ya && zeta <= yb,
      "{} Bad interval [{},{}] for zeta = {}\n", msg, ya, yb, zeta
    );
    UTILS_ASSERT(
      a < b,
      "{} Bad x interval [{},{}]\n", msg,  a, b
    );
    if ( S->type() == SplineType1D::LINEAR ) {
      x = a + (b-a)*(zeta-ya)/(yb-ya);
    } else {
      real_type const * dX = m_Yp[spl];
      real_type const  dya = dX[interval];
      real_type const  dyb = dX[interval+1];
      PolynomialRoots::Cubic const cubic(
        (dyb+dya-2*DY/DX)/(DX*DX), // A
        (3*DY/DX-2*dya-dyb)/DX,    // B
        dya,                       // C
        ya-zeta                    // D
      );
      real_type r[3];
      integer const npr{ cubic.getRealRoots( r ) };
      // cerca radice buona
      bool ok{ false };
      for ( integer i{0}; i < npr && !ok; ++i ) {
        ok = r[i] >= 0 && r[i] <= DX;
        if ( ok ) x = a + r[i];
      }
      UTILS_ASSERT(
        ok,
        "{}\nfailed to find intersection with independent spline at zeta = {}\n",
        msg, zeta
      );
    }
    return S.get();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineSet::eval2(
    integer   const indep,
    real_type const zeta,
    real_type       vals[],
    integer   const incy
  ) const {
    real_type x;
    intersect( indep, zeta, x );
    size_t ii{0};
    for ( integer i{0}; i < m_nspl; ++i, ii += incy )
      vals[ii] = m_splines[i]->eval(x);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineSet::eval2(
    integer   const     indep,
    real_type const     zeta,
    vector<real_type> & vals
  ) const {
    vals.resize( m_nspl );
    this->eval2( indep, zeta, vals.data(), 1 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  SplineSet::eval2( real_type const zeta, integer const indep, integer const spl ) const {
    real_type x;
    intersect( indep, zeta, x );
    return m_splines[spl]->eval(x);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  SplineSet::eval2(
    real_type   const zeta,
    string_view const indep,
    string_view const name
  ) const {
    return this->eval2(
      zeta, this->get_position(indep), this->get_position(name)
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineSet::eval2_D(
    integer   const indep,
    real_type const zeta,
    real_type       vals[],
    integer   const incy
  ) const {
    real_type x;
    Spline    const * S{ intersect( indep, zeta, x ) };
    real_type const   ds{S->D(x)};
    size_t ii{0};
    for ( integer i{0}; i < m_nspl; ++i, ii += incy )
      vals[ii] = m_splines[i]->D(x)/ds;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineSet::eval2_D(
    integer   const     spl,
    real_type const     zeta,
    vector<real_type> & vals
  ) const {
    vals.resize( m_nspl );
    this->eval2_D( spl, zeta, vals.data(), 1 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  SplineSet::eval2_D( real_type const zeta, integer const indep, integer const spl ) const {
    real_type x;
    Spline const * S{ intersect( indep, zeta, x ) };
    return m_splines[spl]->D(x)/S->D(x);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  SplineSet::eval2_D(
    real_type   const zeta,
    string_view const indep,
    string_view const name
  ) const {
    return this->eval2_D(
      zeta, this->get_position(indep), this->get_position(name)
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineSet::eval2_DD(
    integer   const indep,
    real_type const zeta,
    real_type       vals[],
    integer   const incy
  ) const {
    real_type x;
    Spline const * S{ intersect( indep, zeta, x ) };
    real_type const dt{ 1/S->D(x) };
    real_type const dt2{ dt*dt };
    real_type const ddt{ -S->DD(x)*(dt*dt2) };
    size_t ii{0};
    for ( integer i{0}; i < m_nspl; ++i, ii += incy ) {
      auto & Si{ m_splines[i] };
      vals[ii] = Si->DD(x)*dt2 + Si->D(x)*ddt;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineSet::eval2_DD(
    integer   const     indep,
    real_type const     zeta,
    vector<real_type> & vals
  ) const {
    vals.resize( m_nspl );
    this->eval2_DD( indep, zeta, vals.data(), 1 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  SplineSet::eval2_DD( real_type const zeta, integer const indep, integer const spl ) const {
    real_type x;
    Spline const * S{ intersect( indep, zeta, x ) };
    real_type const dt{ 1/S->D(x) };
    real_type const dt2{ dt*dt };
    real_type const ddt{ -S->DD(x)*(dt*dt2) };
    auto & SPL{ m_splines[spl] };
    return SPL->DD(x)*dt2 + SPL->D(x)*ddt;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  SplineSet::eval2_DD(
    real_type   const zeta,
    string_view const indep,
    string_view const name
  ) const {
    return this->eval2_DD( zeta, this->get_position(indep), this->get_position(name) );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineSet::eval2_DDD(
    integer   const indep,
    real_type const zeta,
    real_type       vals[],
    integer   const incy
  ) const {
    real_type x;
    Spline const * S{ intersect( indep, zeta, x ) };
    real_type const dt{ 1/S->D(x) };
    real_type const dt3{ dt*dt*dt };
    real_type const ddt{ -S->DD(x)*dt3 };
    real_type const dddt{ 3*(ddt*ddt)/dt-S->DDD(x)*(dt*dt3) };
    size_t ii{0};
    for ( integer i{0}; i < m_nspl; ++i, ii += incy ) {
      auto & Si{ m_splines[i] };
      vals[ii] = Si->DDD(x)*dt3 + 3*Si->DD(x)*dt*ddt + Si->D(x)*dddt;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineSet::eval2_DDD(
    integer   const     spl,
    real_type const     zeta,
    vector<real_type> & vals
  ) const {
    vals.resize( m_nspl );
    this->eval2_DDD( spl, zeta, vals.data(), 1 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  SplineSet::eval2_DDD( real_type const zeta, integer const indep, integer const spl ) const {
    real_type x;
    Spline const * S{ intersect( indep, zeta, x ) };
    real_type const dt{ 1/S->D(x) };
    real_type const dt3{ dt*dt*dt };
    real_type const ddt{ -S->DD(x)*dt3 };
    real_type const dddt{ 3*(ddt*ddt)/dt-S->DDD(x)*(dt*dt3) };

    auto & SPL{ m_splines[spl] };
    return SPL->DDD(x)*dt3 + 3*SPL->DD(x)*dt*ddt + SPL->D(x)*dddt;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  SplineSet::eval2_DDD(
    real_type   const zeta,
    string_view const indep,
    string_view const name
  ) const {
    return this->eval2_DDD(
      zeta, this->get_position(indep), this->get_position(name)
    );
  }
}
