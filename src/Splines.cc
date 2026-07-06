/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2026                                                      |
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

//!
//! Namespace of Splines library
//!
namespace Splines
{

  SplineType1D string_to_splineType1D( string_view nin )
  {
    string n{ nin };
    std::transform( n.begin(), n.end(), n.begin(), []( unsigned char c ) { return char( std::tolower( c ) ); } );
    if ( n == "constant" ) return SplineType1D::CONSTANT;
    if ( n == "linear" ) return SplineType1D::LINEAR;

    if ( n == "cubic" ) return SplineType1D::CUBIC;
    if ( n == "akima" ) return SplineType1D::AKIMA;
    if ( n == "vanleer" ) return SplineType1D::VANLEER;
    if ( n == "pchip" ) return SplineType1D::PCHIP;

    if ( n == "quintic_cubic" ) return SplineType1D::QUINTIC_CUBIC;
    if ( n == "quintic_akima" ) return SplineType1D::QUINTIC_AKIMA;
    if ( n == "quintic_vanleer" ) return SplineType1D::QUINTIC_VANLEER;
    if ( n == "quintic_pchip" || n == "quintic" ) return SplineType1D::QUINTIC_PCHIP;

    if ( n == "hermite" ) return SplineType1D::HERMITE;
    if ( n == "spline_set" ) return SplineType1D::SPLINE_SET;
    if ( n == "spline_vec" ) return SplineType1D::SPLINE_VEC;
    throw std::runtime_error( fmt::format( "string_to_splineType1D({}) unknown type\n", n ) );
  }

  string to_string( Spline_sub_type t )
  {
    switch ( t )
    {
      case Spline_sub_type::CUBIC: return "CUBIC";
      case Spline_sub_type::PCHIP: return "PCHIP";
      case Spline_sub_type::AKIMA: return "AKIMA";
      case Spline_sub_type::VANLEER: return "VANLEER";
    }
    return "NOTYPE";
  }

  SplineType2D string_to_splineType2D( string_view nin )
  {
    string n{ nin };
    std::transform( n.begin(), n.end(), n.begin(), []( unsigned char c ) { return char( std::tolower( c ) ); } );
    if ( n == "bilinear" ) return SplineType2D::BILINEAR;

    if ( n == "bicubic" || n == "bicubic_cubic" ) return SplineType2D::BICUBIC_CUBIC;
    if ( n == "bicubic_akima" ) return SplineType2D::BICUBIC_AKIMA;
    if ( n == "bicubic_vanleer" ) return SplineType2D::BICUBIC_VANLEER;
    if ( n == "bicubic_pchip" ) return SplineType2D::BICUBIC_PCHIP;

    if ( n == "biquintic" || n == "biquintic_cubic" ) return SplineType2D::BIQUINTIC_CUBIC;
    if ( n == "biquintic_akima" ) return SplineType2D::BIQUINTIC_AKIMA;
    if ( n == "biquintic_vanleer" ) return SplineType2D::BIQUINTIC_VANLEER;
    if ( n == "biquintic_pchip" ) return SplineType2D::BIQUINTIC_PCHIP;

    throw std::runtime_error( fmt::format( "string_to_splineType2D({}) unknown type\n", n ) );
  }

  char const * to_string( SplineType1D const t )
  {
    switch ( t )
    {
      case SplineType1D::CONSTANT: return "SPLINE_CONSTANT";
      case SplineType1D::LINEAR: return "SPLINE_LINEAR";
      case SplineType1D::CUBIC: return "SPLINE_CUBIC";
      case SplineType1D::AKIMA: return "SPLINE_AKIMA";
      case SplineType1D::VANLEER: return "SPLINE_VANLEER";
      case SplineType1D::PCHIP: return "SPLINE_PCHIP";
      case SplineType1D::QUINTIC_CUBIC: return "SPLINE_QUINTIC";
      case SplineType1D::QUINTIC_AKIMA: return "SPLINE_QUINTIC_AKIMA";
      case SplineType1D::QUINTIC_VANLEER: return "SPLINE_QUINTIC_VANLEER";
      case SplineType1D::QUINTIC_PCHIP: return "SPLINE_QUINTIC_PCHIP";
      case SplineType1D::HERMITE: return "SPLINE_HERMITE";
      case SplineType1D::SPLINE_SET: return "SPLINE_SPLINE_SET";
      case SplineType1D::SPLINE_VEC: return "SPLINE_SPLINE_VEC";
    }
    return "NO_TYPE";
  }

  char const * to_string( SplineType2D const t )
  {
    switch ( t )
    {
      case SplineType2D::BILINEAR: return "SPLINE2D_BILINEAR";
      case SplineType2D::BICUBIC_CUBIC: return "SPLINE2D_BICUBIC_CUBIC";
      case SplineType2D::BICUBIC_AKIMA: return "SPLINE2D_BICUBIC_AKIMA";
      case SplineType2D::BICUBIC_VANLEER: return "SPLINE2D_BICUBIC_VANLEER";
      case SplineType2D::BICUBIC_PCHIP: return "SPLINE2D_BICUBIC_PCHIP";
      case SplineType2D::BIQUINTIC_CUBIC: return "SPLINE2D_BIQUINTIC_CUBIC";
      case SplineType2D::BIQUINTIC_AKIMA: return "SPLINE2D_BIQUINTIC_AKIMA";
      case SplineType2D::BIQUINTIC_VANLEER: return "SPLINE2D_BIQUINTIC_VANLEER";
      case SplineType2D::BIQUINTIC_PCHIP: return "SPLINE2D_BIQUINTIC_PCHIP";
    }
    return "NO_TYPE";
  }

  /*\
   |   ____        _ _
   |  / ___| _ __ | (_)_ __   ___
   |  \___ \| '_ \| | | '_ \ / _ \
   |   ___) | |_) | | | | | |  __/
   |  |____/| .__/|_|_|_| |_|\___|
   |        |_|
  \*/

  void Spline::setup( GenericContainer const & gc )
  {
    /*
    // gc["xdata"]
    // gc["ydata"]
    //
    */
    string const where = fmt::format( "Spline[{}]::setup( gc ):", m_name );

    std::set<std::string> keywords;
    for ( auto const & pair : gc.get_map( where ) ) { keywords.insert( pair.first ); }

    GenericContainer const & gc_x = gc( "xdata", where );
    keywords.erase( "xdata" );
    GenericContainer const & gc_y = gc( "ydata", where );
    keywords.erase( "ydata" );
    keywords.erase( "spline_type" );

    vec_real_type x, y;
    {
      string const ff = fmt::format( "{}, field `xdata'", where );
      gc_x.copyto_vec_real( x, ff );
    }
    {
      string const ff = fmt::format( "{}, field `ydata'", where );
      gc_y.copyto_vec_real( y, ff );
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
    build( x, y );
  }

  void Spline::push_back( real_type const x, real_type const y )
  {
    if ( m_npts > 0 )
    {
      UTILS_ASSERT(
        x >= m_X[m_npts - 1],  // ammetto punti doppi
        "Spline[{}]::push_back, non monotone insert at insert N.{}"
        "\nX[{}] = {:.5}\nX[{}] = {:>5}\n",
        m_name,
        m_npts,
        m_npts - 1,
        m_X[m_npts - 1],
        m_npts,
        x );
    }
    if ( m_npts >= m_npts_reserved )
    {
      // Versione "Fall-back" se devi per forza usare il buffer temporaneo
      integer const saved_npts = m_npts;
      Vec           Xsaved( m_npts );
      Vec           Ysaved( m_npts );

      // memcpy è più veloce di copy_n per array raw
      std::memcpy( Xsaved.data(), m_X, m_npts * sizeof( real_type ) );
      std::memcpy( Ysaved.data(), m_Y, m_npts * sizeof( real_type ) );

      reserve( ( m_npts + 1 ) * 2 );  // La tua funzione distruttiva

      m_npts = saved_npts;
      std::memcpy( m_X, Xsaved.data(), m_npts * sizeof( real_type ) );
      std::memcpy( m_Y, Ysaved.data(), m_npts * sizeof( real_type ) );
    }
    m_X[m_npts] = x;
    m_Y[m_npts] = y;
    ++m_npts;
    m_search.must_reset();
  }

  void Spline::set_range( real_type xmin, real_type xmax )
  {
    UTILS_ASSERT( xmax > xmin, "Spline[{}]::set_range({},{}) bad range ", m_name, xmin, xmax );
    UTILS_ASSERT( m_npts > 1, "Spline[{}]::set_range: empty spline or only one point", m_name );

    // Calcolo range attuale (delta)
    real_type const dx_old = m_X[m_npts - 1] - m_X[0];

    // Opzionale: Protezione divisione per zero se la spline è collassata su un punto
    // UTILS_ASSERT( std::abs(dx_old) > epsilon, ... );

    real_type const S  = ( xmax - xmin ) / dx_old;
    real_type const Tx = xmin - S * m_X[0];

    // EIGEN IMPLEMENTATION
    // Usiamo Eigen::Array invece di Matrix perché l'operazione è element-wise (scalare)
    // Map crea una vista sui dati esistenti senza copiare memoria.
    Eigen::Map<Vec> map_X{ m_X, m_npts };

    // Operazione in-place: x[i] = x[i] * S + Tx
    // Eigen utilizzerà istruzioni vettoriali (es. vfmadd su AVX2)
    map_X = map_X * S + Tx;
  }

  void Spline::dump( ostream_type & s, integer const nintervals, string_view header ) const
  {
    s << header << '\n';
    real_type const dx = ( x_max() - x_min() ) / nintervals;
    for ( integer i = 0; i <= nintervals; ++i )
    {
      real_type x = x_min() + i * dx;
      fmt::print( s, "{}\t{}\n", x, this->eval( x ) );
    }
  }

  // Helper function: Determinante 2D (analogo al prodotto vettoriale "cross" in 3D)
  // Restituisce v1.x * v2.y - v1.y * v2.x
  static inline real_type kross( const Vec2 & a, const Vec2 & b )
  { return a.x() * b.y() - a.y() * b.x(); };

  real_type curvature( real_type const s, Spline const & X, Spline const & Y )
  {
    // Eigen caricherà questi valori direttamente nei registri CPU
    Vec2 v( X.D( s ), Y.D( s ) );    // Velocità
    Vec2 a( X.DD( s ), Y.DD( s ) );  // Accelerazione

    const real_type speed2 = v.squaredNorm();  // |v|^2

    // Protezione divisione per zero
    if ( speed2 <= 100 * std::numeric_limits<real_type>::epsilon() ) return 0;

    // Numeratore: v x a
    const real_type num = kross( v, a );

    // Denominatore: |v|^3 = speed2 * sqrt(speed2)
    return num * pow( speed2, -1.5 );
  }

  real_type curvature_D( real_type const s, Spline const & X, Spline const & Y )
  {
    Vec2 v( X.D( s ), Y.D( s ) );
    Vec2 a( X.DD( s ), Y.DD( s ) );
    Vec2 j( X.DDD( s ), Y.DDD( s ) );  // Jerk (Derivata terza)

    const real_type v2 = v.squaredNorm();  // |v|^2
    if ( v2 <= 1e-12 ) return 0;

    // Termini vettoriali comuni
    const real_type v_cross_a = kross( v, a );
    const real_type v_dot_a   = v.dot( a );
    const real_type v_cross_j = kross( v, j );

    // Applicazione formula derivata del quoziente vettoriale
    const real_type num = v2 * v_cross_j - 3 * v_cross_a * v_dot_a;

    // Denominatore: |v|^5
    return num / ( v2 * v2 * std::sqrt( v2 ) );
  }

  real_type curvature_DD( real_type const s, Spline const & X, Spline const & Y )
  {
    // Recupero dati (idealmente la classe Spline dovrebbe avere una funzione eval_all(s, ...) per farlo in una
    // chiamata)
    Vec2 v( X.D( s ), Y.D( s ) );
    Vec2 a( X.DD( s ), Y.DD( s ) );
    Vec2 j( X.DDD( s ), Y.DDD( s ) );
    Vec2 s_snap( X.DDDD( s ), Y.DDDD( s ) );  // Snap/Jounce (Derivata quarta)

    const real_type v2 = v.squaredNorm();
    if ( v2 <= 1e-12 ) return 0;

    // --- Precalcoli Termini vettoriali (Cache nei registri) ---
    const real_type v_dot_a   = v.dot( a );
    const real_type v_cross_a = kross( v, a );
    const real_type v_cross_j = kross( v, j );

    // Termini misti derivata seconda
    // T2 = (a x j) + (v x s)
    const real_type T2 = kross( a, j ) + kross( v, s_snap );

    // D2 = |a|^2 + v . j
    const real_type D2 = a.squaredNorm() + v.dot( j );

    // --- Calcolo Numeratore ---
    // Formula derivata analiticamente:
    // N = v^2 [ v^2 * T2  -  6 * (v.a) * (v x j)  -  3 * (v x a) * D2 ]  +  15 * (v x a) * (v.a)^2

    const real_type term_bracket = v2 * T2 - 6 * v_dot_a * v_cross_j - 3 * v_cross_a * D2;
    const real_type num          = v2 * term_bracket + 15 * v_cross_a * ( v_dot_a * v_dot_a );

    // --- Risultato Finale ---
    // Denominatore: |v|^7
    return num / ( v2 * v2 * v2 * std::sqrt( v2 ) );
  }

}  // namespace Splines
