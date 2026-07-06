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

#pragma once

#ifndef SPLINES_HH
#define SPLINES_HH

#ifdef __GNUC__
#pragma GCC diagnostic push
#endif

#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wpoison-system-directories"
#pragma clang diagnostic ignored "-Wundefined-func-template"
#pragma clang diagnostic ignored "-Wsign-conversion"
#endif

#include "SplinesConfig.hh"

#include "Utils_fmt.hh"
#include "Utils_FD.hh"
#include "Utils_search_intervals2.hh"
#include "Utils_TridiagonalSolver.hh"
#include "Utils_AlgoHNewton.hh"

//!
//! Namespace of Splines library
//!
namespace Splines
{

  using std::basic_istream;
  using std::basic_ostream;
  using std::cerr;
  using std::cin;
  using std::copy_n;
  using std::cout;
  using std::exception;
  using std::lower_bound;
  using std::pair;
  using std::runtime_error;
  using std::string;
  using std::string_view;
  using std::vector;

  using real_type = double;  //!< Floating point type for splines
  using integer   = int;     //!< Signed integer type for splines

  using Malloc_real  = Utils::Malloc<real_type>;
  using ostream_type = basic_ostream<char>;
  using istream_type = basic_istream<char>;

  using Vec  = Eigen::Array<real_type, Eigen::Dynamic, 1>;
  using Mat  = Eigen::Array<real_type, Eigen::Dynamic, Eigen::Dynamic>;
  using MatC = Eigen::Array<real_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

  using Vec2   = Eigen::Matrix<real_type, 2, 1>;
  using Vec4   = Eigen::Matrix<real_type, 4, 1>;
  using Vec6   = Eigen::Matrix<real_type, 6, 1>;
  using Mat2x2 = Eigen::Matrix<real_type, 2, 2>;
  using Mat4x4 = Eigen::Matrix<real_type, 4, 4>;
  using Mat6x6 = Eigen::Matrix<real_type, 6, 6>;

  using GC_namespace::GC_type;
  using GC_namespace::map_type;
  using GC_namespace::mat_real_type;
  using GC_namespace::vec_int_type;
  using GC_namespace::vec_real_type;
  using GC_namespace::vec_string_type;
  using GC_namespace::vector_type;

  void backtrace( ostream_type & );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  //! Associate a number for each type of splines implemented
  using SplineType1D = enum class SplineType1D : integer {
    CONSTANT        = 0,
    LINEAR          = 1,
    CUBIC           = 2,
    AKIMA           = 3,
    VANLEER         = 4,
    PCHIP           = 5,
    QUINTIC_CUBIC   = 6,
    QUINTIC_AKIMA   = 7,
    QUINTIC_VANLEER = 8,
    QUINTIC_PCHIP   = 9,
    HERMITE         = 10,
    SPLINE_SET      = 11,
    SPLINE_VEC      = 12
  };

  using Spline_sub_type = enum class Spline_sub_type : integer { CUBIC = 0, PCHIP = 1, AKIMA = 2, VANLEER = 3 };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  //! Associate a number for each type of splines implemented
  using SplineType2D = enum class SplineType2D : integer {
    BILINEAR          = 0,
    BICUBIC_CUBIC     = 1,
    BICUBIC_AKIMA     = 2,
    BICUBIC_VANLEER   = 3,
    BICUBIC_PCHIP     = 4,
    BIQUINTIC_CUBIC   = 5,
    BIQUINTIC_AKIMA   = 6,
    BIQUINTIC_VANLEER = 7,
    BIQUINTIC_PCHIP   = 8,
  };

  SplineType1D string_to_splineType1D( string_view nin );
  SplineType2D string_to_splineType2D( string_view nin );
  string       to_string( Spline_sub_type t );
  char const * to_string( SplineType1D const t );
  char const * to_string( SplineType2D const t );

  using GC_namespace::GC_type;
  using GC_namespace::GenericContainer;
  using GC_namespace::map_type;
  using GC_namespace::vec_real_type;
  using GC_namespace::vec_string_type;
  using GC_namespace::vector_type;

  /*
  //   _   _                     _ _
  //  | | | | ___ _ __ _ __ ___ (_) |_ ___
  //  | |_| |/ _ \ '__| '_ ` _ \| | __/ _ \
  //  |  _  |  __/ |  | | | | | | | ||  __/
  //  |_| |_|\___|_|  |_| |_| |_|_|\__\___|
  */

#ifdef AUTODIFF_SUPPORT
  template <typename T> inline void Hermite3( T const & x, real_type const H, T base[4] )
  {
    auto X  = x / H;
    base[1] = X * X * ( 3 - 2 * X );
    base[0] = 1 - base[1];
    base[2] = x * ( X * ( X - 2 ) + 1 );
    base[3] = x * X * ( X - 1 );
  }

  template <typename T> inline void Hermite5( T const & x, real_type const H, T base[6] )
  {
    // 1. Normalizzazione: t va da 0 a 1
    // Sostituisce le costose divisioni ripetute (invH, invH4, invH5)
    auto s = 1 / H;
    auto t = x * s;
    auto u = 1 - t;  // u è il complemento (H-x)/H

    // 2. Precalcolo delle potenze di t e u
    // Usiamo t*t invece di pow() per velocità
    auto t2 = t * t;
    auto t3 = t2 * t;

    auto u2 = u * u;
    auto u3 = u2 * u;

    // 3. Calcolo delle funzioni di base
    // Polinomio base: h00 = (1 + 3t + 6t^2) * u^3
    auto poly_t = 1 + t * ( 3 + 6 * t );
    auto poly_u = 1 + u * ( 3 + 6 * u );

    // base[0]: Valore al nodo sinistro (x=0)
    base[0] = u3 * poly_t;

    // base[1]: Valore al nodo destro (x=H)
    // Sfrutta la simmetria: scambia t con u
    base[1] = t3 * poly_u;

    // base[2]: Derivata prima al nodo sinistro
    // Scala originale: H * t * (1+3t) * u^3
    base[2] = H * t * u3 * ( 1 + 3 * t );

    // base[3]: Derivata prima al nodo destro
    // Scala originale: -H * u * (1+3u) * t^3
    base[3] = -H * u * t3 * ( 1 + 3 * u );

    // base[4]: Derivata seconda al nodo sinistro
    // Scala originale: 0.5 * H^2 * t^2 * u^3
    const real_type H2_half = ( H * H ) / 2;
    base[4]                 = H2_half * t2 * u3;

    // base[5]: Derivata seconda al nodo destro
    // Scala originale: 0.5 * H^2 * u^2 * t^3
    base[5] = H2_half * u2 * t3;
  }

#endif

  inline constexpr void Hermite3( real_type const x, real_type const H, real_type base[4] )
  {
    auto X  = x / H;
    base[1] = X * X * ( 3 - 2 * X );
    base[0] = 1 - base[1];
    base[2] = x * ( X * ( X - 2 ) + 1 );
    base[3] = x * X * ( X - 1 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  inline constexpr void Hermite3_D( real_type const x, real_type const H, real_type base_D[4] )
  {
    auto X    = x / H;
    base_D[0] = 6 * X * ( X - 1 ) / H;
    base_D[1] = -base_D[0];
    base_D[2] = ( ( 3 * X - 4 ) * X + 1 );
    base_D[3] = X * ( 3 * X - 2 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  inline constexpr void Hermite3_DD( real_type const x, real_type const H, real_type base_DD[4] )
  {
    auto X     = x / H;
    base_DD[0] = ( 12 * X - 6 ) / ( H * H );
    base_DD[1] = -base_DD[0];
    base_DD[2] = ( 6 * X - 4 ) / H;
    base_DD[3] = ( 6 * X - 2 ) / H;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  inline constexpr void Hermite3_DDD( real_type, real_type const H, real_type base_DDD[4] )
  {
    base_DDD[0] = 12 / ( H * H * H );
    base_DDD[1] = -base_DDD[0];
    base_DDD[2] = 6 / ( H * H );
    base_DDD[3] = base_DDD[2];
  }

  // --------------------------------------------------------------------------

  inline constexpr void Hermite5( real_type x, real_type H, real_type base[6] )
  {
    // 1. Normalizzazione: t va da 0 a 1
    // Sostituisce le costose divisioni ripetute (invH, invH4, invH5)
    const real_type s = 1 / H;
    const real_type t = x * s;
    const real_type u = 1 - t;  // u è il complemento (H-x)/H

    // 2. Precalcolo delle potenze di t e u
    // Usiamo t*t invece di pow() per velocità
    const real_type t2 = t * t;
    const real_type t3 = t2 * t;

    const real_type u2 = u * u;
    const real_type u3 = u2 * u;

    // 3. Calcolo delle funzioni di base
    // Polinomio base: h00 = (1 + 3t + 6t^2) * u^3
    const real_type poly_t = 1 + t * ( 3 + 6 * t );
    const real_type poly_u = 1 + u * ( 3 + 6 * u );

    // base[0]: Valore al nodo sinistro (x=0)
    base[0] = u3 * poly_t;

    // base[1]: Valore al nodo destro (x=H)
    // Sfrutta la simmetria: scambia t con u
    base[1] = t3 * poly_u;

    // base[2]: Derivata prima al nodo sinistro
    // Scala originale: H * t * (1+3t) * u^3
    base[2] = H * t * u3 * ( 1 + 3 * t );

    // base[3]: Derivata prima al nodo destro
    // Scala originale: -H * u * (1+3u) * t^3
    base[3] = -H * u * t3 * ( 1 + 3 * u );

    // base[4]: Derivata seconda al nodo sinistro
    // Scala originale: 0.5 * H^2 * t^2 * u^3
    const real_type H2_half = ( H * H ) / 2;
    base[4]                 = H2_half * t2 * u3;

    // base[5]: Derivata seconda al nodo destro
    // Scala originale: 0.5 * H^2 * u^2 * t^3
    base[5] = H2_half * u2 * t3;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  inline constexpr void Hermite5_D( real_type x, real_type H, real_type base_D[6] )
  {
    // 1. Normalizzazione
    // t = x/H, u = (H-x)/H = 1-t
    const real_type s = 1 / H;
    const real_type t = x * s;
    const real_type u = 1 - t;

    // 2. Precalcolo potenze
    const real_type t2 = t * t;
    const real_type u2 = u * u;

    // 3. Calcolo Derivate
    // Nota: La derivata rispetto a x scala di (1/H) rispetto alla derivata su t.

    // --- Basi 0 e 1 (Valori ai nodi) ---
    // Derivata del Quintic Step: 30 * t^2 * u^2 / H
    // common = 30 * t^2 * u^2 * s
    const real_type common = ( 30 * s ) * t2 * u2;

    base_D[0] = -common;
    base_D[1] = common;

    // --- Basi 2 e 3 (Derivate Prime ai nodi) ---
    // Polinomio derivato base: p(t) = (1 - 3t) * (1 + 5t)
    // base_D[2] = u^2 * p(t)
    // base_D[3] = t^2 * p(u)  <-- Simmetria perfetta

    // Espressione fattorizzata per stabilità e velocità
    base_D[2] = u2 * ( 1 - 3 * t ) * ( 1 + 5 * t );

    base_D[3] = t2 * ( 1 - 3 * u ) * ( 1 + 5 * u );

    // --- Basi 4 e 5 (Derivate Seconde ai nodi) ---
    // Polinomio derivato base: q(t) = t * (2 - 5t)
    // Scala fisica: 0.5 * H

    const real_type h_half = H / 2;

    base_D[4] = h_half * u2 * t * ( 2 - 5 * t );

    // Nota il segno meno per simmetria speculare sulla derivata
    base_D[5] = -h_half * t2 * u * ( 2 - 5 * u );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  inline constexpr void Hermite5_DD( real_type x, real_type H, real_type base_DD[6] )
  {
    // 1. Normalizzazione
    // s = 1/H, s2 = 1/H^2
    const real_type s  = 1 / H;
    const real_type s2 = s * s;

    // t = x/H (da 0 a 1), u = 1 - t
    const real_type t = x * s;
    const real_type u = 1 - t;

    // 2. Termini comuni
    // Molte basi condividono il termine t * u
    const real_type tu = t * u;

    // --- Basi 0 e 1 (Valori ai nodi) ---
    // Formula originale normalizzata: 60/H^2 * t * u * (1 - 2t)
    // Nota: (1 - 2t) è equivalente a (u - t)

    const real_type common = ( 60 * s2 ) * tu * ( u - t );

    base_DD[0] = -common;
    base_DD[1] = common;

    // --- Basi 2 e 3 (Derivate Prime ai nodi) ---
    // Scala fisica: 1/H
    // Base 2: 12/H * t * u * (5t - 3)
    // Base 3: 12/H * t * u * (5t - 2)

    const real_type k_tan = ( 12 * s ) * tu;

    base_DD[2] = k_tan * ( 5 * t - 3 );
    base_DD[3] = k_tan * ( 5 * t - 2 );

    // --- Basi 4 e 5 (Derivate Seconde ai nodi) ---
    // Scala fisica: 1 (Adimensionale rispetto a H nella derivata seconda)
    // Base 4: u * (1 - 8t + 10t^2)
    // Base 5: t * (3 - 12t + 10t^2)

    // Sfruttiamo la simmetria: Base5(t) = Base4(u)
    // Calcoliamo il polinomio 'p' usando t per base 4, e u per base 5.
    // p(v) = 1 - 8v + 10v^2

    // Base 4 (lato sinistro, usa t per il polinomio interno)
    // Nota: L'espressione u * (...) deriva dalla fattorizzazione
    base_DD[4] = u * ( 1 + t * ( 10 * t - 8 ) );

    // Base 5 (lato destro, usa u per il polinomio interno per simmetria)
    // Equivale a: t * (3 - 12t + 10t^2)
    base_DD[5] = t * ( 1 + u * ( 10 * u - 8 ) );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  inline constexpr void Hermite5_DDD( real_type x, real_type H, real_type base_DDD[6] )
  {
    // 1. Normalizzazione
    const real_type s = 1 / H;
    const real_type t = x * s;

    // t^2 per valutare le quadratiche
    const real_type t2 = t * t;

    // 2. Fattori di scala precalcolati
    // Base 0,1 scalano con 1/H^3
    // Base 2,3 scalano con 1/H^2
    // Base 4,5 scalano con 1/H
    const real_type s2 = s * s;
    const real_type s3 = s2 * s;

    // --- Basi 0 e 1 ---
    // Polinomio: 60 * (6t - 6t^2 - 1)
    const real_type k0     = 60 * s3;
    const real_type common = k0 * ( 6 * ( t - t2 ) - 1 );

    base_DDD[0] = common;
    base_DDD[1] = -common;

    // --- Basi 2 e 3 ---
    // Base 2: 12 * (16t - 15t^2 - 3)
    // Base 3: 12 * (14t - 15t^2 - 2)
    const real_type k1 = real_type( 12 ) * s2;

    // Fattorizzazione di Horner per risparmiare 1 mul: t*(16 - 15t) - 3
    base_DDD[2] = k1 * ( t * ( 16 - 15 * t ) - real_type( 3 ) );

    // Fattorizzazione: t*(14 - 15t) - 2
    base_DDD[3] = k1 * ( t * ( 14 - 15 * t ) - real_type( 2 ) );

    // --- Basi 4 e 5 ---
    // Base 4: 3 * (12t - 10t^2 - 3)
    // Base 5: 3 * (10t^2 - 8t + 1)
    const real_type k2 = 3 * s;

    base_DDD[4] = k2 * ( t * ( 12 - 10 * t ) - 3 );

    base_DDD[5] = k2 * ( t * ( 10 * t - 8 ) + 1 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  inline constexpr void Hermite5_DDDD( real_type x, real_type H, real_type base_DDDD[6] )
  {
    // 1. Normalizzazione
    // s = 1/H, t = x/H
    const real_type s = 1 / H;
    const real_type t = x * s;

    // 2. Fattori di scala
    // Base 0,1 scalano con 1/H^4
    // Base 2,3 scalano con 1/H^3
    // Base 4,5 scalano con 1/H^2
    const real_type s2 = s * s;
    const real_type s3 = s2 * s;
    const real_type s4 = s2 * s2;

    // --- Basi 0 e 1 ---
    // Originale: (360*H - 720*x) / H^5
    // Normalizzato: 360/H^4 * (1 - 2t)
    const real_type k0    = 360 * s4;
    const real_type term0 = k0 * ( 1 - 2 * t );

    base_DDDD[0] = term0;
    base_DDDD[1] = -term0;

    // --- Basi 2 e 3 ---
    // Costante comune: 24 / H^3
    const real_type k1 = 24 * s3;

    // Base 2: k1 * (8 - 15t)
    base_DDDD[2] = k1 * ( 8 - 15 * t );

    // Base 3: k1 * (7 - 15t)
    base_DDDD[3] = k1 * ( 7 - 15 * t );

    // --- Basi 4 e 5 ---
    // Costante comune: 12 / H^2
    const real_type k2 = 12 * s2;

    // Base 4: k2 * (3 - 5t)
    base_DDDD[4] = k2 * ( 3 - 5 * t );

    // Base 5: k2 * (5t - 2)
    base_DDDD[5] = k2 * ( 5 * t - 2 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  inline constexpr void Hermite5_DDDDD( real_type /*x*/, real_type H, real_type base_DDDDD[6] )
  {
    // 1. Unica divisione necessaria
    const real_type s = 1 / H;

    // 2. Calcolo delle potenze di s (moltiplicazioni)
    // s3 = 1/H^3, s4 = 1/H^4, s5 = 1/H^5
    const real_type s2 = s * s;
    const real_type s3 = s2 * s;
    const real_type s4 = s2 * s2;  // Oppure s3 * s
    const real_type s5 = s4 * s;

    // 3. Calcolo delle costanti
    // Nota: Le derivate sono costanti, quindi calcoliamo solo i valori scalari.

    // --- Basi 0 e 1 (Valore) ---
    // Derivata 5a scala con 1/H^5
    const real_type c0 = 720 * s5;

    base_DDDDD[0] = -c0;
    base_DDDDD[1] = c0;

    // --- Basi 2 e 3 (Derivata) ---
    // Derivata 5a scala con 1/H^4 (poiché H^1 * H^-5 = H^-4)
    const real_type c1 = 360 * s4;

    // Entrambe negative perché il coefficiente di x^5 nelle basi originali
    // ha lo stesso segno relativo per le derivate (antisimmetria nello sviluppo)
    base_DDDDD[2] = -c1;
    base_DDDDD[3] = -c1;

    // --- Basi 4 e 5 (Curvatura) ---
    // Derivata 5a scala con 1/H^3 (poiché H^2 * H^-5 = H^-3)
    const real_type c2 = 60 * s3;

    base_DDDDD[4] = -c2;
    base_DDDDD[5] = c2;
  }

  //!
  //! Convert polynomial defined using Hermite base
  //!
  //! \f[ p(x) = p_0 h_0(t) + p_1 h_1(t) +  p'_0 h_2(t) + p'_1 h_3(t) \f]
  //!
  //! to standard form
  //!
  //! \f[ p(x) = A t^3 + B t^2 + C t + D \f]
  //!
  //!
  //!
  static inline constexpr void Hermite3_to_poly(
    real_type   H,
    real_type   P0,
    real_type   P1,
    real_type   DP0,
    real_type   DP1,
    real_type & A,
    real_type & B,
    real_type & C,
    real_type & D )
  {
    const real_type invH  = real_type( 1 ) / H;
    const real_type slope = ( P1 - P0 ) * invH;  // Pendenza media (Secante)

    // Precalcolo inverso quadrato
    const real_type invH2 = invH * invH;

    // A = (Sum(derivs) - 2*slope) / H^2
    // Corrisponde a: (DP0 + DP1 - 2*(P1-P0)/H) / H^2
    A = ( ( DP0 + DP1 ) - 2 * slope ) * invH2;

    // B = (3*slope - (2*DP0 + DP1)) / H
    // Corrisponde a: (3*(P1-P0)/H - 2*DP0 - DP1) / H
    B = ( 3 * slope - ( 2 * DP0 + DP1 ) ) * invH;

    // C = Derivata iniziale
    C = DP0;

    // D = Posizione iniziale
    D = P0;
  }

  //!
  //! Convert polynomial defined using Hermite base
  //!
  //! \f[
  //! p(x) = p_0 h_0(t) + p_1 h_1(t) +
  //! p'_0 h_2(t) + p'_1 h_3(t) +
  //! p''_0 h_4(t) + p''_1 h_5(t)
  //! \f]
  //!
  //! to standard form
  //!
  //! \f[ p(x) = A t^5 + B t^4 + C t^3 + D t^2 + E t + F \f]
  //!
  //!
  static inline constexpr void Hermite5_to_poly(
    real_type   h,
    real_type   P0,
    real_type   P1,
    real_type   DP0,
    real_type   DP1,
    real_type   DDP0,
    real_type   DDP1,
    real_type & A,
    real_type & B,
    real_type & C,
    real_type & D,
    real_type & E,
    real_type & F )
  {
    // 1. Inverso di h (type-safe)
    const real_type invH = 1 / h;

    // 2. Precalcoli
    // Pendenza media (Velocità secante)
    const real_type slope = ( P1 - P0 ) * invH;

    // Accelerazioni dimezzate (moltiplicazione vs divisione)
    const real_type half_DDP0 = DDP0 / 2;
    const real_type half_DDP1 = DDP1 / 2;

    // 3. Potenze dell'inverso
    const real_type invH2 = invH * invH;
    const real_type invH3 = invH2 * invH;

    // 4. Calcolo Coefficienti

    // Grado 0, 1, 2
    F = P0;
    E = DP0;
    D = half_DDP0;

    // Grado 3 (C)
    // Correzione: Parentesi aggiunte per moltiplicare TUTTO per invH
    // Struttura: ( [Accelerazione] + [Velocità]*invH ) * invH
    // Risultato dimensionale: Accelerazione/H -> P/H^3
    C = ( ( half_DDP1 - 3 * half_DDP0 ) + invH * ( 10 * slope - ( 6 * DP0 + 4 * DP1 ) ) ) * invH;

    // Grado 4 (B)
    // Struttura: ( [Accelerazione] + [Velocità]*invH ) * invH^2
    // Risultato dimensionale: Accelerazione/H^2 -> P/H^4
    B = ( ( 3 * half_DDP0 - 2 * half_DDP1 ) + invH * ( ( 8 * DP0 + 7 * DP1 ) - 15 * slope ) ) * invH2;

    // Grado 5 (A)
    // Struttura: ( [Accelerazione] + [Velocità]*invH ) * invH^3
    // Risultato dimensionale: Accelerazione/H^3 -> P/H^5
    A = ( ( half_DDP1 - half_DDP0 ) + invH * ( 6 * slope - 3 * ( DP0 + DP1 ) ) ) * invH3;
  }

}  // namespace Splines

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#include "Splines/SplineUtils.hxx"
#include "Splines/SplineParametrization.hxx"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

namespace Splines
{

  /*\
   |   ____        _ _
   |  / ___| _ __ | (_)_ __   ___
   |  \___ \| '_ \| | | '_ \ / _ \
   |   ___) | |_) | | | | | |  __/
   |  |____/| .__/|_|_|_| |_|\___|
   |        |_|
  \*/
  //!
  //! Spline Management Class
  //!
  class Spline
  {
    friend class SplineSet;

  protected:
    string const m_name;
    bool         m_curve_is_closed         = false;
    bool         m_curve_can_extend        = true;
    bool         m_curve_extended_constant = false;

    integer     m_npts          = 0;
    integer     m_npts_reserved = 0;
    real_type * m_X             = nullptr;  // allocated in the derived class!
    real_type * m_Y             = nullptr;  // allocated in the derived class!

    Utils::SearchInterval<real_type, integer> m_search;

  protected:
    void copy_flags( Spline const & S )
    {
      m_curve_is_closed         = S.m_curve_is_closed;
      m_curve_can_extend        = S.m_curve_can_extend;
      m_curve_extended_constant = S.m_curve_extended_constant;
    }

  public:
    Spline( Spline const & )                   = delete;
    Spline const & operator=( Spline const & ) = delete;

    //! \name Constructors
    ///@{

    //!
    //! spline constructor
    //!
    explicit Spline( string_view const name = "Spline" ) : m_name( name )
    { m_search.setup( &m_name, &m_npts, &m_X, &m_curve_is_closed, &m_curve_can_extend ); }

    //!
    //! spline destructor
    //!
    virtual ~Spline() = default;

    ///@}

    //!
    //! \name Open/Close
    //!
    ///@{

    //!
    //! \return string with the name of the spline
    //!
    [[nodiscard]] string_view name() const noexcept { return m_name; }

    //! \return `true` if spline is a closed spline
    [[nodiscard]] bool is_closed() const noexcept { return m_curve_is_closed; }
    //!
    //! Set spline as a closed spline.
    //! When evaluated if parameter is outside the domain
    //! is wrapped cyclically before evalation.
    //!
    void make_closed() { m_curve_is_closed = true; }
    //!
    //! Set spline as an opened spline.
    //! When evaluated if parameter is outside the domain
    //! an error is produced.
    //!
    void make_opened() { m_curve_is_closed = false; }

    //!
    //! \return `true` if spline cannot extend outside interval of definition
    //!
    [[nodiscard]] bool is_bounded() const noexcept { return !m_curve_can_extend; }
    //!
    //! Set spline as unbounded.
    //! When evaluated if parameter is outside the domain
    //! an extrapolated value is used.
    //!
    void make_unbounded() { m_curve_can_extend = true; }
    //!
    //! Set spline as bounded.
    //! When evaluated if parameter is outside the domain
    //! an error is issued.
    //!
    void make_bounded() { m_curve_can_extend = false; }

    //!
    //! \return `true` if the spline extend with a constant value
    //!
    [[nodiscard]] bool is_extended_constant() const noexcept { return m_curve_extended_constant; }
    //!
    //! Set spline to extend constant.
    //! When evaluated if parameter is outside the domain
    //! the value returned is the value of the closed border.
    //!
    void make_extended_constant() { m_curve_extended_constant = true; }
    //!
    //! Set spline to extend NOT constant.
    //! When evaluated if parameter is outside the domain
    //! teh value returned is extrapolated using the last spline polynomial.
    //!
    void make_extended_not_constant() { m_curve_extended_constant = false; }

    ///@}

    //! \name Spline Data Info
    ///@{

    //!
    //! return true if spline is empty (no points)
    //!
    [[nodiscard]] bool empty() const noexcept { return m_npts == 0; }

    //!
    //! the number of support points of the spline.
    //!
    [[nodiscard]] integer num_points() const noexcept { return m_npts; }

    //!
    //! Return the pointer of values of x-nodes.
    //!
    [[nodiscard]] real_type const * x_nodes() const noexcept { return m_X; }

    //!
    //! Return the pointer of values of y-nodes.
    //!
    [[nodiscard]] real_type const * y_nodes() const noexcept { return m_Y; }

    //!
    //! the i-th node of the spline (x component).
    //!
    [[nodiscard]] real_type x_node( integer const i ) const noexcept { return m_X[i]; }

    //!
    //! the i-th node of the spline (y component).
    //!
    [[nodiscard]] real_type y_node( integer const i ) const noexcept { return m_Y[i]; }

    //!
    //! first node of the spline (x component).
    //!
    [[nodiscard]] real_type x_begin() const noexcept { return m_X[0]; }

    //!
    //! first node of the spline (y component).
    //!
    [[nodiscard]] real_type y_begin() const noexcept { return m_Y[0]; }

    //!
    //! last node of the spline (x component).
    //!
    [[nodiscard]] real_type x_end() const noexcept { return m_X[m_npts - 1]; }

    //!
    //! last node of the spline (y component).
    //!
    [[nodiscard]] real_type y_end() const noexcept { return m_Y[m_npts - 1]; }

    //!
    //! x-minumum spline value
    //!
    [[nodiscard]] real_type x_min() const noexcept { return m_X[0]; }

    //!
    //! x-maximum spline value
    //!
    [[nodiscard]] real_type x_max() const noexcept { return m_X[m_npts - 1]; }

    //!
    //! y-minumum spline value
    //!
    [[nodiscard]] real_type y_min() const
    {
      integer N = ( type() == SplineType1D::CONSTANT ) ? m_npts - 1 : m_npts;
      return *std::min_element( m_Y, m_Y + N );
    }

    //!
    //! return y-maximum spline value
    //!
    [[nodiscard]] real_type y_max() const
    {
      integer N = ( type() == SplineType1D::CONSTANT ) ? m_npts - 1 : m_npts;
      return *std::max_element( m_Y, m_Y + N );
    }

    //!
    //! Search the max and min values of `y` along the spline
    //! with the corresponding `x` position
    //!
    //! \param[out] i_min_pos interval where is the minimum
    //! \param[out] x_min_pos where is the minimum
    //! \param[out] y_min     the minimum value
    //! \param[out] i_max_pos interval where is the maximum
    //! \param[out] x_max_pos where is the maximum
    //! \param[out] y_max     the maximum value
    //!
    virtual void y_min_max(
      integer &   i_min_pos,
      real_type & x_min_pos,
      real_type & y_min,
      integer &   i_max_pos,
      real_type & x_max_pos,
      real_type & y_max ) const
    {
      i_min_pos = i_max_pos = 0;
      x_min_pos = y_min = x_max_pos = y_max = 0;
      UTILS_ERROR( "In spline: {} y_min_max not implemented\n", info() );
    }

    //!
    //! Search the max and min values of `y` along the spline
    //! with the corresponding `x` position
    //!
    //! \param[out] i_min_pos interval where is the minimum
    //! \param[out] x_min_pos where is the minimum
    //! \param[out] y_min     the minimum value
    //! \param[out] i_max_pos interval where is the maximum
    //! \param[out] x_max_pos where is the maximum
    //! \param[out] y_max     the maximum value
    //!
    virtual void y_min_max(
      vector<integer> &   i_min_pos,
      vector<real_type> & x_min_pos,
      vector<real_type> & y_min,
      vector<integer> &   i_max_pos,
      vector<real_type> & x_max_pos,
      vector<real_type> & y_max ) const
    {
      i_min_pos.clear();
      i_max_pos.clear();
      x_min_pos.clear();
      x_max_pos.clear();
      y_min.clear();
      y_max.clear();
      UTILS_ERROR( "In spline: {} y_min_max not implemented\n", info() );
    }
    ///@}

    //! \name Build
    ///@{

    //!
    //! Build a spline using data in `GenericContainer`
    //!
    void build( GenericContainer const & gc ) { setup( gc ); }

    //!
    //! Build a spline using data in a file `file_name`
    //!
    void build( string const & file_name ) { setup( file_name ); }

    //!
    //! Build a spline.
    //!
    //! \param x    vector of x-coordinates
    //! \param incx access elements as `x[0]`, `x[incx]`, `x[2*incx]`,...
    //! \param y    vector of y-coordinates
    //! \param incy access elements as `y[0]`, `y[incx]`, `y[2*incx]`,...
    //! \param n    total number of points
    //!
    template <typename VectorX, typename VectorY>
    void build( VectorX const & x, integer const incx, VectorY const & y, integer const incy, integer n )
    {
      reserve( n );
      for ( integer i = 0; i < n; ++i )
      {
        m_X[i] = x[i * incx];
        m_Y[i] = y[i * incy];
      }
      m_npts = n;
      this->build();
    }

    //!
    //! Build a spline.
    //!
    //! \param x vector of x-coordinates
    //! \param y vector of y-coordinates
    //! \param n total number of points
    //!
    template <typename VectorX, typename VectorY> void build( VectorX const & x, VectorY const & y, integer n )
    {
      reserve( n );
      for ( integer i = 0; i < n; ++i )
      {
        m_X[i] = x[i];
        m_Y[i] = y[i];
      }
      m_npts = n;
      this->build();
    }

    //!
    //! Build a spline.
    //!
    //! \param x vector of x-coordinates
    //! \param y vector of y-coordinates
    //! \param n total number of points
    //!
    template <typename VectorX, typename VectorY> void build( VectorX const & x, VectorY const & y )
    {
      integer const nx = static_cast<integer>( x.size() );
      integer const ny = static_cast<integer>( y.size() );
      integer const n  = std::min( nx, ny );
      reserve( n );
      for ( integer i = 0; i < n; ++i )
      {
        m_X[i] = x[i];
        m_Y[i] = y[i];
      }
      m_npts = n;
      this->build();
    }

    //!
    //! Build a spline using internal stored data
    //!
    virtual void build() = 0;

    //!
    //! Setup a spline using a `GenericContainer`
    //!
    //! - gc("xdata") vector with the `x` coordinate of the data
    //! - gc("ydata") vector with the `y` coordinate of the data
    //!
    virtual void setup( GenericContainer const & gc );

    //!
    //! Setup a spline using a `GenericContainer` readed from file
    //!
    //! - file_name file name of the file with data
    //!
    //! - `.json`
    //! - `.yaml` or `.yml`
    //! - `.toml`
    //!
    void setup( string const & file_name )
    {
      GenericContainer gc;
      UTILS_ASSERT( gc.from_file( file_name ), "Spline::setup( '{}' ) failed to read\n", file_name );
      setup( gc );
    }

    ///@}

    //! \name Incremental Build
    ///@{

    //!
    //! Allocate memory for `npts` points
    //!
    virtual void reserve( integer const npts ) = 0;

    //!
    //! Add a support point (x,y) to the spline.
    //!
    void push_back( real_type const x, real_type const y );

    //!
    //! Drop last inserted point of the spline.
    //!
    void drop_back()
    {
      if ( m_npts > 0 ) --m_npts;
      m_search.must_reset();
    }

    //!
    //! Delete the support points, empty the spline.
    //!
    virtual void clear() = 0;

    ///@}

    //! \name Manipulate
    ///@{

    //!
    //! change X-origin of the spline
    //!
    void set_origin( real_type const x0 ) const
    {
      real_type const Tx = x0 - m_X[0];
      real_type *     ix = m_X;
      while ( ix < m_X + m_npts ) *ix++ += Tx;
    }

    //!
    //! change X-range of the spline
    //!
    void set_range( real_type xmin, real_type xmax );

    ///@}

    //! \name Dump Data
    ///@{

    //!
    //! dump a sample of the spline
    //!
    void dump( ostream_type & s, integer const nintervals, string_view header = "x\ty" ) const;

    //!
    //! dump a sample of the spline
    //!
    void dump( string_view fname, integer const nintervals, string_view header = "x\ty" ) const
    {
      std::ofstream file{ std::string{ fname } };
      this->dump( file, nintervals, header );
      file.close();
    }

    //!
    //! Print spline coefficients
    //!
    virtual void write_to_stream( ostream_type & s ) const = 0;

    ///@}

    //! \name Evaluation
    ///@{

    [[nodiscard]] virtual real_type eval( real_type const x ) const = 0;

#ifdef AUTODIFF_SUPPORT
    //!
    //! Evaluate spline value
    //!
    [[nodiscard]] virtual autodiff::dual1st eval( autodiff::dual1st const & x ) const = 0;
    [[nodiscard]] virtual autodiff::dual2nd eval( autodiff::dual2nd const & x ) const = 0;

    // Template unificato per tutti i tipi
    template <typename T> [[nodiscard]] auto eval( T const & x ) const
    {
      if constexpr ( std::is_arithmetic<T>::value )
      {
        // Se T è un tipo numerico (int, float, double, etc.), promuovi a real_type
        return eval( static_cast<real_type>( x ) );
      }
      else
      {
        // Altrimenti deduce automaticamente il tipo duale appropriato
        return eval( autodiff::detail::to_dual( x ) );
      }
    }

    template <typename T> [[nodiscard]] auto operator()( T const & x ) const -> decltype( eval( x ) ) { return eval( x ); }
#endif

    //!
    //! First derivative
    //!
    [[nodiscard]] virtual real_type D( real_type const x ) const = 0;

    //!
    //! Second derivative
    //!
    [[nodiscard]] virtual real_type DD( real_type const x ) const = 0;

    //!
    //! Third derivative
    //!
    [[nodiscard]] virtual real_type DDD( real_type const x ) const = 0;

    //!
    //! 4th derivative
    //!
    [[nodiscard]] virtual real_type DDDD( real_type const ) const { return real_type( 0 ); }

    //!
    //! 5th derivative
    //!
    [[nodiscard]] virtual real_type DDDDD( real_type const ) const { return real_type( 0 ); }

    virtual void D( real_type const x, real_type dd[2] ) const  = 0;
    virtual void DD( real_type const x, real_type dd[3] ) const = 0;

    //!
    //! \name Batch evaluation
    //!
    //! Evaluate at every abscissa in `x`, writing to `y` (which must hold at
    //! least `x.size()` elements). These amortize the virtual dispatch of the
    //! scalar `eval`/`D`/`DD` over the whole array -- one dynamic call for the
    //! batch instead of one per point -- while producing results that are
    //! bit-identical to calling the scalar form on each element (same interval
    //! search, same per-segment evaluator). Handy for resampling/plotting.
    //!
    ///@{
    void eval( std::span<real_type const> x, std::span<real_type> y ) const
    { eval_batch( x, y, [this]( integer const ni, real_type const xx ) { return this->id_eval( ni, xx ); } ); }

    void D( std::span<real_type const> x, std::span<real_type> y ) const
    { eval_batch( x, y, [this]( integer const ni, real_type const xx ) { return this->id_D( ni, xx ); } ); }

    void DD( std::span<real_type const> x, std::span<real_type> y ) const
    { eval_batch( x, y, [this]( integer const ni, real_type const xx ) { return this->id_DD( ni, xx ); } ); }
    ///@}

    ///@}

    //!
    //! \name Evaluation Aliases
    //!
    ///@{
    //!
    //! Alias for `real_type eval( real_type x )`
    //!
    [[nodiscard]] real_type operator()( real_type const x ) const { return this->eval( x ); }
    //!
    //! Alias for `real_type D( real_type x )`
    //!
    [[nodiscard]] real_type eval_D( real_type const x ) const { return this->D( x ); }
    //!
    //! Alias for `real_type DD( real_type x )`
    //!
    [[nodiscard]] real_type eval_DD( real_type const x ) const { return this->DD( x ); }
    //!
    //! Alias for `real_type DDD( real_type x )`
    //!
    [[nodiscard]] real_type eval_DDD( real_type const x ) const { return this->DDD( x ); }
    //!
    //! Alias for `real_type DDDD( real_type x )`
    //!
    [[nodiscard]] real_type eval_DDDD( real_type const x ) const { return this->DDDD( x ); }
    //!
    //! Alias for `real_type DDDD( real_type x )`
    //!
    [[nodiscard]] real_type eval_DDDDD( real_type const x ) const { return this->DDDDD( x ); }
    ///@}

    //!
    //! \name Evaluation when segment is known
    //!
    ///@{

    //!
    //! Evaluate spline value
    //!
    [[nodiscard]] virtual real_type id_eval( integer const ni, real_type const x ) const = 0;

    //!
    //! First derivative
    //!
    [[nodiscard]] virtual real_type id_D( integer const ni, real_type const x ) const = 0;

    //!
    //! Second derivative
    //!
    [[nodiscard]] virtual real_type id_DD( integer const ni, real_type const x ) const = 0;

    //!
    //! Third derivative
    //!
    [[nodiscard]] virtual real_type id_DDD( integer const ni, real_type const x ) const = 0;

    //!
    //! 4th derivative
    //!
    [[nodiscard]] virtual real_type id_DDDD( integer const, real_type const ) const { return real_type( 0 ); }

    //!
    //! 5th derivative
    //!
    [[nodiscard]] virtual real_type id_DDDDD( integer const, real_type const ) const { return real_type( 0 ); }

    ///@}

    //! \name Get Info
    ///@{

    //!
    //! get the piecewise polinomials of the spline
    //!
    virtual integer  // order
    coeffs( real_type cfs[], real_type nodes[], bool transpose = false ) const = 0;

    //!
    //! Spline order of the piecewise polynomial
    //!
    [[nodiscard]] virtual integer order() const = 0;

    //!
    //! spline type returned as a string
    //!
    [[nodiscard]] char const * type_name() const { return to_string( type() ); }

    //!
    //! spline type returned as integer
    //!
    [[nodiscard]] virtual SplineType1D type() const = 0;

    //!
    //! String information of the kind and order of the spline
    //!
    [[nodiscard]] string info() const
    {
      string res = fmt::format( "Spline `{}` of type: {} of order: {}", m_name, type_name(), order() );
      if ( m_npts > 0 )
        res += fmt::format( "\nx_min={:.5} x_max={:.5} y_min={:.5} y_max={:.5}", x_min(), x_max(), y_min(), y_max() );
      return res;
    }

    //!
    //! Print information of the kind and order of the spline
    //!
    void info( ostream_type & stream ) const { stream << this->info() << '\n'; }

    ///@}

#ifdef SPLINES_BACK_COMPATIBILITY
    void pushBack( real_type x, real_type y ) { push_back( x, y ); }
    void dropBack() { drop_back(); }
    void setOrigin( real_type x0 ) { set_origin( x0 ); }
    void setRange( real_type xmin, real_type xmax ) { set_range( xmin, xmax ); }
    void writeToStream( ostream_type & s ) const { write_to_stream( s ); }
#endif

  private:
    //!
    //! Shared batch loop: for each abscissa run the interval search once and
    //! apply `f( segment_index, wrapped_abscissa )`. Kept in one place so the
    //! `eval`/`D`/`DD` batch overloads stay bit-identical to their scalar
    //! forms (identical search, identical per-segment evaluator).
    //!
    template <typename PointEval>
    void eval_batch( std::span<real_type const> x, std::span<real_type> y, PointEval f ) const
    {
      UTILS_ASSERT(
        y.size() >= x.size(),
        "Spline[{}]::batch eval: output span ({}) smaller than input ({})",
        m_name, y.size(), x.size()
      );
      std::size_t const n{ x.size() };
      if ( n == 0 ) return;

      std::pair<integer, real_type> res{ 0, real_type( 0 ) };

      // Both branches below are bit-identical to the scalar path; the only
      // difference is speed, so we pick per input shape. `find()` is the
      // dominant per-point cost, so for monotone non-decreasing input we skip
      // it whenever a query stays in the previous segment (the common
      // resample/plot case -- measured 3-4x faster). For unsorted input that
      // reuse guard almost never hits and its branch just mispredicts, so we
      // fall back to a plain find-per-point loop that matches the scalar cost.
      // The one up-front is_sorted() scan is O(n), branch-predictable, and
      // negligible next to n spline evaluations.
      if ( std::is_sorted( x.begin(), x.end() ) )
      {
        // Interval-reuse fast path. The reuse guard is the strict half-open
        // [X[cur], X[cur+1]) -- exactly the segment find() returns for an
        // in-range x (largest k with X[k] <= x) -- so at an exact upper knot,
        // or anything out of the current segment (out-of-domain, closed-curve
        // wrapping), we fall back to find() and use its possibly-wrapped
        // abscissa. benchmarks/bench_eval.cc asserts the resulting
        // equivalence across knots, out-of-domain points, and closed splines.
        integer cur{ 0 };
        bool    have_cur{ false };
        for ( std::size_t i{ 0 }; i < n; ++i )
        {
          real_type const xi{ x[i] };
          if ( have_cur && xi >= m_X[cur] && xi < m_X[cur + 1] )
          {
            y[i] = f( cur, xi );
          }
          else
          {
            res.first  = 0;
            res.second = xi;
            m_search.find( res );
            cur        = res.first;
            have_cur   = true;
            y[i]       = f( res.first, res.second );
          }
        }
      }
      else
      {
        for ( std::size_t i{ 0 }; i < n; ++i )
        {
          res.first  = 0;
          res.second = x[i];
          m_search.find( res );
          y[i] = f( res.first, res.second );
        }
      }
    }
  };

  //!
  //! Compute curvature of a planar curve
  //! Formula: k = (v x a) / |v|^3
  //!
  [[nodiscard]] real_type curvature( real_type const s, Spline const & X, Spline const & Y );

  //!
  //! Compute curvature derivative
  //! Formula vettoriale compatta: k' = ( |v|^2 (v x j) - 3 (v x a) (v . a) ) / |v|^5
  //!
  [[nodiscard]] real_type curvature_D( real_type const s, Spline const & X, Spline const & Y );

  //!
  //! Compute curvature second derivative
  //! Sostituisce l'espansione polinomiale manuale con algebra vettoriale esatta
  //! Riduce drasticamente il numero di moltiplicazioni e migliora la precisione.
  //!
  [[nodiscard]] real_type curvature_DD( real_type const s, Spline const & X, Spline const & Y );

}  // namespace Splines

#include "Splines/SplineCubicBase.hxx"
#include "Splines/SplineCubic.hxx"

#include "Splines/SplineBuild.hxx"

#include "Splines/SplineHermite.hxx"
#include "Splines/SplineAkima.hxx"
#include "Splines/SplineVanLeer.hxx"
#include "Splines/SplineConstant.hxx"
#include "Splines/SplineLinear.hxx"
#include "Splines/SplinePchip.hxx"
#include "Splines/SplineQuinticBase.hxx"
#include "Splines/SplineQuintic.hxx"

#include "Splines/SplineSurf.hxx"

#include "Splines/SplineBilinear.hxx"

#include "Splines/SplineBiCubicBase.hxx"
#include "Splines/SplineBiCubic.hxx"

#include "Splines/SplineBiQuinticBase.hxx"
#include "Splines/SplineBiQuintic.hxx"

#include "Splines/SplineVec.hxx"
#include "Splines/SplineSet.hxx"
#include "Splines/Splines1D.hxx"
#include "Splines/Splines2D.hxx"

#include "Splines/Splines1Dblend.hxx"
#include "Splines/Splines2Dblend.hxx"

namespace SplinesLoad
{

  using Splines::AkimaSpline;
  using Splines::ConstantSpline;
  using Splines::CubicSpline;
  using Splines::CubicSplineBase;
  using Splines::LinearSpline;
  using Splines::PchipSpline;
  using Splines::QuinticSpline;
  using Splines::Spline;
  using Splines::Spline1D;
  using Splines::VanLeerSpline;

  using Splines::BiCubicSpline;
  using Splines::BilinearSpline;
  using Splines::BiQuinticSpline;
  using Splines::Spline2D;

  using Splines::SplineSet;
  using Splines::SplineVec;

  using Splines::SplineType1D;
  using Splines::SplineType2D;
}  // namespace SplinesLoad

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
#ifdef __clang__
#pragma clang diagnostic pop
#endif

#endif
