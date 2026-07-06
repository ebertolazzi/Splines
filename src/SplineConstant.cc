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
 |    ____                _              _       ____        _ _
 |   / ___|___  _ __  ___| |_ __ _ _ __ | |_ ___/ ___| _ __ | (_)_ __   ___
 |  | |   / _ \| '_ \/ __| __/ _` | '_ \| __/ __\___ \| '_ \| | | '_ \ / _ \
 |  | |__| (_) | | | \__ \ || (_| | | | | |_\__ \___) | |_) | | | | | |  __/
 |   \____\___/|_| |_|___/\__\__,_|_| |_|\__|___/____/| .__/|_|_|_| |_|\___|
 |                                                    |_|
\*/

#include "Splines/Splines.hh"

namespace Splines
{

  void ConstantSpline::write_to_stream( ostream_type & s ) const
  {
    integer const nseg = m_npts > 0 ? m_npts - 1 : 0;
    for ( integer i = 0; i < nseg; ++i )
      fmt::print( s, "segment N. {:4} X:[{:.5},{:.5}] Y:{:.5}\n", i, m_X[i], m_X[i + 1], m_Y[i] );
  }

  // --------------------------- VIRTUALS -----------------------------------

  void ConstantSpline::reserve( integer const npts )
  {
    if ( m_external_alloc && npts <= m_npts_reserved )
    {
      // nothing to do!, already allocated
    }
    else
    {
      m_mem_constant.reallocate( 2 * npts );
      m_npts_reserved  = npts;
      m_external_alloc = false;
      m_X              = m_mem_constant( npts );
      m_Y              = m_mem_constant( npts );
    }
    m_npts = 0;
  }

  void ConstantSpline::clear()
  {
    if ( !m_external_alloc ) m_mem_constant.free();
    m_npts = m_npts_reserved = 0;
    m_external_alloc         = false;
    m_X = m_Y = nullptr;
  }

  void ConstantSpline::setup( GenericContainer const & gc )
  {
    /*
    // gc["xdata"]
    // gc["ydata"]
    //
    */
    string const where = fmt::format( "ConstantSpline[{}]::setup( gc ):", m_name );

    std::set<std::string> keywords;
    for ( auto const & pair : gc.get_map( where ) ) { keywords.insert( pair.first ); }
    keywords.erase( "spline_type" );

    GenericContainer const & gc_x = gc( "xdata", where );
    keywords.erase( "xdata" );
    GenericContainer const & gc_y = gc( "ydata", where );
    keywords.erase( "ydata" );

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

    this->build( x, y );
  }


  void ConstantSpline::y_min_max(
    integer &   i_min_pos,
    real_type & x_min_pos,
    real_type & y_min,
    integer &   i_max_pos,
    real_type & x_max_pos,
    real_type & y_max ) const
  {
    UTILS_ASSERT( m_npts > 0, "ConstantSpline[{}]::y_min_max() empty spline!", m_name );
    // find max min alongh the nodes
    i_min_pos = i_max_pos = 0;
    x_min_pos = x_max_pos = m_X[0];
    y_min = y_max = m_Y[0];
    for ( integer i = 1; i < m_npts - 1; ++i )
    {
      real_type const & P1 = m_Y[i];
      if ( P1 > y_max )
      {
        y_max     = P1;
        x_max_pos = m_X[i];
        i_max_pos = i;
      }
      else if ( P1 < y_min )
      {
        y_min     = P1;
        x_min_pos = m_X[i];
        i_min_pos = i;
      }
    }
  }

  void ConstantSpline::y_min_max(
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
    UTILS_ASSERT( m_npts > 0, "ConstantSpline[{}]::y_min_max() empty spline!", m_name );
    // find max min along the nodes
    for ( integer i = 1; i < m_npts - 1; ++i )
    {
      real_type const & P0 = m_Y[i - 1];
      real_type const & P1 = m_Y[i];
      real_type const & P2 = m_Y[i + 1];
      if ( P1 > P0 && P1 > P2 )
      {
        y_max.emplace_back( P1 );
        x_max_pos.emplace_back( m_X[i] );
        i_max_pos.emplace_back( i );
      }
      else if ( P1 < P0 && P1 < P2 )
      {
        y_min.emplace_back( P1 );
        x_min_pos.emplace_back( m_X[i] );
        i_min_pos.emplace_back( i );
      }
    }
  }

  void ConstantSpline::copy_spline( ConstantSpline const & S )
  {
    // 1. Protezione contro l'auto-assegnamento (Self-assignment check)
    // Se stiamo copiando l'oggetto su se stesso, non facciamo nulla.
    if ( this == &S ) return;

    // 2. Gestione casi vuoti
    if ( S.m_npts <= 0 )
    {
      m_npts = 0;
      return;
    }

    // 3. Allocazione memoria
    // Assumiamo che reserve gestisca la riallocazione solo se necessario.
    this->reserve( S.m_npts );
    m_npts = S.m_npts;

    // 4. Copia rapida della memoria (Block Copy)
    std::memcpy( m_X, S.m_X, m_npts * sizeof( real_type ) );

    // 5. Copia sicura dei valori Y (Intervalli)
    if ( m_npts > 1 ) std::memcpy( m_Y, S.m_Y, ( m_npts - 1 ) * sizeof( real_type ) );

    // 6. Copia flag e metadati
    copy_flags( S );
  }

}  // namespace Splines

//
// EOF: SplineConstant.cc
//
