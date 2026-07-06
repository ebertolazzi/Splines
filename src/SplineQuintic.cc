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
 |    ___        _       _   _      ____        _ _
 |   / _ \ _   _(_)_ __ | |_(_) ___/ ___| _ __ | (_)_ __   ___
 |  | | | | | | | | '_ \| __| |/ __\___ \| '_ \| | | '_ \ / _ \
 |  | |_| | |_| | | | | | |_| | (__ ___) | |_) | | | | | |  __/
 |   \__\_\\__,_|_|_| |_|\__|_|\___|____/| .__/|_|_|_| |_|\___|
 |                                       |_|
 |
\*/

#include "Splines/Splines.hh"

namespace Splines
{

  void QuinticSpline::build()
  {
    string msg = fmt::format( "QuinticSpline[{}]::build():", m_name );
    UTILS_ASSERT( m_npts > 1, "{} npts = {} not enought points\n", msg, m_npts );
    Utils::check_NaN( m_X, msg + " X", m_npts, __LINE__, __FILE__ );
    Utils::check_NaN( m_Y, msg + " Y", m_npts, __LINE__, __FILE__ );

    integer ibegin = 0;
    integer iend   = 0;
    do
    {
      // cerca intervallo monotono strettamente crescente
      for ( ++iend; iend < m_npts && m_X[iend - 1] < m_X[iend]; ++iend ) {}
      Quintic_build2( m_sub_type, m_X + ibegin, m_Y + ibegin, m_Yp + ibegin, m_Ypp + ibegin, iend - ibegin );
      ibegin = iend;
    } while ( iend < m_npts );

    Utils::check_NaN( m_Yp, msg + " Yp", m_npts, __LINE__, __FILE__ );
    Utils::check_NaN( m_Ypp, msg + " Ypp", m_npts, __LINE__, __FILE__ );
    m_search.must_reset();
  }

  //! Build a Monotone quintic spline from data from `gc`
  void QuinticSpline::setup( GenericContainer const & gc )
  {
    string const where = fmt::format( "QuinticSpline[{}]::setup( gc ):", m_name );

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

    if ( gc.exists( "spline_sub_type" ) )
    {
      string_view st = gc.get_map_string( "spline_sub_type", where );
      keywords.erase( "spline_sub_type" );

      // OTTIMIZZAZIONE 5: Switch più efficiente del if-else chain
      switch ( st[0] )  // First character hint
      {
        case 'c':  // cubic
          if ( st == "cubic" )
            m_sub_type = Spline_sub_type::CUBIC;
          else
            goto unknown_type;
          break;
        case 'p':  // pchip
          if ( st == "pchip" )
            m_sub_type = Spline_sub_type::PCHIP;
          else
            goto unknown_type;
          break;
        case 'a':  // akima
          if ( st == "akima" )
            m_sub_type = Spline_sub_type::AKIMA;
          else
            goto unknown_type;
          break;
        case 'v':  // van Leer
          if ( st == "vanleer" )
            m_sub_type = Spline_sub_type::VANLEER;
          else
            goto unknown_type;
          break;
        default:
        unknown_type:
          UTILS_ERROR( "{} unknow sub type: {}\n", where, st );
      }
    }
    else
    {
      UTILS_WARNING( false, "{}, missing field `spline_sub_type` using `cubic` as default value\n", where );
    }

    UTILS_WARNING(
      keywords.empty(),
      "{}: unused keys\n{}\n",
      where,
      [&keywords]() -> string
      {
        string res;
        res.reserve( keywords.size() * 10 );  // OTTIMIZZAZIONE 6: Pre-allocazione
        for ( auto const & it : keywords )
        {
          res += it;
          res += ' ';
        }
        return res;
      }() );

    this->build( x, y );
  }
}  // namespace Splines

//
// EOF: SplineQuintic.cc
//
