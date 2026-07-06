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

/*\
 |      _    _    _                   ____        _ _
 |     / \  | | _(_)_ __ ___   __ _  / ___| _ __ | (_)_ __   ___
 |    / _ \ | |/ / | '_ ` _ \ / _` | \___ \| '_ \| | | '_ \ / _ \
 |   / ___ \|   <| | | | | | | (_| |  ___) | |_) | | | | | |  __/
 |  /_/   \_\_|\_\_|_| |_| |_|\__,_| |____/| .__/|_|_|_| |_|\___|
 |                                         |_|
\*/

namespace Splines
{

  void AkimaSpline::build()
  {
    string const msg = fmt::format( "AkimaSpline[{}]::build():", m_name );
    UTILS_ASSERT( m_npts > 1, "{} npts = {} not enought points\n", msg, m_npts );
    Utils::check_NaN( m_X, msg + " X ", m_npts, __LINE__, __FILE__ );
    Utils::check_NaN( m_Y, msg + " Y ", m_npts, __LINE__, __FILE__ );
    integer ibegin = 0;
    integer iend   = 0;

    Vec work( m_npts );

    do
    {
      // cerca intervallo monotono strettamente crescente
      for ( ++iend; iend < m_npts && m_X[iend - 1] < m_X[iend]; ++iend ) {}
      Akima_build( m_X + ibegin, m_Y + ibegin, m_Yp + ibegin, work.data(), iend - ibegin );
      ibegin = iend;
    } while ( iend < m_npts );

    Utils::check_NaN( m_Yp, msg + " Yp", m_npts, __LINE__, __FILE__ );
    m_search.must_reset();
  }

  void AkimaSpline::setup( GenericContainer const & gc )
  {
    /*
    // gc["xdata"]
    // gc["ydata"]
    //
    */
    string const where = fmt::format( "AkimaSpline[{}]::setup():", m_name );

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
}  // namespace Splines

//
// EOF: SplineAkima.cc
//
