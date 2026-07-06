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
   |    _   _                     _ _       ____        _ _
   |   | | | | ___ _ __ _ __ ___ (_) |_ ___/ ___| _ __ | (_)_ __   ___
   |   | |_| |/ _ \ '__| '_ ` _ \| | __/ _ \___ \| '_ \| | | '_ \ / _ \
   |   |  _  |  __/ |  | | | | | | | ||  __/___) | |_) | | | | | |  __/
   |   |_| |_|\___|_|  |_| |_| |_|_|\__\___|____/| .__/|_|_|_| |_|\___|
   |                                             |_|
  \*/

  void HermiteSpline::setup( GenericContainer const & gc )
  {
    /*
    // gc["xdata"]
    // gc["ydata"]
    //
    */
    string const where = fmt::format( "HermiteSpline[{}]::setup( gc ):", m_name );

    std::set<std::string> keywords;
    for ( auto const & pair : gc.get_map( where ) ) { keywords.insert( pair.first ); }
    keywords.erase( "spline_type" );

    GenericContainer const & gc_x = gc( "xdata", where );
    keywords.erase( "xdata" );
    GenericContainer const & gc_y = gc( "ydata", where );
    keywords.erase( "ydata" );
    GenericContainer const & gc_yp = gc( "ypdata", where );
    keywords.erase( "ypdata" );

    vec_real_type x, y, yp;
    {
      string const ff = fmt::format( "{}, field `xdata'", where );
      gc_x.copyto_vec_real( x, ff );
    }
    {
      string const ff = fmt::format( "{}, field `ydata'", where );
      gc_y.copyto_vec_real( y, ff );
    }
    {
      string const ff = fmt::format( "{}, field `ypdata'", where );
      gc_yp.copyto_vec_real( yp, ff );
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

    this->build( x, y, yp );
  }

}  // namespace Splines
