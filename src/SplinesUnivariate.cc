/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 1998                                                      |
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
#include "SplinesUtils.hh"
#include <cmath>
#include <iomanip>

/**
 * 
 */

namespace Splines {
  //! Associate a number for each type of splines implemented

  /*\
   |       _               _    _   _       _   _
   |   ___| |__   ___  ___| | _| \ | | __ _| \ | |
   |  / __| '_ \ / _ \/ __| |/ /  \| |/ _` |  \| |
   | | (__| | | |  __/ (__|   <| |\  | (_| | |\  |
   |  \___|_| |_|\___|\___|_|\_\_| \_|\__,_|_| \_|
  \*/
  /*!
   | check if the vector `pv` os size `DIM` contains only regular floats.
   | If not an error is issued
  \*/
  void
  checkNaN(
    real_type const pv[],
    char      const v_name[],
    integer         DIM
  ) {
    for ( size_t i = 0; i < size_t(DIM); ++i ) {
      if ( isNaN(pv[i]) ) {
        std::ostringstream ost;
        ost << "\nfound NaN at " << v_name << "[" << i << "]";
        throw std::runtime_error(ost.str());
      }
      if ( isInfinite(pv[i]) ) {
        std::ostringstream ost;
        ost << "\nfound Infinity at " << v_name << "[" << i << "]";
        throw std::runtime_error(ost.str());
      }
    }
  }
}
