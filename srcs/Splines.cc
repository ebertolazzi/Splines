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
/*!
 *
 * \date     October 28, 2015
 * \version  5.1
 * \note     first release Jan 12, 1998
 *
 * \author   Enrico Bertolazzi
 *
 * \par      Affiliation:
 *           Department of Industrial Engineering<br>
 *           University of Trento<br>
 *           Via Sommarive 9, I -- 38123 Trento, Italy <br>
 *           enrico.bertolazzi@unitn.it
 *
 */

#include "Splines.hh"

//! Various kind of splines
namespace Splines {

  //! Check if cubic spline with this data is monotone, -2 non monotone data, -1 no, 0 yes, 1 strictly monotone
  indexType
  checkCubicSplineMonotonicity( valueType const X[],
                                valueType const Y[],
                                valueType const Yp[],
                                indexType       npts ) {
    // check monotonicity of data: (assuming X monotone)
    indexType flag = 1 ;
    for ( indexType i = 1 ; i < npts ; ++i ) {
      if ( Y[i-1] > Y[i] ) return -2 ; // non monotone data
      if ( Y[i-1] == Y[i] && X[i-1] < X[i] ) flag = 0 ; // non strict monotone
    }
    // pag 146 Methods of Shape-Preserving Spline Approximation, K
    for ( indexType i = 1 ; i < npts ; ++i ) {
      if ( X[i] <= X[i-1] ) continue ; // skip duplicate points
      valueType dd = (Y[i]-Y[i-1])/(X[i]-X[i-1]) ;
      valueType m0 = Yp[i-1]/dd ;
      valueType m1 = Yp[i]/dd ;
      if ( m0 < 0 || m1 < 0 ) return -1 ; // non monotone
      if ( m0 <= 3 && m1 <= 3 ) {
        if ( flag > 0 && i > 1      && (m0 == 0 || m0 == 3) ) flag = 0 ;
        if ( flag > 0 && i < npts-1 && (m1 == 0 || m1 == 3) ) flag = 0 ;
      } else {
        valueType tmp1 = 2*m0+m1-3 ;
        valueType tmp2 = 2*(m0+m1-2) ;
        valueType tmp3 = m0*tmp2-(tmp1*tmp1)  ;
        if ( tmp2 >= 0 ) {
          if ( tmp3 < 0 ) return -1 ; // non monotone spline
        } else {
          if ( tmp3 > 0 ) return -1 ;
        }
        if ( tmp3 == 0 ) flag = 0 ;
      }
    }
    return flag ; // passed all check
  }
}
