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
 |      Universita` degli Studi di Trento                                   |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

/*\
 | file: SplinesCinterface.h
 */

/*!
 * \addtogroup Splines C interface
 */
/* @{ */
#pragma once

#ifndef SPLINES_C_INTERFACE_H
#define SPLINES_C_INTERFACE_H

#ifdef __cplusplus
extern "C" {
#endif

  /*! Create a new ``Spline`` object ``id`` */
  int
  SPLINE_new(
    char const * id,
    char const * type
  );

  /*! Select a ``Spline`` object ``id`` */
  int SPLINE_select( char const * id );

  /*! Delete the ``Spline`` object ``id`` */
  int SPLINE_delete( char const * id );

  /*! Print the actual ``Spline`` */
  int SPLINE_print();

  /*! Get type of actual pointed element of ``Spline`` */
  char const * SPLINE_get_type_name();

  /*! Get pointer to the internal ``Spline`` object 'id' */
  void * SPLINE_mem_ptr( char const * id );

  /*! Set actual pointed element of ``Spline`` to an empty spline. */
  int SPLINE_init();

  /*! Push `(x,y)` interpolation point to a spline `bool` */
  int SPLINE_push( double x, double y );

  /*! Build a Spline with the points previously inserted */
  int SPLINE_build();

  /*! Build a Spline with the points passed as arguments */
  int
  SPLINE_build2(
    double const * x,
    double const * y,
    int            n
  );

  /*! Evaluate spline at ``x`` */
  double SPLINE_eval( double x );

  /*! Evaluate spline first derivative at ``x`` */
  double SPLINE_eval_D( double x );

  /*! Evaluate spline second derivative at ``x`` */
  double SPLINE_eval_DD( double x );

  /*! Evaluate spline third derivative at ``x`` */
  double SPLINE_eval_DDD( double x );

  /*! Evaluate spline 4th derivative at ``x`` */
  double SPLINE_eval_DDDD( double x );

  /*! Evaluate spline 5th derivative at ``x`` */
  double SPLINE_eval_DDDDD( double x );

#ifdef __cplusplus
}
#endif

#endif

/* @} */

/*
// eof: SplineCinterface.hh
*/
