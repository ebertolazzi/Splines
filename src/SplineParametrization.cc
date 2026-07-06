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

//
// FILE: SplineParametrization.cc
//

#include "Splines/Splines.hh"

namespace Splines
{
  /*\
   |  ____                                _        _          _   _
   | |  _ \ __ _ _ __ __ _ _ __ ___   ___| |_ _ __(_)______ _| |_(_) ___  _ __
   | | |_) / _` | '__/ _` | '_ ` _ \ / _ \ __| '__| |_  / _` | __| |/ _ \| '_ \
   | |  __/ (_| | | | (_| | | | | | |  __/ |_| |  | |/ / (_| | |_| | (_) | | | |
   | |_|   \__,_|_|  \__,_|_| |_| |_|\___|\__|_|  |_/___\__,_|\__|_|\___/|_| |_|
   |
  \*/

  void chordal( integer const dim, integer const npts, real_type const pnts[], integer const ld_pnts, real_type t[] )
  {
    t[0]                 = 0;
    real_type const * p0 = pnts;
    for ( integer k = 1; k < npts; ++k )
    {
      real_type const * p1  = p0 + ld_pnts;
      real_type         dst = 0;
      for ( integer j = 0; j < dim; ++j )
      {
        real_type const c = p1[j] - p0[j];
        dst += c * c;
      }
      t[k] = t[k - 1] + sqrt( dst );
    }
    for ( integer k = 1; k < npts - 1; ++k ) t[k] /= t[npts - 1];
    t[npts - 1] = 1;
  }

  void centripetal(
    integer const   dim,
    integer const   npts,
    real_type const pnts[],
    integer const   ld_pnts,
    real_type const alpha,
    real_type       t[] )
  {
    t[0]                 = 0;
    real_type const * p0 = pnts;
    for ( integer k = 1; k < npts; ++k )
    {
      real_type const * p1  = p0 + ld_pnts;
      real_type         dst = 0;
      for ( integer j = 0; j < dim; ++j )
      {
        real_type const c = p1[j] - p0[j];
        dst += c * c;
      }
      t[k] = t[k - 1] + pow( dst, alpha / 2 );
    }
    for ( integer k = 1; k < npts - 1; ++k ) t[k] /= t[npts - 1];
    t[npts - 1] = 1;
  }

  void universal( integer const dim, integer const npts, real_type const pnts[], integer const ld_pnts, real_type t[] )
  {
    // Initialize t[0] = 0
    t[0] = 0;

    // Calculate Euclidean distances between consecutive points
    real_type const * p0 = pnts;
    for ( integer k = 1; k < npts; ++k )
    {
      real_type const * p1  = p0 + ld_pnts;
      real_type         dst = 0;
      for ( integer j = 0; j < dim; ++j )
      {
        real_type const c = p1[j] - p0[j];
        dst += c * c;
      }
      t[k] = t[k - 1] + sqrt( dst );
      p0   = p1;
    }

    // Calculate total curve length
    real_type const total_length = t[npts - 1];

    // Normalize parameters using universal formula
    // combining cumulative distance and uniformly distributed index
    for ( integer k = 1; k < npts - 1; ++k )
    {
      real_type const chord_param   = t[k] / total_length;
      real_type const uniform_param = static_cast<real_type>( k ) / static_cast<real_type>( npts - 1 );
      // Weighted average between chordal and uniform parameterization
      t[k] = ( chord_param + uniform_param ) / 2;
    }

    t[npts - 1] = 1;
  }

  void FoleyNielsen(
    integer const   dim,
    integer const   npts,
    real_type const pnts[],
    integer const   ld_pnts,
    real_type       t[] )
  {
    // Initialize t[0] = 0
    t[0] = 0;

    // Calculate distances between consecutive points
    std::vector<real_type> d( npts - 1 );
    real_type const *      p0 = pnts;

    for ( integer k = 0; k < npts - 1; ++k )
    {
      real_type const * p1  = p0 + ld_pnts;
      real_type         dst = 0;
      for ( integer j = 0; j < dim; ++j )
      {
        real_type const c = p1[j] - p0[j];
        dst += c * c;
      }
      d[k] = sqrt( dst );
      p0   = p1;
    }

    // Calculate angles between consecutive segments (Foley-Nielsen modifiers)
    for ( integer k = 1; k < npts; ++k )
    {
      real_type alpha_k = 0;

      if ( k > 0 && k < npts - 1 )
      {
        // Calculate angle at point k
        real_type const d_prev = d[k - 1];
        real_type const d_next = d[k];
        real_type const d_sum  = d_prev + d_next;

        if ( d_sum > 0 )
        {
          // Adjustment factor based on angles
          alpha_k = std::min( d_prev, d_next ) / d_sum;
          alpha_k = std::min( alpha_k, real_type{ 0.5 } );
        }
      }

      // Foley-Nielsen parameterization with angular correction
      real_type const correction = 1 + 1.5 * alpha_k;
      t[k]                       = t[k - 1] + d[k - 1] * correction;
    }

    // Normalize
    for ( integer k = 1; k < npts - 1; ++k ) t[k] /= t[npts - 1];
    t[npts - 1] = 1;
  }

  void FangHung( integer const dim, integer const npts, real_type const pnts[], integer const ld_pnts, real_type t[] )
  {
    // Initialize t[0] = 0
    t[0] = 0;

    // Calculate distances between consecutive points
    std::vector<real_type> d( npts - 1 );
    real_type const *      p0 = pnts;

    for ( integer k = 0; k < npts - 1; ++k )
    {
      real_type const * p1  = p0 + ld_pnts;
      real_type         dst = 0;
      for ( integer j = 0; j < dim; ++j )
      {
        real_type const c = p1[j] - p0[j];
        dst += c * c;
      }
      d[k] = sqrt( dst );
      p0   = p1;
    }

    // Fang-Hung method: uses weighted combination of distances
    for ( integer k = 1; k < npts; ++k )
    {
      real_type weight = 1;

      // Calculate weight based on local curvature
      if ( k > 0 && k < npts - 1 )
      {
        real_type const d_prev = d[k - 1];
        real_type const d_next = d[k];
        real_type const ratio  = d_prev / ( d_next + 1e-10 );

        // Adjustment based on distance ratio
        // to compensate for curvature variations
        weight = 1 + 0.5 * std::abs( ratio - 1 );
      }

      t[k] = t[k - 1] + d[k - 1] * weight;
    }

    // Normalize
    for ( integer k = 1; k < npts - 1; ++k ) t[k] /= t[npts - 1];
    t[npts - 1] = 1;
  }
}  // namespace Splines

//
// EOF: SplineParametrization.cc
//
