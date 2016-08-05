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

#include "Splines.hh"
#include <limits>
#include <cmath>

/**
 * 
 */

// interface with GenericContainer
#ifdef SPLINES_USE_GENERIC_CONTAINER

namespace Splines {

  //! Evaluate all the splines at `x` and fill a map of values in a GenericContainer
  void
  SplineSet::eval( valueType x, GenericContainer & gc ) const {
    map_type & vals = gc.set_map() ;
    for ( map<string,sizeType>::const_iterator s_to_pos = header_to_position.begin() ;
          s_to_pos != header_to_position.end() ; ++s_to_pos ) {
      vals[s_to_pos->first] = splines[s_to_pos->second]->eval(x) ;
    }
  }

  //! Evaluate all the splines at `x` values contained in vec and fill a map of vector in a GenericContainer
  void
  SplineSet::eval( vec_real_type const & vec, GenericContainer & gc ) const {
    sizeType npts = sizeType(vec.size()) ;
    map_type & vals = gc.set_map() ;
    for ( map<string,sizeType>::const_iterator s_to_pos = header_to_position.begin() ;
          s_to_pos != header_to_position.end() ; ++s_to_pos ) {
      vec_real_type & v = vals[s_to_pos->first].set_vec_real(npts) ;
      Spline const * p_spl = splines[s_to_pos->second] ;
      for ( sizeType i = 0 ; i < npts ; ++i ) v[i] = p_spl->eval(vec[i]) ;
    }
  }

  //! Evaluate all the splines at `x` and fill a map of values in a GenericContainer with keys in `columns`
  void
  SplineSet::eval( valueType x, vec_string_type const & columns, GenericContainer & gc ) const {
    map_type & vals = gc.set_map() ;
    for ( vec_string_type::const_iterator is = columns.begin() ;
          is != columns.end() ; ++is ) {
      Spline const * p_spl = getSpline( is->c_str() ) ;
      vals[*is] = p_spl->eval(x) ;
    }
  }

  //! Evaluate all the splines at `x` values contained in vec and fill a map of vector in a GenericContainer with keys in `columns`
  void
  SplineSet::eval( vec_real_type const & vec, vec_string_type const & columns, GenericContainer & gc ) const {
    sizeType npts = sizeType(vec.size()) ;
    map_type & vals = gc.set_map() ;
    for ( vec_string_type::const_iterator is = columns.begin() ;
          is != columns.end() ; ++is ) {
      Spline const * p_spl = getSpline( is->c_str() ) ;
      vec_real_type & v = vals[*is].set_vec_real(npts) ;
      for ( sizeType i = 0 ; i < npts ; ++i ) v[i] = p_spl->eval(vec[i]) ;
    }
  }

  //! Evaluate all the splines at `zeta` and fill a map of values in a GenericContainer and `indep` as independent spline
  void
  SplineSet::eval2( valueType zeta, sizeType indep, GenericContainer & gc ) const {
    map_type & vals = gc.set_map() ;
    valueType x ;
    intersect( indep, zeta, x ) ;
    for ( map<string,sizeType>::const_iterator s_to_pos = header_to_position.begin() ;
          s_to_pos != header_to_position.end() ; ++s_to_pos ) {
      vals[s_to_pos->first] = splines[s_to_pos->second]->eval(x) ;
    }
  }

  //! Evaluate all the splines at `zeta` values contained in vec and fill a map of vector in a GenericContainer and `indep` as independent spline
  void
  SplineSet::eval2( vec_real_type const & zetas, sizeType indep, GenericContainer & gc ) const {
    sizeType npts = sizeType(zetas.size()) ;
    map_type & vals = gc.set_map() ;

    // preallocation
    for ( map<string,sizeType>::const_iterator s_to_pos = header_to_position.begin() ;
          s_to_pos != header_to_position.end() ; ++s_to_pos )
      vals[s_to_pos->first].set_vec_real(npts) ;

    for ( sizeType i = 0 ; i < npts ; ++i ) {
      valueType x ;
      intersect( indep, zetas[i], x ) ;
      for ( map<string,sizeType>::const_iterator s_to_pos = header_to_position.begin() ;
            s_to_pos != header_to_position.end() ; ++s_to_pos ) {
        vec_real_type & v = vals[s_to_pos->first].get_vec_real() ;
        v[i] = splines[s_to_pos->second]->eval(x) ;
      }
    }
  }

  //! Evaluate all the splines at `zeta` and fill a map of values in a GenericContainer with keys in `columns` and `indep` as independent spline
  void
  SplineSet::eval2( valueType zeta, sizeType indep, vec_string_type const & columns, GenericContainer & gc ) const {
    map_type & vals = gc.set_map() ;
    valueType x ;
    intersect( indep, zeta, x ) ;
    for ( vec_string_type::const_iterator is = columns.begin() ;
          is != columns.end() ; ++is ) {
      Spline const * p_spl = getSpline( is->c_str() ) ;
      vals[*is] = p_spl->eval(x) ;
    }
  }

  //! Evaluate all the splines at `zeta` values contained in vec and fill a map of vector in a GenericContainer with keys in `columns` and `indep` as independent spline
  void
  SplineSet::eval2( vec_real_type const & zetas, sizeType indep, vec_string_type const & columns, GenericContainer & gc ) const {
    sizeType npts = sizeType(zetas.size()) ;
    map_type & vals = gc.set_map() ;

    // preallocation
    for ( vec_string_type::const_iterator is = columns.begin() ;
          is != columns.end() ; ++is )
      vals[*is].set_vec_real(npts) ;

    for ( sizeType i = 0 ; i < npts ; ++i ) {
      valueType x ;
      intersect( indep, zetas[i], x ) ;
      for ( vec_string_type::const_iterator is = columns.begin() ;
            is != columns.end() ; ++is ) {
        vec_real_type & v = vals[*is].get_vec_real() ;
        Spline const * p_spl = getSpline( is->c_str() ) ;
        v[i] = p_spl->eval(x) ;
      }
    }
  }

  /*
  //    __ _          _         _           _            _   _
  //   / _(_)_ __ ___| |_    __| | ___ _ __(_)_   ____ _| |_(_)_   _____
  //  | |_| | '__/ __| __|  / _` |/ _ \ '__| \ \ / / _` | __| \ \ / / _ \
  //  |  _| | |  \__ \ |_  | (_| |  __/ |  | |\ V / (_| | |_| |\ V /  __/
  //  |_| |_|_|  |___/\__|  \__,_|\___|_|  |_| \_/ \__,_|\__|_| \_/ \___|
  */
  void
  SplineSet::eval_D( valueType x, GenericContainer & gc ) const {
    map_type & vals = gc.set_map() ;
    for ( map<string,sizeType>::const_iterator s_to_pos = header_to_position.begin() ;
          s_to_pos != header_to_position.end() ; ++s_to_pos ) {
      vals[s_to_pos->first] = splines[s_to_pos->second]->eval_D(x) ;
    }
  }

  //! Evaluate all the splines at `x` values contained in vec and fill a map of vector in a GenericContainer
  void
  SplineSet::eval_D( vec_real_type const & vec, GenericContainer & gc ) const {
    sizeType npts = sizeType(vec.size()) ;
    map_type & vals = gc.set_map() ;
    for ( map<string,sizeType>::const_iterator s_to_pos = header_to_position.begin() ;
          s_to_pos != header_to_position.end() ; ++s_to_pos ) {
      vec_real_type & v = vals[s_to_pos->first].set_vec_real(npts) ;
      Spline const * p_spl = splines[s_to_pos->second] ;
      for ( sizeType i = 0 ; i < npts ; ++i ) v[i] = p_spl->eval_D(vec[i]) ;
    }
  }

  //! Evaluate all the splines at `x` and fill a map of values in a GenericContainer with keys in `columns`
  void
  SplineSet::eval_D( valueType x, vec_string_type const & columns, GenericContainer & gc ) const {
    map_type & vals = gc.set_map() ;
    for ( vec_string_type::const_iterator is = columns.begin() ;
          is != columns.end() ; ++is ) {
      Spline const * p_spl = getSpline( is->c_str() ) ;
      vals[*is] = p_spl->eval_D(x) ;
    }
  }

  //! Evaluate all the splines at `x` values contained in vec and fill a map of vector in a GenericContainer with keys in `columns`
  void
  SplineSet::eval_D( vec_real_type const & vec, vec_string_type const & columns, GenericContainer & gc ) const {
    sizeType npts = sizeType(vec.size()) ;
    map_type & vals = gc.set_map() ;
    for ( vec_string_type::const_iterator is = columns.begin() ;
          is != columns.end() ; ++is ) {
      Spline const * p_spl = getSpline( is->c_str() ) ;
      vec_real_type & v = vals[*is].set_vec_real(npts) ;
      for ( sizeType i = 0 ; i < npts ; ++i ) v[i] = p_spl->eval_D(vec[i]) ;
    }
  }

  //! Evaluate all the splines at `zeta` and fill a map of values in a GenericContainer and `indep` as independent spline
  void
  SplineSet::eval2_D( valueType zeta, sizeType indep, GenericContainer & gc ) const {
    map_type & vals = gc.set_map() ;
    valueType x ;
    intersect( indep, zeta, x ) ;
    for ( map<string,sizeType>::const_iterator s_to_pos = header_to_position.begin() ;
          s_to_pos != header_to_position.end() ; ++s_to_pos ) {
      vals[s_to_pos->first] = splines[s_to_pos->second]->eval_D(x) ;
    }
  }

  //! Evaluate all the splines at `zeta` values contained in vec and fill a map of vector in a GenericContainer and `indep` as independent spline
  void
  SplineSet::eval2_D( vec_real_type const & zetas, sizeType indep, GenericContainer & gc ) const {
    sizeType npts = sizeType(zetas.size()) ;
    map_type & vals = gc.set_map() ;

    // preallocation
    for ( map<string,sizeType>::const_iterator s_to_pos = header_to_position.begin() ;
          s_to_pos != header_to_position.end() ; ++s_to_pos )
      vals[s_to_pos->first].set_vec_real(npts) ;

    for ( sizeType i = 0 ; i < npts ; ++i ) {
      valueType x ;
      intersect( indep, zetas[i], x ) ;
      for ( map<string,sizeType>::const_iterator s_to_pos = header_to_position.begin() ;
            s_to_pos != header_to_position.end() ; ++s_to_pos ) {
        vec_real_type & v = vals[s_to_pos->first].get_vec_real() ;
        v[i] = splines[s_to_pos->second]->eval_D(x) ;
      }
    }
  }

  //! Evaluate all the splines at `zeta` and fill a map of values in a GenericContainer with keys in `columns` and `indep` as independent spline
  void
  SplineSet::eval2_D( valueType zeta, sizeType indep, vec_string_type const & columns, GenericContainer & gc ) const {
    map_type & vals = gc.set_map() ;
    valueType x ;
    intersect( indep, zeta, x ) ;
    for ( vec_string_type::const_iterator is = columns.begin() ;
          is != columns.end() ; ++is ) {
      Spline const * p_spl = getSpline( is->c_str() ) ;
      vals[*is] = p_spl->eval_D(x) ;
    }
  }

  //! Evaluate all the splines at `zeta` values contained in vec and fill a map of vector in a GenericContainer with keys in `columns` and `indep` as independent spline
  void
  SplineSet::eval2_D( vec_real_type const & zetas, sizeType indep, vec_string_type const & columns, GenericContainer & gc ) const {
    sizeType npts = sizeType(zetas.size()) ;
    map_type & vals = gc.set_map() ;

    // preallocation
    for ( vec_string_type::const_iterator is = columns.begin() ;
          is != columns.end() ; ++is )
      vals[*is].set_vec_real(npts) ;

    for ( sizeType i = 0 ; i < npts ; ++i ) {
      valueType x ;
      intersect( indep, zetas[i], x ) ;
      for ( vec_string_type::const_iterator is = columns.begin() ;
            is != columns.end() ; ++is ) {
        vec_real_type & v = vals[*is].get_vec_real() ;
        Spline const * p_spl = getSpline( is->c_str() ) ;
        v[i] = p_spl->eval_D(x) ;
      }
    }
  }

  /*                                _       _           _            _   _
  //   ___  ___  ___ ___  _ __   __| |   __| | ___ _ __(_)_   ____ _| |_(_)_   _____
  //  / __|/ _ \/ __/ _ \| '_ \ / _` |  / _` |/ _ \ '__| \ \ / / _` | __| \ \ / / _ \
  //  \__ \  __/ (_| (_) | | | | (_| | | (_| |  __/ |  | |\ V / (_| | |_| |\ V /  __/
  //  |___/\___|\___\___/|_| |_|\__,_|  \__,_|\___|_|  |_| \_/ \__,_|\__|_| \_/ \___|
  */
  void
  SplineSet::eval_DD( valueType x, GenericContainer & gc ) const {
    map_type & vals = gc.set_map() ;
    for ( map<string,sizeType>::const_iterator s_to_pos = header_to_position.begin() ;
          s_to_pos != header_to_position.end() ; ++s_to_pos ) {
      vals[s_to_pos->first] = splines[s_to_pos->second]->eval_DD(x) ;
    }
  }

  //! Evaluate all the splines at `x` values contained in vec and fill a map of vector in a GenericContainer
  void
  SplineSet::eval_DD( vec_real_type const & vec, GenericContainer & gc ) const {
    sizeType npts = sizeType(vec.size()) ;
    map_type & vals = gc.set_map() ;
    for ( map<string,sizeType>::const_iterator s_to_pos = header_to_position.begin() ;
          s_to_pos != header_to_position.end() ; ++s_to_pos ) {
      vec_real_type & v = vals[s_to_pos->first].set_vec_real(npts) ;
      Spline const * p_spl = splines[s_to_pos->second] ;
      for ( sizeType i = 0 ; i < npts ; ++i ) v[i] = p_spl->eval_DD(vec[i]) ;
    }
  }

  //! Evaluate all the splines at `x` and fill a map of values in a GenericContainer with keys in `columns`
  void
  SplineSet::eval_DD( valueType x, vec_string_type const & columns, GenericContainer & gc ) const {
    map_type & vals = gc.set_map() ;
    for ( vec_string_type::const_iterator is = columns.begin() ;
          is != columns.end() ; ++is ) {
      Spline const * p_spl = getSpline( is->c_str() ) ;
      vals[*is] = p_spl->eval_DD(x) ;
    }
  }

  //! Evaluate all the splines at `x` values contained in vec and fill a map of vector in a GenericContainer with keys in `columns`
  void
  SplineSet::eval_DD( vec_real_type const & vec, vec_string_type const & columns, GenericContainer & gc ) const {
    sizeType npts = sizeType(vec.size()) ;
    map_type & vals = gc.set_map() ;
    for ( vec_string_type::const_iterator is = columns.begin() ;
          is != columns.end() ; ++is ) {
      Spline const * p_spl = getSpline( is->c_str() ) ;
      vec_real_type & v = vals[*is].set_vec_real(npts) ;
      for ( sizeType i = 0 ; i < npts ; ++i ) v[i] = p_spl->eval_DD(vec[i]) ;
    }
  }

  //! Evaluate all the splines at `zeta` and fill a map of values in a GenericContainer and `indep` as independent spline
  void
  SplineSet::eval2_DD( valueType zeta, sizeType indep, GenericContainer & gc ) const {
    map_type & vals = gc.set_map() ;
    valueType x ;
    intersect( indep, zeta, x ) ;
    for ( map<string,sizeType>::const_iterator s_to_pos = header_to_position.begin() ;
          s_to_pos != header_to_position.end() ; ++s_to_pos ) {
      vals[s_to_pos->first] = splines[s_to_pos->second]->eval_DD(x) ;
    }
  }

  //! Evaluate all the splines at `zeta` values contained in vec and fill a map of vector in a GenericContainer and `indep` as independent spline
  void
  SplineSet::eval2_DD( vec_real_type const & zetas, sizeType indep, GenericContainer & gc ) const {
    sizeType npts = sizeType(zetas.size()) ;
    map_type & vals = gc.set_map() ;

    // preallocation
    for ( map<string,sizeType>::const_iterator s_to_pos = header_to_position.begin() ;
          s_to_pos != header_to_position.end() ; ++s_to_pos )
      vals[s_to_pos->first].set_vec_real(npts) ;

    for ( sizeType i = 0 ; i < npts ; ++i ) {
      valueType x ;
      intersect( indep, zetas[i], x ) ;
      for ( map<string,sizeType>::const_iterator s_to_pos = header_to_position.begin() ;
            s_to_pos != header_to_position.end() ; ++s_to_pos ) {
        vec_real_type & v = vals[s_to_pos->first].get_vec_real() ;
        v[i] = splines[s_to_pos->second]->eval_DD(x) ;
      }
    }
  }

  //! Evaluate all the splines at `zeta` and fill a map of values in a GenericContainer with keys in `columns` and `indep` as independent spline
  void
  SplineSet::eval2_DD( valueType zeta, sizeType indep, vec_string_type const & columns, GenericContainer & gc ) const {
    map_type & vals = gc.set_map() ;
    valueType x ;
    intersect( indep, zeta, x ) ;
    for ( vec_string_type::const_iterator is = columns.begin() ;
          is != columns.end() ; ++is ) {
      Spline const * p_spl = getSpline( is->c_str() ) ;
      vals[*is] = p_spl->eval_DD(x) ;
    }
  }

  //! Evaluate all the splines at `zeta` values contained in vec and fill a map of vector in a GenericContainer with keys in `columns` and `indep` as independent spline
  void
  SplineSet::eval2_DD( vec_real_type const & zetas, sizeType indep, vec_string_type const & columns, GenericContainer & gc ) const {
    sizeType npts = sizeType(zetas.size()) ;
    map_type & vals = gc.set_map() ;

    // preallocation
    for ( vec_string_type::const_iterator is = columns.begin() ;
          is != columns.end() ; ++is )
      vals[*is].set_vec_real(npts) ;

    for ( sizeType i = 0 ; i < npts ; ++i ) {
      valueType x ;
      intersect( indep, zetas[i], x ) ;
      for ( vec_string_type::const_iterator is = columns.begin() ;
            is != columns.end() ; ++is ) {
        vec_real_type & v = vals[*is].get_vec_real() ;
        Spline const * p_spl = getSpline( is->c_str() ) ;
        v[i] = p_spl->eval_DD(x) ;
      }
    }
  }

  /*
  //   _   _     _         _       _           _            _   _
  //  | |_| |__ (_)_ __ __| |   __| | ___ _ __(_)_   ____ _| |_(_)_   _____
  //  | __| '_ \| | '__/ _` |  / _` |/ _ \ '__| \ \ / / _` | __| \ \ / / _ \
  //  | |_| | | | | | | (_| | | (_| |  __/ |  | |\ V / (_| | |_| |\ V /  __/
  //   \__|_| |_|_|_|  \__,_|  \__,_|\___|_|  |_| \_/ \__,_|\__|_| \_/ \___|
  */
  void
  SplineSet::eval_DDD( valueType x, GenericContainer & gc ) const {
    map_type & vals = gc.set_map() ;
    for ( map<string,sizeType>::const_iterator s_to_pos = header_to_position.begin() ;
          s_to_pos != header_to_position.end() ; ++s_to_pos ) {
      vals[s_to_pos->first] = splines[s_to_pos->second]->eval_DDD(x) ;
    }
  }

  //! Evaluate all the splines at `x` values contained in vec and fill a map of vector in a GenericContainer
  void
  SplineSet::eval_DDD( vec_real_type const & vec, GenericContainer & gc ) const {
    sizeType npts = sizeType(vec.size()) ;
    map_type & vals = gc.set_map() ;
    for ( map<string,sizeType>::const_iterator s_to_pos = header_to_position.begin() ;
          s_to_pos != header_to_position.end() ; ++s_to_pos ) {
      vec_real_type & v = vals[s_to_pos->first].set_vec_real(npts) ;
      Spline const * p_spl = splines[s_to_pos->second] ;
      for ( sizeType i = 0 ; i < npts ; ++i ) v[i] = p_spl->eval_DDD(vec[i]) ;
    }
  }

  //! Evaluate all the splines at `x` and fill a map of values in a GenericContainer with keys in `columns`
  void
  SplineSet::eval_DDD( valueType x, vec_string_type const & columns, GenericContainer & gc ) const {
    map_type & vals = gc.set_map() ;
    for ( vec_string_type::const_iterator is = columns.begin() ;
          is != columns.end() ; ++is ) {
      Spline const * p_spl = getSpline( is->c_str() ) ;
      vals[*is] = p_spl->eval_DDD(x) ;
    }
  }

  //! Evaluate all the splines at `x` values contained in vec and fill a map of vector in a GenericContainer with keys in `columns`
  void
  SplineSet::eval_DDD( vec_real_type const & vec, vec_string_type const & columns, GenericContainer & gc ) const {
    sizeType npts = sizeType(vec.size()) ;
    map_type & vals = gc.set_map() ;
    for ( vec_string_type::const_iterator is = columns.begin() ;
          is != columns.end() ; ++is ) {
      Spline const * p_spl = getSpline( is->c_str() ) ;
      vec_real_type & v = vals[*is].set_vec_real(npts) ;
      for ( sizeType i = 0 ; i < npts ; ++i ) v[i] = p_spl->eval_DDD(vec[i]) ;
    }
  }

  //! Evaluate all the splines at `zeta` and fill a map of values in a GenericContainer and `indep` as independent spline
  void
  SplineSet::eval2_DDD( valueType zeta, sizeType indep, GenericContainer & gc ) const {
    map_type & vals = gc.set_map() ;
    valueType x ;
    intersect( indep, zeta, x ) ;
    for ( map<string,sizeType>::const_iterator s_to_pos = header_to_position.begin() ;
          s_to_pos != header_to_position.end() ; ++s_to_pos ) {
      vals[s_to_pos->first] = splines[s_to_pos->second]->eval_DDD(x) ;
    }
  }

  //! Evaluate all the splines at `zeta` values contained in vec and fill a map of vector in a GenericContainer and `indep` as independent spline
  void
  SplineSet::eval2_DDD( vec_real_type const & zetas, sizeType indep, GenericContainer & gc ) const {
    sizeType npts = sizeType(zetas.size()) ;
    map_type & vals = gc.set_map() ;

    // preallocation
    for ( map<string,sizeType>::const_iterator s_to_pos = header_to_position.begin() ;
          s_to_pos != header_to_position.end() ; ++s_to_pos )
      vals[s_to_pos->first].set_vec_real(npts) ;

    for ( sizeType i = 0 ; i < npts ; ++i ) {
      valueType x ;
      intersect( indep, zetas[i], x ) ;
      for ( map<string,sizeType>::const_iterator s_to_pos = header_to_position.begin() ;
            s_to_pos != header_to_position.end() ; ++s_to_pos ) {
        vec_real_type & v = vals[s_to_pos->first].get_vec_real() ;
        v[i] = splines[s_to_pos->second]->eval_DDD(x) ;
      }
    }
  }

  //! Evaluate all the splines at `zeta` and fill a map of values in a GenericContainer with keys in `columns` and `indep` as independent spline
  void
  SplineSet::eval2_DDD( valueType zeta, sizeType indep, vec_string_type const & columns, GenericContainer & gc ) const {
    map_type & vals = gc.set_map() ;
    valueType x ;
    intersect( indep, zeta, x ) ;
    for ( vec_string_type::const_iterator is = columns.begin() ;
          is != columns.end() ; ++is ) {
      Spline const * p_spl = getSpline( is->c_str() ) ;
      vals[*is] = p_spl->eval_DDD(x) ;
    }
  }

  //! Evaluate all the splines at `zeta` values contained in vec and fill a map of vector in a GenericContainer with keys in `columns` and `indep` as independent spline
  void
  SplineSet::eval2_DDD( vec_real_type const & zetas, sizeType indep, vec_string_type const & columns, GenericContainer & gc ) const {
    sizeType npts = sizeType(zetas.size()) ;
    map_type & vals = gc.set_map() ;

    // preallocation
    for ( vec_string_type::const_iterator is = columns.begin() ;
          is != columns.end() ; ++is )
      vals[*is].set_vec_real(npts) ;

    for ( sizeType i = 0 ; i < npts ; ++i ) {
      valueType x ;
      intersect( indep, zetas[i], x ) ;
      for ( vec_string_type::const_iterator is = columns.begin() ;
            is != columns.end() ; ++is ) {
        vec_real_type & v = vals[*is].get_vec_real() ;
        Spline const * p_spl = getSpline( is->c_str() ) ;
        v[i] = p_spl->eval_DDD(x) ;
      }
    }
  }
}

#endif
