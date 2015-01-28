/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2015                                                      |
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

#ifndef _MALLOC_HH
#define _MALLOC_HH

#include <string>
#include <iostream>
#include <sstream>
#include <stdexcept>

/*
//   __  __       _ _            
//  |  \/  | __ _| | | ___   ___ 
//  | |\/| |/ _` | | |/ _ \ / __|
//  | |  | | (_| | | | (_) | (__ 
//  |_|  |_|\__,_|_|_|\___/ \___|
*/

//! Allocate memory
template <typename T>
class Malloc {
  typedef T                 valueType           ;
  typedef valueType*        valuePointer        ;  
  typedef const valueType*  valueConstPointer   ; 
  typedef valueType&        valueReference      ;
  typedef const valueType&  valueConstReference ;

  typedef long              indexType           ;
  typedef indexType*        indexPointer        ;  
  typedef const indexType*  indexConstPointer   ; 
  typedef indexType&        indexReference      ;
  typedef const indexType&  indexConstReference ;

private:

  std::string  _name ;
  size_t       numTotValues ;
  size_t       numTotReserved ;
  size_t       numAllocated ;
  valuePointer pMalloc ;

  Malloc(Malloc<T> const &) ; // block copy constructor
  Malloc<T> const & operator = (Malloc<T> &) const ; // block copy constructor

public:

  //! malloc object constructor
  explicit
  Malloc( std::string const & __name )
  : _name(__name)
  , numTotValues(0)
  , numTotReserved(0)
  , numAllocated(0)
  , pMalloc(nullptr)
  {}

  //! malloc object destructor
  ~Malloc()
  { free() ; }

  //! allocate memory for `n` objects
  void
  allocate( size_t n ) {
    try {
      if ( n > numTotReserved ) {
        delete [] pMalloc ;
        numTotValues   = n ;
        numTotReserved = n + (n>>3) ; // 12% more values
        pMalloc = new T[numTotReserved] ;
      }
    }
    catch ( std::exception const & exc ) {
      std::cerr << "Memory allocation failed: " << exc.what()
                << "\nTry to allocate " << n << " bytes for " << _name
                << '\n' ;
      exit(0) ;
    }
    catch (...) {
      std::cerr << "Malloc allocation failed for " << _name << ": memory exausted\n"
                << "Requesting " << n << " blocks\n";
      exit(0) ;
    }
    numTotValues = n ;
    numAllocated = 0 ;
  }

  //! free memory
  void
  free(void) {
    if ( pMalloc != 0 ) {
      delete [] pMalloc ;
      numTotValues   = 0 ;
      numTotReserved = 0 ;
      numAllocated   = 0 ;
      pMalloc = nullptr ;
    }
  }

  //! number of objects allocated 
  indexType size(void) const { return numTotValues ; }

  //! get pointer of allocated memory for `sz` objets
  T * operator () ( size_t sz ) {
    size_t offs = numAllocated ;
    numAllocated += sz ;
    if ( numAllocated > numTotValues ) {
      std::ostringstream ost ;
      ost << "\nMalloc<" << _name << ">::operator () (" << sz << ") -- Malloc EXAUSTED\n"
          << "request = " << numAllocated << " > " << numTotValues << " = available\n" ;
      throw std::runtime_error(ost.str()) ;
    }
    return pMalloc + offs ;
  }
} ;

#endif
