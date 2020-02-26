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
/****************************************************************************\
Copyright (c) 2016, Enrico Bertolazzi
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the FreeBSD Project.
\****************************************************************************/
#pragma once

#ifndef SPLINES_CONFIG_HH
#define SPLINES_CONFIG_HH

// Uncomment this if you want to enable debugging
// #define DEBUG

#include "GenericContainer.hh"

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <utility>      // std::pair
#include <algorithm>
#include <thread>
#include <mutex>

//
// file: Splines
//
// if C++ < C++11 define nullptr
// select computer architecture
#if defined(_WIN32) || defined(WIN32) || defined(_WIN64) || defined(WIN64)
  // windows architecture
  #define SPLINES_OS_WINDOWS 1
  #if defined(_MSC_VER) && _MSC_VER < 1900
    #error "Splines need a C++11 compiler"
  #endif
  #if defined(__cplusplus) && __cplusplus < 199711L
    #error "Splines need a C++11 compiler"
  #endif
  #ifdef _MSC_VER
    #include <math.h>
  #endif
#else
  #if defined(__APPLE__) && defined(__MACH__)
    // osx architecture
    #define SPLINES_OS_OSX 1
  #elif defined(__unix__)
    // linux architecture
    #define SPLINES_OS_LINUX 1
  #endif
  #include <cmath>
  #include <cfloat>
  #if __cplusplus < 199711L
    #error "Splines need a C++11 compiler"
  #endif
#endif

#ifdef DEBUG
  #define SPLINE_CHECK_NAN( PTR, MSG, DIM ) Splines::checkNaN( PTR, MSG, DIM )
#else
  #define SPLINE_CHECK_NAN( PTR, MSG, DIM )
#endif

#endif
