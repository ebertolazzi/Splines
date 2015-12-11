/****************************************************************************\
  Copyright (c) Enrico Bertolazzi 2014
  All Rights Reserved.

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation;

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
\****************************************************************************/

#include "Splines.hh"
#include "mex.h"
#include <map>
#include <string>

using namespace std ;
using Splines::valueType ;
using Splines::sizeType ;

#define arg_spline_name   prhs[0]
#define arg_spline_type   prhs[1]

#define out_z             plhs[0]
#define out_Dx            plhs[1]
#define out_Dy            plhs[2]
#define out_Dxx           plhs[3]
#define out_Dxy           plhs[4]
#define out_Dyy           plhs[5]

#define ASSERT(COND,MSG)                 \
  if ( !(COND) ) {                       \
    std::ostringstream ost ;             \
    ost << "spline2d: " << MSG << '\n' ; \
    mexErrMsgTxt(ost.str().c_str()) ;    \
  }

static
void
mexErrorMessage() {
  // Check for proper number of arguments, etc
  mexErrMsgTxt(
"%======================================================================%\n"
"% spline2d:  Compute clothoid parameters along a Clothoid curve        %\n"
"%                                                                      %\n"
"% USAGE: [X,Y,TH,K] = evalClothoid( x0, y0, theta0, k, dk, s ) ;       %\n"
"%                                                                      %\n"
"% On input:                                                            %\n"
"%                                                                      %\n"
"%  x0, y0 = coodinate of initial point                                 %\n"
"%  theta0 = orientation (angle) of the clothoid at initial point       %\n"
"%  k      = curvature at initial point                                 %\n"
"%  dk     = derivative of curvature respect to arclength               %\n"
"%  s      = vector of curvilinear coordinate where to compute clothoid %\n"
"%                                                                      %\n"
"% On output:                                                           %\n"
"%                                                                      %\n"
"%  X     = X coordinate of the points of the clothoid at s coordinate  %\n"
"%  Y     = Y coordinate of the points of the clothoid at s coordinate  %\n"
"%  TH    = angle of the clothoid at s coordinate                       %\n"
"%  K     = curvature of the clothoid at s coordinate                   %\n"
"%======================================================================%\n"
"%                                                                      %\n"
"%  Autor: Enrico Bertolazzi                                            %\n"
"%         Department of Industrial Engineering                         %\n"
"%         University of Trento                                         %\n"
"%         enrico.bertolazzi@unitn.it                                   %\n"
"%                                                                      %\n"
"%======================================================================%\n" ) ;
}

extern "C"
void
mexFunction( int nlhs, mxArray       *plhs[],
             int nrhs, mxArray const *prhs[] ) {

  static std::map<string,Splines::SplineSurf*> pSplines ;

  // the first argument must be a string
  
  ASSERT( mxIsChar(arg_spline_name), "First argument must be a string" ) ;
  string sname = mxArrayToString(arg_spline_name) ;

  if ( mxIsChar(arg_spline_type) ) { // constructor
    string tname = mxArrayToString(arg_spline_type) ;
    
    Splines::SplineSurf * p_spline = nullptr ;

    if ( tname == "bilinear" ) {
      p_spline = new Splines::BilinearSpline() ;
    } else if ( tname == "bicubic" ) {
      p_spline = new Splines::BiCubicSpline() ;
    } else if ( tname == "akima" ) {
      p_spline = new Splines::Akima2Dspline() ;
    } else if ( tname == "biquintic" ) {
      p_spline = new Splines::BiQuinticSpline() ;
    } else {
      ASSERT(false,"Second argument must be a string: 'bilinear', 'bicubic', 'akima', 'biquintic'" ) ;
    }
    pSplines[sname] = p_spline ;

    mxArray const * arg_x = prhs[2] ;
    mxArray const * arg_y = prhs[3] ;
    mxArray const * arg_z = prhs[4] ;

    mwSize number_of_dimensions = mxGetNumberOfDimensions(arg_x) ;
    ASSERT( number_of_dimensions == 2, "Expect vector as third argument" ) ;
    mwSize const * dims_x = mxGetDimensions(arg_x) ;
    ASSERT( dims_x[0] == 1 || dims_x[1] == 1, "Expect (1 x n or n x 1) matrix as third argument, found " << dims_x[0] << " x " << dims_x[1] ) ;
    double const * x = mxGetPr(arg_x) ;
    sizeType  nx = dims_x[0]*dims_x[1] ;

    number_of_dimensions = mxGetNumberOfDimensions(arg_y) ;
    ASSERT( number_of_dimensions == 2, "Expect vector as 4th argument" ) ;
    mwSize const * dims_y = mxGetDimensions(arg_y) ;
    ASSERT( dims_y[0] == 1 || dims_y[1] == 1, "Expect (1 x n or n x 1) matrix as 4th argument, found " << dims_y[0] << " x " << dims_y[1] ) ;
    double const * y = mxGetPr(arg_y) ;
    sizeType  ny = dims_y[0]*dims_y[1] ;

    number_of_dimensions = mxGetNumberOfDimensions(arg_z) ;
    ASSERT( number_of_dimensions == 2, "Expect matrix as 5th argument" ) ;
    mwSize const * dims_z = mxGetDimensions(arg_z) ;
    ASSERT( dims_z[0] == nx || dims_z[1] == ny, "Expect (" << nx << " x " << ny << " matrix as 5th argument, found " << dims_z[0] << " x " << dims_z[1] ) ;
    double const * z = mxGetPr(arg_z) ;
    sizeType ldZ = nx ;

    //for ( mwSize i = 0 ; i < nx ; ++i ) mexPrintf("x[%d] = %g\n", i, x[i]) ;
    //for ( mwSize j = 0 ; j < ny ; ++j ) mexPrintf("y[%d] = %g\n", j, y[j]) ;
    //for ( mwSize j = 0 ; j < nx*ny ; ++j ) mexPrintf("z[%d] = %g\n", j, z[j]) ;

    p_spline -> build ( x, 1, y, 1, z, ldZ, nx, ny, true, false ) ;

  } else { // eval
    // check if spline exists
    std::map<string,Splines::SplineSurf*>::iterator it = pSplines.find(sname) ;
    ASSERT( it != pSplines.end(), "Spline: ``" << sname << "'' not defined" );
    
    // valutazione
    if ( nrhs == 2 ) { // spline2d('nome',XY)
      //arg_xy
      mxArray const * arg_xy = prhs[1] ;

      mwSize number_of_dimensions = mxGetNumberOfDimensions(arg_xy) ;
      ASSERT( number_of_dimensions == 2, "Expect (2 x n) matrix as second argument" ) ;
      mwSize const * dims = mxGetDimensions(arg_xy) ;
      ASSERT( dims[0] == 2, "Expect (2 x n) matrix as second argument, found " << dims[0] << " x " << dims[1] ) ;

      valueType * x = mxGetPr(arg_xy) ;
      
      ASSERT( nlhs >= 1, "Expect at least one output argument" ) ;
      
      Splines::SplineSurf * p_spline = it->second ;
      //mexPrintf("SPLINE: %s\n", p_spline->type_name()) ;
      //mexPrintf("NX: %d\n", p_spline->numPointX()) ;
      //mexPrintf("NY: %d\n", p_spline->numPointY()) ;

      if ( nlhs == 1 ) {
        //mexPrintf("dims: %d %d\n", dims[0], dims[1]) ;
        out_z = mxCreateNumericMatrix(1, dims[1], mxDOUBLE_CLASS, mxREAL);
        valueType * z = mxGetPr(out_z) ;
        for ( mwSize i = 0 ; i < dims[1] ; ++i ) {
          z[i] = (*p_spline)(x[2*i],x[2*i+1]) ;
          //mexPrintf("%g = S(%g,%g)\n", z[i], x[2*i],x[2*i+1]) ;
        }
      } else if ( nlhs == 3 ) {
        out_z  = mxCreateNumericMatrix(1, dims[1], mxDOUBLE_CLASS, mxREAL);
        out_Dx = mxCreateNumericMatrix(1, dims[1], mxDOUBLE_CLASS, mxREAL);
        out_Dy = mxCreateNumericMatrix(1, dims[1], mxDOUBLE_CLASS, mxREAL);
        valueType * z  = mxGetPr(out_z) ;
        valueType * dx = mxGetPr(out_Dx) ;
        valueType * dy = mxGetPr(out_Dy) ;
        for ( mwSize i = 0 ; i < dims[1] ; ++i ) {
          valueType vals[3] ;
          p_spline->D(x[2*i],x[2*i+1],vals) ;
          z[i] = vals[0] ;
          dx[i] = vals[1] ;
          dy[i] = vals[2] ;
        }
      } else if ( nlhs == 6 ) {
        out_z   = mxCreateNumericMatrix(1, dims[1], mxDOUBLE_CLASS, mxREAL);
        out_Dx  = mxCreateNumericMatrix(1, dims[1], mxDOUBLE_CLASS, mxREAL);
        out_Dy  = mxCreateNumericMatrix(1, dims[1], mxDOUBLE_CLASS, mxREAL);
        out_Dxx = mxCreateNumericMatrix(1, dims[1], mxDOUBLE_CLASS, mxREAL);
        out_Dxy = mxCreateNumericMatrix(1, dims[1], mxDOUBLE_CLASS, mxREAL);
        out_Dyy = mxCreateNumericMatrix(1, dims[1], mxDOUBLE_CLASS, mxREAL);
        valueType * z   = mxGetPr(out_z) ;
        valueType * dx  = mxGetPr(out_Dx) ;
        valueType * dy  = mxGetPr(out_Dy) ;
        valueType * dxx = mxGetPr(out_Dxx) ;
        valueType * dxy = mxGetPr(out_Dxy) ;
        valueType * dyy = mxGetPr(out_Dyy) ;
        for ( mwSize i = 0 ; i < dims[1] ; ++i ) {
          valueType vals[6] ;
          p_spline->DD(x[2*i],x[2*i+1],vals) ;
          z[i]   = vals[0] ;
          dx[i]  = vals[1] ;
          dy[i]  = vals[2] ;
          dxx[i] = vals[3] ;
          dxy[i] = vals[4] ;
          dyy[i] = vals[5] ;
        }
      } else {
        ASSERT( false, "Expect 1, 3, or 6 output arguments" ) ;
      }
    } else {
      mexErrorMessage() ;
    }
  }
}
