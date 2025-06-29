/****************************************************************************\
  Copyright (c) Enrico Bertolazzi 2014
  See file license.txt
\****************************************************************************/

#include "GenericContainer/GenericContainerInterface_matlab.hh"

#include <iostream>
#include <cstdint>

#ifndef DOXYGEN_SHOULD_SKIP_THIS
using namespace std;
#endif

namespace GC_namespace {

  // ===========================================================================
  void
  mexPrint( GenericContainer const & gc ) {
    std::ostringstream ss;
    gc.print(ss);
    mexPrintf("%s\n", ss.str().data());
  }

  // ===========================================================================

  static
  void
  mx_to_vec_bool( mxArray const * mx, GenericContainer & gc ) {
    if ( mxGetNumberOfElements(mx) == 1 ) {
      gc = mxIsLogicalScalarTrue(mx);
    } else {
      mxLogical const * pr{ mxGetLogicals(mx) };
      unsigned total_num_of_elements{ static_cast<unsigned>(mxGetNumberOfElements(mx)) };
      vector_type & vec { gc.set_vector( total_num_of_elements ) };
      for ( unsigned idx{0}; idx < total_num_of_elements; ++idx )
        vec[idx] = *pr++;
    }
  }

  template <typename T>
  static
  void
  mx_to_vec_int( mxArray const * mx, GenericContainer & gc ) {
    mwSize const * dims{ mxGetDimensions(mx) };
    T * pr{ static_cast<T*>(mxGetData(mx)) };
    if ( mxGetNumberOfElements(mx) == 1 ) {
      gc.set_int( int_type(*pr) );
    } else {
      unsigned total_num_of_elements { static_cast<unsigned>(mxGetNumberOfElements(mx)) };
      unsigned number_of_dimensions  { static_cast<unsigned>(mxGetNumberOfDimensions(mx)) };
      switch ( number_of_dimensions ) {
      case 1:
        { vec_int_type & vec{ gc.set_vec_int( total_num_of_elements ) };
          for ( unsigned idx{0}; idx < total_num_of_elements; ++idx )
            vec[idx] = int_type(*pr++);
        }
        break;
      case 2:
        {
          unsigned nx{ static_cast<unsigned>(dims[0]) };
          unsigned ny{ static_cast<unsigned>(dims[1]) };
          if ( nx == 1 ) {
            vec_int_type & vec = gc.set_vec_int( total_num_of_elements );
            for ( unsigned idx{0}; idx < total_num_of_elements; ++idx )
              vec[idx] = int_type(*pr++);
          } else {
            mat_int_type & mat{ gc.set_mat_int( nx, ny ) };
            for ( unsigned j{0}; j < ny; ++j )
              for ( unsigned i{0}; i < nx; ++i )
                mat(i,j) = int_type(*pr++);
          }
        }
        break;
      default:
        mexPrintf("number_of_dimensions = %d\n", number_of_dimensions );
        break;
      }
    }
  }

  template <typename T>
  static
  void
  mx_to_vec_long( mxArray const * mx, GenericContainer & gc ) {
    mwSize const * dims{ mxGetDimensions(mx) };
    T * pr{ static_cast<T*>(mxGetData(mx)) };
    if ( mxGetNumberOfElements(mx) == 1 ) {
      gc.set_long( long_type(*pr) );
    } else {
      unsigned number_of_dimensions  { static_cast<unsigned>(mxGetNumberOfDimensions(mx)) };
      unsigned total_num_of_elements { static_cast<unsigned>(mxGetNumberOfElements(mx)) };
      switch ( number_of_dimensions ) {
      case 1:
        { vec_long_type & vec{ gc.set_vec_long( total_num_of_elements ) };
          for ( unsigned idx{0}; idx < total_num_of_elements; ++idx )
            vec[idx] = long_type(*pr++);
        }
        break;
      case 2:
        {
          unsigned nx{ static_cast<unsigned>(dims[0]) };
          unsigned ny{ static_cast<unsigned>(dims[1]) };
          if ( nx == 1 ) {
            vec_long_type & vec{ gc.set_vec_long( total_num_of_elements ) };
            for ( unsigned idx{0}; idx < total_num_of_elements; ++idx )
              vec[idx] = long_type(*pr++);
          } else {
            mat_long_type & mat{ gc.set_mat_long( nx, ny ) };
            for ( unsigned j{0}; j < ny; ++j )
              for ( unsigned i{0}; i < nx; ++i )
                mat(i,j) = long_type(*pr++);
          }
        }
        break;
      default:
        mexPrintf("number_of_dimensions = %d\n", number_of_dimensions );
        break;
      }
    }
  }

  // ===========================================================================

  template <typename T>
  static
  void
  mx_to_vec_real( mxArray const * mx, GenericContainer & gc ) {
    mwSize const * dims = mxGetDimensions(mx);
    T * pr = static_cast<T*>(mxGetData(mx));
    if ( mxGetNumberOfElements(mx) == 1 ) {
      gc = real_type(*pr);
    } else {
      unsigned total_num_of_elements{ static_cast<unsigned>(mxGetNumberOfElements(mx)) };
      unsigned number_of_dimensions  { static_cast<unsigned>(mxGetNumberOfDimensions(mx)) };
      switch ( number_of_dimensions ) {
      case 1:
        { vec_real_type & vec{ gc.set_vec_real( total_num_of_elements ) };
          for ( mwSize idx = 0; idx < total_num_of_elements; ++idx )
            vec[idx] = real_type(*pr++);
        }
        break;
      case 2:
        {
          unsigned nx{ static_cast<unsigned>(dims[0]) };
          unsigned ny{ static_cast<unsigned>(dims[1]) };
          if ( nx == 1 ) {
            vec_real_type & vec{ gc.set_vec_real( total_num_of_elements ) };
            for ( unsigned idx{0}; idx < total_num_of_elements; ++idx )
              vec[idx] = real_type(*pr++);
          } else {
            mat_real_type & mat = gc.set_mat_real( nx, ny );
            for ( unsigned idx{0}; idx < total_num_of_elements; ++idx )
              mat[idx] = real_type(*pr++);
          }
        }
        break;
      default:
        mexPrintf("number_of_dimensions = %d\n", number_of_dimensions );
        break;
      }
    }
  }
  
  #if MX_HAS_INTERLEAVED_COMPLEX

  template <typename T>
  static
  void
  mx_to_vec_complex(mxArray const* mx, GenericContainer& gc) {
    mwSize const* dims = mxGetDimensions(mx);

    // Dati complessi interleaved
    mxComplexDouble const* data = reinterpret_cast<mxComplexDouble const*>(mxGetData(mx));

    if (mxGetNumberOfElements(mx) == 1) {
      gc = complex_type(real_type(data[0].real), real_type(data[0].imag));
    } else {
      unsigned number_of_dimensions  { static_cast<unsigned>(mxGetNumberOfDimensions(mx)) };
      unsigned total_num_of_elements { static_cast<unsigned>(mxGetNumberOfElements(mx)) };

      switch (number_of_dimensions) {
        case 1: {
          vec_complex_type& vec = gc.set_vec_complex(total_num_of_elements);
          for ( unsigned idx{0}; idx < total_num_of_elements; ++idx) {
            vec[idx] = complex_type(real_type(data[idx].real), real_type(data[idx].imag));
          }
          break;
        }
        case 2:
          {
            unsigned nx{ static_cast<unsigned>(dims[0]) };
            unsigned ny{ static_cast<unsigned>(dims[1]) };
            if ( nx == 1) {
              vec_complex_type& vec = gc.set_vec_complex(total_num_of_elements);
              for ( unsigned idx{0}; idx < total_num_of_elements; ++idx) {
                vec[idx] = complex_type(real_type(data[idx].real), real_type(data[idx].imag));
              }
            } else {
              mat_complex_type & mat{ gc.set_mat_complex( nx, ny ) };
              for ( unsigned idx{0}; idx < total_num_of_elements; ++idx) {
                mat[idx] = complex_type(real_type(data[idx].real), real_type(data[idx].imag));
              }
            }
          }
          break;
        default:
          mexPrintf("number_of_dimensions = %d\n", number_of_dimensions);
          break;
      }
    }
  }

  #else

  // ===========================================================================

  template <typename T>
  static
  void
  mx_to_vec_complex( mxArray const * mx, GenericContainer & gc ) {
    mwSize const * dims{ mxGetDimensions(mx) };
    T * pr{ static_cast<T*>(mxGetData(mx)) };
    T * pi{ static_cast<T*>(mxGetImagData(mx)) };

    if ( mxGetNumberOfElements(mx) == 1 ) {
      gc = complex_type(real_type(*pr),real_type(*pi));
    } else {
      unsigned number_of_dimensions  { static_cast<unsigned>(mxGetNumberOfDimensions(mx)) };
      unsigned total_num_of_elements { static_cast<unsigned>(mxGetNumberOfElements(mx)) };
      switch ( number_of_dimensions ) {
      case 1:
        { vec_complex_type & vec{ gc.set_vec_complex( total_num_of_elements ) };
          for ( unsigned idx{0}; idx < total_num_of_elements; ++idx )
            vec[idx] = complex_type(real_type(*pr++),real_type(*pi++));
        }
        break;
      case 2:
        if ( nx == 1 ) {
          vec_complex_type & vec{ gc.set_vec_complex( total_num_of_elements ) };
          for ( unsigned idx{0}; idx < total_num_of_elements; ++idx )
            vec[idx] = complex_type(real_type(*pr++),real_type(*pi++));
        } else {
          mat_complex_type & mat{ gc.set_mat_complex( nx, ny ) };
          for ( unsigned idx{0}; idx < total_num_of_elements; ++idx )
            mat[idx] = complex_type(real_type(*pr++),real_type(*pi++));
        }
        break;
      default:
        mexPrintf("number_of_dimensions = %d\n", number_of_dimensions );
        break;
      }
    }
  }

  #endif

  // ===========================================================================

  static
  void
  mx_to_vector( mxArray const * mx, GenericContainer & gc ) {
    unsigned ne{ static_cast<unsigned>( mxGetNumberOfElements(mx) ) };
    gc.set_vector(ne);
    for ( unsigned idx{0}; idx < ne; ++idx ) {
      GenericContainer & gc1 = gc[idx];
      mxArray const * cell = mxGetCell(mx, idx);
      mxArray_to_GenericContainer( cell, gc1 );
    }
  }

  // ===========================================================================

  static
  void
  mx_to_map( mxArray const * mx, GenericContainer & gc ) {
    gc.set_map();
    int numFields{ static_cast<int>( mxGetNumberOfFields(mx) ) };
    for ( int ifield{0}; ifield < numFields; ++ifield ) {
      char const * field_name{ static_cast<char const *>( mxGetFieldNameByNumber(mx,ifield) ) };
      GenericContainer & gc1 = gc[field_name];
      unsigned nr{ static_cast<unsigned>( mxGetM(mx) ) };
      unsigned nc{ static_cast<unsigned>( mxGetN(mx) ) };
      if ( nc == 1 && nr == 1 ) {
        mxArray const * mxField = mxGetFieldByNumber(mx,0,ifield);
        mxArray_to_GenericContainer( mxField, gc1 );
      } else {
        gc1.set_vector(nr*nc);
        for ( unsigned i{0}; i < nr*nc; ++i ) {
          mxArray const * mxField{ mxGetFieldByNumber(mx,i,ifield) };
          mxArray_to_GenericContainer( mxField, gc1[i] );
        }
      }
    }
  }

  // ===========================================================================

  void
  mxSparse_to_GenericContainer( mxArray const * mx, GenericContainer & gc ) {

    mwIndex const * irs = mxGetIr(mx);
    mwIndex const * jcs = mxGetJc(mx);
    size_t          nc  = mxGetN(mx);
    mwIndex         nnz = jcs[nc];

    vec_int_type & jc = gc["jc"].set_vec_int( nc+1 );
    vec_int_type & ir = gc["ir"].set_vec_int( nnz );

    for ( unsigned i{0}; i <= nc; ++i ) jc[i] = jcs[i];
    for ( unsigned i{0}; i < nnz; ++i ) ir[i] = irs[i];

    real_type * sr = mxGetPr(mx);
    if ( mxIsComplex( mx ) ) {
      #if MX_HAS_INTERLEAVED_COMPLEX
        mxComplexDouble* data = reinterpret_cast<mxComplexDouble*>(mxGetData(mx));
        vec_complex_type& val = gc["values"].set_vec_complex(nnz);
        for ( unsigned i{0}; i < nnz; ++i ) {
          val[i] = complex_type(data[i].real, data[i].imag);
        }
      #else
        real_type* si = mxGetPi(mx);
        vec_complex_type& val = gc["values"].set_vec_complex(nnz);
        for ( unsigned i{0}; i < nnz; ++i ) {
          val[i] = complex_type(sr[i], si[i]);
        }
      #endif
    } else {
      vec_real_type & val = gc["values"].set_vec_real( nnz );
      for ( unsigned i{0}; i < nnz; ++i ) val[i] = sr[i];
    }
  }

  // ===========================================================================

  void
  mxArray_to_GenericContainer( mxArray const * mx, GenericContainer & gc ) {
    gc.clear();
    if ( mx == nullptr ) return;
    mxClassID category = mxGetClassID(mx);
    //mexPrintf("\n\n\n%s\n\n\n",mxGetClassName(mx));
    if ( category == mxCELL_CLASS ) {
      mx_to_vector(mx,gc);
    } else if ( category == mxSTRUCT_CLASS ) {
      mx_to_map(mx,gc);
    } else if ( category == mxCHAR_CLASS ) {
      gc = mxArrayToString(mx);
    } else if ( mxIsSparse(mx) ) {
      mxSparse_to_GenericContainer( mx, gc );
    } else if ( mxIsClass(mx,"string") ) {
      // https://it.mathworks.com/matlabcentral/answers/330929-how-to-access-matlab-string-data-in-mex-c
      // Matlab's String class is encapsulated,
      // use Matlab call to convert it to char array
      mxArray * string_class[1];
      mxArray * char_array[1];
      string_class[0] = const_cast<mxArray*>(mx);
      mexCallMATLAB(1, char_array, 1, string_class, "char");
      // Parse the char array to create an std::string
      int buflen = mxGetN(char_array[0])*sizeof(mxChar)+1;
      char* buf = new char[buflen];
      mxGetString(char_array[0],buf,buflen);
      gc = buf;
      delete [] buf;
    } else {
      if ( mxIsComplex(mx) ) {
        switch (category)  {
          case mxLOGICAL_CLASS: mx_to_vec_bool(mx,gc);              break;
          case mxINT8_CLASS:    mx_to_vec_complex<int8_t>(mx,gc);   break;
          case mxUINT8_CLASS:   mx_to_vec_complex<uint8_t>(mx,gc);  break;
          case mxINT16_CLASS:   mx_to_vec_complex<int16_t>(mx,gc);  break;
          case mxUINT16_CLASS:  mx_to_vec_complex<uint16_t>(mx,gc); break;
          case mxINT32_CLASS:   mx_to_vec_complex<int32_t>(mx,gc);  break;
          case mxUINT32_CLASS:  mx_to_vec_complex<uint32_t>(mx,gc); break;
          case mxINT64_CLASS:   mx_to_vec_complex<int64_t>(mx,gc);  break;
          case mxUINT64_CLASS:  mx_to_vec_complex<uint64_t>(mx,gc); break;
          case mxSINGLE_CLASS:  mx_to_vec_complex<float>(mx,gc);    break;
          case mxDOUBLE_CLASS:  mx_to_vec_complex<double>(mx,gc);   break;
          default:
            mexPrintf("Complex Class ID = %d not converted!\n", category );
          break;
        }
      } else {
        switch (category)  {
          case mxLOGICAL_CLASS: mx_to_vec_bool(mx,gc);           break;
          case mxINT8_CLASS:    mx_to_vec_int<int8_t>(mx,gc);    break;
          case mxUINT8_CLASS:   mx_to_vec_int<uint8_t>(mx,gc);   break;
          case mxINT16_CLASS:   mx_to_vec_int<int16_t>(mx,gc);   break;
          case mxUINT16_CLASS:  mx_to_vec_int<uint16_t>(mx,gc);  break;
          case mxINT32_CLASS:   mx_to_vec_int<int32_t>(mx,gc);   break;
          case mxUINT32_CLASS:  mx_to_vec_int<uint32_t>(mx,gc);  break;
          case mxINT64_CLASS:   mx_to_vec_long<int64_t>(mx,gc);  break;
          case mxUINT64_CLASS:  mx_to_vec_long<uint64_t>(mx,gc); break;
          case mxSINGLE_CLASS:  mx_to_vec_real<float>(mx,gc);    break;
          case mxDOUBLE_CLASS:  mx_to_vec_real<double>(mx,gc);   break;
          default:
            mexPrintf("Class ID = %d not converted!\n", category );
          break;
        }
      }
    }
  }

  // ===========================================================================

  void
  to_mxArray( bool const & val_in, mxArray * & mx ) {
    mxLogical val = val_in ? 1 : 0;
    mx = mxCreateLogicalScalar(val);
  }

  void
  to_mxArray( int_type const & val, mxArray * & mx ) {
    mwSize dims[2] = {1,1};
    mx = mxCreateNumericArray(2,dims,mxINT64_CLASS,mxREAL);
    *static_cast<mwSize *>(mxGetData(mx)) = mwSize(val);
  }

  void
  to_mxArray( long_type const & val, mxArray * & mx ) {
    mwSize dims[2] = {1,1};
    mx = mxCreateNumericArray(2,dims,mxINT64_CLASS,mxREAL);
    *static_cast<int64_t*>(mxGetData(mx)) = int64_t(val);
  }

  void
  to_mxArray( real_type const & val, mxArray * & mx ) {
    mx = mxCreateDoubleScalar(val);
  }

  void
  to_mxArray(complex_type const& val, mxArray*& mx) {
    mwSize dims[2] = {1, 1};
    mx = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxCOMPLEX);
  
  #if MX_HAS_INTERLEAVED_COMPLEX
    mxComplexDouble* data = reinterpret_cast<mxComplexDouble*>(mxGetData(mx));
    data[0].real = val.real();
    data[0].imag = val.imag();
  #else
    *mxGetPr(mx) = val.real();
    *mxGetPi(mx) = val.imag();
  #endif
  }

  void
  to_mxArray( string_view val, mxArray * & mx ) {
    mx = mxCreateString( val.data() );
  }

  void
  to_mxArray( vec_bool_type const & val, mxArray * & mx ) {
    mwSize dims[2] = { mwSize(val.size()), 1 };
    mx = mxCreateNumericArray(2,dims,mxLOGICAL_CLASS,mxREAL);
    mxLogical * ptr = static_cast<mxLogical*>(mxGetData(mx));
    for ( mwSize i = 0; i < dims[0]; ++i ) ptr[i] = val[i];
  }

  void
  to_mxArray( vec_int_type const & val, mxArray * & mx ) {
  cout << "in  to_mxArray( vec_int_type\n";
    mwSize dims[2] = { mwSize(val.size()), 1 };
    mx = mxCreateNumericArray(2,dims,mxINT32_CLASS,mxREAL);
    int32_t * ptr = static_cast<int32_t*>(mxGetData(mx));
    for ( mwSize i = 0; i < dims[0]; ++i ) ptr[i] = int32_t(val[i]);
  }

  void
  to_mxArray( vec_long_type const & val, mxArray * & mx ) {
    mwSize dims[2] = { mwSize(val.size()), 1 };
    mx = mxCreateNumericArray(2,dims,mxINT64_CLASS,mxREAL);
    int64_t * ptr = static_cast<int64_t*>(mxGetData(mx));
    for ( mwSize i = 0; i < dims[0]; ++i ) ptr[i] = int64_t(val[i]);
  }

  void
  to_mxArray( vec_real_type const & val, mxArray * & mx ) {
    mwSize dims[2] = { mwSize(val.size()), 1 };
    mx = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxREAL);
    real_type * ptr = mxGetPr(mx);
    for ( mwSize i = 0; i < dims[0]; ++i ) ptr[i] = val[i];
  }


  void
  to_mxArray(vec_complex_type const& val, mxArray*& mx) {
    mwSize dims[2] = { mwSize(val.size()), 1 };
    mx = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxCOMPLEX);

  #if MX_HAS_INTERLEAVED_COMPLEX
    mxComplexDouble* data = reinterpret_cast<mxComplexDouble*>(mxGetData(mx));
    for (mwSize i = 0; i < dims[0]; ++i) {
      data[i].real = val[i].real();
      data[i].imag = val[i].imag();
    }
  #else
    real_type* ptr = mxGetPr(mx);
    real_type* pti = mxGetPi(mx);
    for (mwSize i = 0; i < dims[0]; ++i) {
      ptr[i] = val[i].real();
      pti[i] = val[i].imag();
    }
  #endif
  }

  void
  to_mxArray( vec_string_type const & val, mxArray * & mx ) {
    mwSize dims[2] = { mwSize(val.size()), 1 };
    mx = mxCreateCellMatrix(dims[0], dims[1]);
    for( mwSize i = 0; i < dims[0]; ++i )
      mxSetCell(mx,i,mxCreateString( val[i].data()) );
  }

  void
  to_mxArray( mat_int_type const & val, mxArray * & mx ) {
    mwSize dims[2] = { mwSize(val.num_rows()), mwSize(val.num_cols()) };
    mx = mxCreateNumericArray(2,dims,mxINT32_CLASS,mxREAL);
    int32_t * ptr = static_cast<int32_t*>(mxGetData(mx));
    mwSize k = 0;
    for ( mwSize j{0}; j < dims[1]; ++j )
      for ( mwSize i{0}; i < dims[0]; ++i )
        ptr[k++] = val(i,j);
  }

  void
  to_mxArray( mat_long_type const & val, mxArray * & mx ) {
    mwSize dims[2] = { mwSize(val.num_rows()), mwSize(val.num_cols()) };
    mx = mxCreateNumericArray(2,dims,mxINT64_CLASS,mxREAL);
    int64_t * ptr = static_cast<int64_t*>(mxGetData(mx));
    mwSize k = 0;
    for ( mwSize j{0}; j < dims[1]; ++j )
      for ( mwSize i{0}; i < dims[0]; ++i )
        ptr[k++] = val(i,j);
  }

  void
  to_mxArray( mat_real_type const & val, mxArray * & mx ) {
    mwSize dims[2] = { mwSize(val.num_rows()), mwSize(val.num_cols()) };
    mx = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxREAL);
    real_type * ptr = mxGetPr(mx);
    mwSize k = 0;
    for ( mwSize j{0}; j < dims[1]; ++j )
      for ( mwSize i{0}; i < dims[0]; ++i )
        ptr[k++] = val(i,j);
  }

  void
  to_mxArray(mat_complex_type const& val, mxArray*& mx) {
    mwSize dims[2] = { mwSize(val.num_rows()), mwSize(val.num_cols()) };
    mx = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxCOMPLEX);
  
  #if MX_HAS_INTERLEAVED_COMPLEX
    mxComplexDouble* data = reinterpret_cast<mxComplexDouble*>(mxGetData(mx));
    mwSize k = 0;
    for ( mwSize j{0}; j < dims[1]; ++j ) {
      for ( mwSize i{0}; i < dims[0]; ++i ) {
        complex_type const & vij{ val(i, j) };
        data[k].real = vij.real();
        data[k].imag = vij.imag();
        ++k;
      }
    }
  #else
    real_type* ptr = mxGetPr(mx);
    real_type* pti = mxGetPi(mx);
    mwSize k = 0;
    for ( mwSize j{0}; j < dims[1]; ++j ) {
      for ( mwSize i{0}; i < dims[0]; ++i ) {
        complex_type const & vij{ val(i, j) };
        ptr[k] = vij.real();
        pti[k] = vij.imag();
        ++k;
      }
    }
  #endif
  }

  // ===========================================================================

  void
  GenericContainer_to_mxArray( GenericContainer const & gc, mxArray * & mx ) {
    static string_view where{"in GenericContainer_to_mxArray: "};
    mwSize dims[2]{1,1};
    switch ( gc.get_type() ) {
    case GC_type::NOTYPE:
    case GC_type::POINTER:
    case GC_type::VEC_POINTER:
      mx = mxCreateDoubleMatrix(0,0,mxREAL);
      break;
    case GC_type::BOOL:
      {
        mxLogical val = gc.get_bool() ? 1 : 0;
        mx = mxCreateLogicalScalar(val);
      }
      break;
    case GC_type::INTEGER:
      mx = mxCreateNumericArray(2,dims,mxINT64_CLASS,mxREAL);
      *static_cast<mwSize *>(mxGetData(mx)) = gc.get_int();
      break;
    case GC_type::LONG:
      mx = mxCreateNumericArray(2,dims,mxINT64_CLASS,mxREAL);
      *static_cast<mwSize *>(mxGetData(mx)) = gc.get_long();
      break;
    case GC_type::REAL:
      mx = mxCreateDoubleScalar(gc.get_real());
      break;
    case GC_type::COMPLEX:
      {
        mx = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxCOMPLEX);
        real_type re, im;
        gc.get_complex_number(re, im);
      
      #if MX_HAS_INTERLEAVED_COMPLEX
        mxComplexDouble* data = reinterpret_cast<mxComplexDouble*>(mxGetData(mx));
        data[0].real = re;
        data[0].imag = im;
      #else
        *mxGetPr(mx) = re;
        *mxGetPi(mx) = im;
      #endif
      }
      break;
    case GC_type::STRING:
      mx = mxCreateString( gc.get_string().data() );
      break;
    case GC_type::VEC_BOOL:
      {
        dims[1] = gc.get_num_elements();
        mx = mxCreateNumericArray(2,dims,mxLOGICAL_CLASS,mxREAL);
        mxLogical * ptr = static_cast<mxLogical*>(mxGetData(mx));
        for ( mwSize i = 0; i < dims[1]; ++i ) ptr[i] = gc.get_bool_at(i,where);
      }
      break;
    case GC_type::VEC_INTEGER:
      {
        dims[1] = gc.get_num_elements();
        mx = mxCreateNumericArray(2,dims,mxINT64_CLASS,mxREAL);
        mwSize * ptr = static_cast<mwSize*>(mxGetData(mx));
        for ( mwSize i = 0; i < dims[1]; ++i ) ptr[i] = gc.get_int_at(i,where);
      }
      break;
    case GC_type::VEC_LONG:
      {
        dims[1] = gc.get_num_elements();
        mx = mxCreateNumericArray(2,dims,mxINT64_CLASS,mxREAL);
        mwSize * ptr = static_cast<mwSize*>(mxGetData(mx));
        for ( mwSize i = 0; i < dims[1]; ++i ) ptr[i] = gc.get_long_at(i,where);
      }
      break;
    case GC_type::VEC_REAL:
      {
        dims[1] = gc.get_num_elements();
        mx = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxREAL);
        real_type * ptr = mxGetPr(mx);
        for ( mwSize i = 0; i < dims[1]; ++i ) ptr[i] = gc.get_real_at(i,where);
      }
      break;
    case GC_type::VEC_COMPLEX:
      {
        dims[1] = gc.get_num_elements();
        mx = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxCOMPLEX);
      
      #if MX_HAS_INTERLEAVED_COMPLEX
        mxComplexDouble* data = reinterpret_cast<mxComplexDouble*>(mxGetData(mx));
        for (mwSize i{0}; i < dims[1]; ++i) {
          real_type re, im;
          gc.get_complex_number_at(i, re, im);
          data[i].real = re;
          data[i].imag = im;
        }
      #else
        real_type* ptr = mxGetPr(mx);
        real_type* pti = mxGetPi(mx);
        for (mwSize i = 0; i < dims[1]; ++i) {
          gc.get_complex_number_at(i, ptr[i], pti[i]);
        }
      #endif
      }
      break;
    case GC_type::VEC_STRING:
      {
        dims[1] = gc.get_num_elements();
        mx = mxCreateCellMatrix(dims[0], dims[1]);
        for( mwSize i{0}; i < dims[1]; ++i )
          mxSetCell(mx,i,mxCreateString( gc.get_string_at(i,where).data()));
      }
      break;
    case GC_type::MAT_INTEGER:
      {
        dims[0] = gc.num_rows();
        dims[1] = gc.num_cols();
        mx = mxCreateNumericArray(2,dims,mxINT32_CLASS,mxREAL);
        int_type * ptr = static_cast<int_type *>(mxGetData(mx));
        mwSize k = 0;
        for ( mwSize j{0}; j < dims[1]; ++j )
          for ( mwSize i{0}; i < dims[0]; ++i )
            ptr[k++] = gc.get_int_at(i,j,where);
      }
      break;
    case GC_type::MAT_LONG:
      {
        dims[0] = gc.num_rows();
        dims[1] = gc.num_cols();
        mx = mxCreateNumericArray(2,dims,mxINT64_CLASS,mxREAL);
        long_type * ptr = static_cast<long_type *>(mxGetData(mx));
        mwSize k = 0;
        for ( mwSize j{0}; j < dims[1]; ++j )
          for ( mwSize i{0}; i < dims[0]; ++i )
            ptr[k++] = gc.get_long_at(i,j,where);
        }
      break;
    case GC_type::MAT_REAL:
      {
        dims[0] = gc.num_rows();
        dims[1] = gc.num_cols();
        mx = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxREAL);
        real_type * ptr = mxGetPr(mx);
        mwSize k = 0;
        for ( mwSize j{0}; j < dims[1]; ++j )
          for ( mwSize i{0}; i < dims[0]; ++i )
            ptr[k++] = gc.get_real_at(i,j,where);
      }
      break;
    case GC_type::MAT_COMPLEX:
      {
        dims[0] = gc.num_rows();
        dims[1] = gc.num_cols();
        mx = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxCOMPLEX);
      
      #if MX_HAS_INTERLEAVED_COMPLEX
        mxComplexDouble* data = reinterpret_cast<mxComplexDouble*>(mxGetData(mx));
        mwSize k = 0;
        for ( mwSize j{0}; j < dims[1]; ++j ) {
          for ( mwSize i{0}; i < dims[0]; ++i ) {
            complex_type val = gc.get_complex_at(i, j, where);
            data[k].real = val.real();
            data[k].imag = val.imag();
            ++k;
          }
        }
      #else
        real_type* ptr = mxGetPr(mx);
        real_type* pti = mxGetPi(mx);
        mwSize k = 0;
        for ( mwSize j{0}; j < dims[1]; ++j ) {
          for ( mwSize i{0}; i < dims[0]; ++i ) {
            complex_type val = gc.get_complex_at(i, j, where);
            ptr[k] = val.real();
            pti[k] = val.imag();
            ++k;
          }
        }
      #endif
      }
      break;
    case GC_type::VECTOR:
      {
        dims[1] = gc.get_num_elements();
        mx = mxCreateCellMatrix(dims[0], dims[1]);
        for( mwSize i = 0; i < dims[1]; ++i ) {
          mxArray * mxi = nullptr;
          GenericContainer_to_mxArray( gc[i], mxi );
          if ( mxi != nullptr ) mxSetCell( mx, i, mxi );
        }
      }
      break;
    case GC_type::MAP:
      {
        map_type const & mappa = gc.get_map();
        vector<char const *> fieldnames;
        int nfield = mappa.size();
        fieldnames.reserve(nfield);
        for ( map_type::const_iterator im = mappa.begin(); im != mappa.end(); ++im )
          fieldnames.push_back(im->first.data());

        mx = mxCreateStructMatrix(1,1,nfield,&fieldnames.front());

        int ifield = 0;
        for ( map_type::const_iterator im = mappa.begin(); im != mappa.end(); ++im, ++ifield ) {
          mxArray * mxi = nullptr;
          GenericContainer_to_mxArray( im->second, mxi );
          if ( mxi != nullptr ) mxSetFieldByNumber( mx, 0, ifield, mxi );
        }
      }
      break;
    }
  }

}
