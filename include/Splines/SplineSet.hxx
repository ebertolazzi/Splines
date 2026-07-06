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

#pragma once

#ifndef SPLINESET_HXX
#define SPLINESET_HXX

namespace Splines
{

  using std::abs;
  using std::sqrt;

  /*\
   |   ____        _ _            ____       _
   |  / ___| _ __ | (_)_ __   ___/ ___|  ___| |_
   |  \___ \| '_ \| | | '_ \ / _ \___ \ / _ \ __|
   |   ___) | |_) | | | | | |  __/___) |  __/ |_
   |  |____/| .__/|_|_|_| |_|\___|____/ \___|\__|
   |        |_|
  \*/

  /**
   * \class SplineSet
   * \brief A collection of splines sharing the same independent variable (x-nodes)
   *
   * The SplineSet class manages multiple 1D splines that share the same
   * independent variable values (x-nodes). This is particularly useful for
   * representing multiple dependent variables that are all functions of the
   * same independent variable.
   *
   * \par Key Features:
   * - Manages multiple splines of different types (linear, cubic, quintic, etc.)
   * - All splines share the same x-nodes
   * - Efficient evaluation of all splines at once
   * - Support for using any spline as the independent variable
   * - Automatic monotonicity checking
   * - Integration with GenericContainer for I/O
   * - Automatic differentiation support
   *
   * \par Architecture:
   * The class stores splines as unique pointers in a vector and maintains
   * additional data structures for fast lookup by name and efficient evaluation.
   *
   * \par Memory Management:
   * Uses custom allocators for efficient memory handling and cache locality.
   * All spline data is stored in contiguous memory blocks.
   *
   * \par Use Cases:
   * - Multiple physical quantities as functions of time
   * - Curve families with shared parameterization
   * - Multi-dimensional function approximation
   * - Inverse function evaluation (using spline as independent variable)
   *
   * \note All splines in the set must have the same number of points and
   *       share the same x-values.
   *
   * \see Spline, ConstantSpline, LinearSpline, CubicSpline, AkimaSpline,
   *      VanLeerSpline, PchipSpline, HermiteSpline, QuinticSpline
   */
  class SplineSet
  {
    SplineSet( SplineSet const & )                   = delete;
    SplineSet const & operator=( SplineSet const & ) = delete;

    /**
     * \class BinarySearch
     * \brief Internal helper class for efficient name-to-index lookup
     *
     * This class implements a binary search data structure for mapping
     * spline names to their indices in the spline set. It maintains
     * a sorted vector of (name, index) pairs for O(log n) lookups.
     *
     * \par Design:
     * - Uses std::vector for storage with pre-allocated capacity
     * - Maintains sorted order by spline name for binary search
     * - Supports insertion and lookup operations
     *
     * \par Complexity:
     * - Lookup: O(log n)
     * - Insertion: O(n) due to shifting elements
     *
     * \note This is an internal class and not intended for direct use.
     */
    class BinarySearch
    {
    public:
      /// Data type stored in the search structure: (name, index) pair
      using DATA_TYPE = std::pair<std::string, integer>;

    private:
      /// Sorted vector of (name, index) pairs
      mutable std::vector<DATA_TYPE> data;

    public:
      /**
       * \brief Default constructor
       *
       * Initializes an empty search structure with pre-allocated capacity.
       */
      BinarySearch()
      {
        data.clear();
        data.reserve( 256 );
      }

      /**
       * \brief Destructor
       *
       * Cleans up the internal data structure.
       */
      ~BinarySearch() { data.clear(); }

      /**
       * \brief Clear all entries from the search structure
       *
       * Removes all stored (name, index) pairs and resets capacity.
       */
      void clear()
      {
        data.clear();
        data.reserve( 256 );
      }

      /**
       * \brief Get the number of elements in the search structure
       *
       * \return Number of stored (name, index) pairs
       */
      integer n_elem() const { return integer( data.size() ); }

      /**
       * \brief Get element at specified position
       *
       * \param[in] i Index of element to retrieve
       * \return Constant reference to the (name, index) pair at position i
       *
       * \throw May throw std::out_of_range if i is out of bounds
       */
      DATA_TYPE const & get_elem( integer i ) const { return data[i]; }

      /**
       * \brief Search for a spline by name using binary search
       *
       * Finds the index of a spline given its name. The search is
       * case-sensitive and performs binary search on the sorted name list.
       *
       * \param[in] id Name of the spline to find
       * \return Index of the spline if found, -1 otherwise
       *
       * \par Algorithm:
       * 1. Perform binary search on sorted name list
       * 2. Compare strings lexicographically
       * 3. Return associated index if found, -1 otherwise
       *
       * \par Complexity:
       * Time: O(log n) where n is number of splines
       * Space: O(1)
       */
      integer search( string_view id ) const
      {
        size_t U = data.size();
        size_t L = 0;
        while ( U - L > 1 )
        {
          size_t const pos    = ( L + U ) >> 1;  // se L=U+1 --> (L+U)/2 ==> L
          string_view  id_pos = data[pos].first;
          if ( id_pos < id )
            L = pos;
          else
            U = pos;
        }
        if ( data[L].first == id ) return data[L].second;
        if ( U < data.size() )
          if ( data[U].first == id ) return data[U].second;
        return -1;  // non trovato
      }

      /**
       * \brief Insert a new (name, index) pair into the search structure
       *
       * Inserts a new entry while maintaining the sorted order by name.
       * If the name already exists, it will be updated with the new index.
       *
       * \param[in] id Name of the spline
       * \param[in] position Index associated with the name
       *
       * \par Algorithm:
       * 1. Append new element at the end
       * 2. Bubble up while previous element has greater name
       * 3. Overwrite existing entry if name already exists
       *
       * \par Complexity:
       * Time: O(n) worst-case due to element shifting
       * Space: O(1) additional
       */
      void insert( string_view const id, integer const position )
      {
        size_t pos = data.size();
        data.emplace_back( id, position );
        while ( pos > 0 )
        {
          size_t const pos1 = pos - 1;
          data[pos].first   = data[pos1].first;
          data[pos].second  = data[pos1].second;
          if ( data[pos1].first < id ) break;
          pos = pos1;
        }
        data[pos] = DATA_TYPE( id, position );
      }
    };

  protected:
    /// Name of the spline set for identification and error messages
    string const m_name;

    /// Custom memory allocator for real-type values
    Utils::Malloc<real_type> m_mem;

    /// Custom memory allocator for real-type pointers
    Utils::Malloc<real_type *> m_mem_p;

    /// Custom memory allocator for integer values
    Utils::Malloc<int> m_mem_int;

    /// Vector of unique pointers to individual splines
    vector<std::unique_ptr<Spline>> m_splines;

    /// Number of points in each spline (shared across all splines)
    integer m_npts = 0;

    /// Number of splines in the set
    integer m_nspl = 0;

    /// Array of x-nodes (independent variable values), size = m_npts
    real_type * m_X = nullptr;

    /// Array of pointers to y-values for each spline, size = m_nspl × m_npts
    real_type ** m_Y = nullptr;

    /// Array of pointers to first derivatives for each spline (if available)
    real_type ** m_Yp = nullptr;

    /// Array of pointers to second derivatives for each spline (if available)
    real_type ** m_Ypp = nullptr;

    /// Minimum y-value for each spline, size = m_nspl
    real_type * m_Ymin = nullptr;

    /// Maximum y-value for each spline, size = m_nspl
    real_type * m_Ymax = nullptr;

    /// Monotonicity information for each spline, size = m_nspl
    /// - -2: non-monotone data
    /// - -1: not monotone
    /// - 0: monotone (non-strict)
    /// - 1: strictly monotone
    int * m_is_monotone = nullptr;

    /// Map from spline name to its index for O(1) lookup
    std::map<string, integer> m_header_to_position;

  private:
    /**
     * \brief Find x-value where monotone spline intersects given y-value
     *
     * For a monotone spline, finds the x-value such that spline(x) = zeta.
     * Uses a combination of interval search and Halley's method for root finding.
     *
     * \param[in] spl Index of the spline to use as independent variable
     * \param[in] zeta Target y-value to intersect
     * \param[out] x Computed x-value such that spline[spl](x) = zeta
     * \return Pointer to the spline used as independent variable
     *
     * \throw UTILS_ASSERT if spline index is invalid
     * \throw UTILS_ASSERT if spline is not monotone
     * \throw UTILS_ASSERT if zeta is outside spline range
     *
     * \par Algorithm:
     * 1. Verify spline is monotone
     * 2. Find interval containing zeta using binary search
     * 3. Apply Halley's method (AlgoHNewton) to solve spline(x) - zeta = 0
     *
     * \note Only works with monotone splines (m_is_monotone[spl] > 0)
     *
     * \see Utils::AlgoHNewton
     */
    Spline const * intersect( integer const spl, real_type const zeta, real_type & x ) const
    {
      string msg = fmt::format( "SplineSet[{}]::intersect(...):", m_name );
      UTILS_ASSERT( spl >= 0 && spl < m_nspl, "{}\nSpline n.{} is not in SplineSet", msg, spl );
      UTILS_ASSERT(
        m_is_monotone[spl] > 0,
        "{}\nSpline n.{} is not monotone and can't be used as independent",
        msg,
        spl );
      Spline * S = m_splines[spl].get();
      // cerco intervallo intersezione
      real_type const * Y = m_Y[spl];
      UTILS_ASSERT(
        zeta >= Y[0] && zeta <= Y[m_npts - 1],
        "{} evaluation at zeta = {} is out of range: [{},{}]\n",
        msg,
        zeta,
        Y[0],
        Y[m_npts - 1] );

      integer interval = static_cast<integer>( std::lower_bound( Y, Y + m_npts, zeta ) - Y );
      if ( interval > 0 ) --interval;
      if ( Utils::is_zero( Y[interval] - Y[interval + 1] ) ) ++interval;  // degenerate interval for duplicated nodes
      if ( interval >= m_npts - 1 ) interval = m_npts - 2;

      Utils::AlgoHNewton<real_type> HN;
      class HNfun : public Utils::AlgoHNewton_base_fun<real_type>
      {
        Spline *  m_S;
        real_type m_zeta;

      public:
        explicit HNfun( Spline * S, real_type zeta ) : m_S( S ), m_zeta( zeta ) {}
        real_type eval( real_type x ) const { return m_S->eval( x ) - m_zeta; }
        real_type D( real_type x ) const { return m_S->D( x ); }
      };

      HNfun fun( S, zeta );

      // compute intersection
      x = HN.eval( m_X[interval], m_X[interval + 1], &fun );
      return S;
    }

  public:
    //!
    //! \name Constructors
    //!
    ///@{

    /**
     * \brief Construct an empty SplineSet
     *
     * Creates an uninitialized spline set with no splines. Use build()
     * or setup() methods to populate the spline set with data.
     *
     * \param[in] name Identifier for the spline set (default: "SplineSet")
     *
     * \par Memory:
     * Initializes custom allocators but allocates no memory until build() is called.
     */
    explicit SplineSet( string_view name = "SplineSet" )
      : m_name( name )
      , m_mem( fmt::format( "SplineSet[{}]::m_mem", name ) )
      , m_mem_p( fmt::format( "SplineSet[{}]::m_mem_p", name ) )
      , m_mem_int( fmt::format( "SplineSet[{}]::m_mem_int", name ) )
    {
    }

    /**
     * \brief Destructor
     *
     * Releases all allocated memory including spline objects and internal arrays.
     * Automatically handles cleanup of all resources.
     */
    virtual ~SplineSet()
    {
      m_mem.free();
      m_mem_p.free();
      m_mem_int.free();
    }
    ///@}

    //!
    //! \name Info
    //!
    ///@{

    /**
     * \brief Get the name of the spline set
     *
     * \return String view containing the spline set identifier
     */
    string_view name() const { return m_name; }

    /**
     * \brief Get the name of a specific spline in the set
     *
     * \param[in] i Index of the spline (0 ≤ i < m_nspl)
     * \return Name of the i-th spline
     *
     * \throw UTILS_ASSERT if i is out of bounds
     */
    string_view header( integer const i ) const { return m_splines[i]->name(); }

    /**
     * \brief Get all spline names as a vector of strings
     *
     * \param[out] names Vector to be filled with spline names
     *
     * \par Post-conditions:
     * - names.size() == m_nspl
     * - names[i] == header(i) for all i
     */
    void get_headers( std::vector<std::string> & names ) const
    {
      names.clear();
      names.reserve( m_nspl );
      for ( integer i = 0; i < m_nspl; ++i ) names.emplace_back( m_splines[i]->name() );
    }

    /**
     * \brief Get a formatted string listing all spline names
     *
     * \return String in format: "[ 'name1' 'name2' ... 'nameN' ]"
     *
     * \par Example:
     * \code
     * "[ 'velocity' 'acceleration' 'position' ]"
     * \endcode
     */
    string name_list() const
    {
      string tmp = "[ ";
      for ( integer i = 0; i < m_nspl; ++i ) tmp += fmt::format( "'{}' ", m_splines[i]->name() );
      tmp += "]";
      return tmp;
    }

    /**
     * \brief Get monotonicity information for a spline
     *
     * \param[in] i Index of the spline
     * \return Monotonicity code:
     *         - +1 = strictly monotone
     *         - 0 = weakly monotone (constant segments allowed)
     *         - -1 = non-monotone
     *
     * \note For cubic splines, this is determined by checking if
     *       the derivative doesn't change sign.
     */
    int is_monotone( integer const i ) const { return m_is_monotone[i]; }

    /**
     * \brief Get number of points in each spline
     *
     * \return Number of x-nodes (same for all splines in the set)
     */
    integer num_points() const { return m_npts; }

    /**
     * \brief Get number of splines in the set
     *
     * \return Number of splines managed by this SplineSet
     */
    integer num_splines() const { return m_nspl; }

    /**
     * \brief Find index of spline by name
     *
     * \param[in] hdr Name of the spline to find
     * \return Index of the spline (0 ≤ index < m_nspl)
     *
     * \throw UTILS_ASSERT if spline name not found
     */
    integer get_position( string_view hdr ) const
    {
      integer const pos = m_header_to_position.at( hdr.data() );
      UTILS_ASSERT(
        pos >= 0 && pos < m_nspl,
        "SplineSet[{}]::get_position(\"{}\") not found!\n"
        "available keys: {}\n",
        m_name,
        hdr,
        name_list() );
      return pos;
    }

    /**
     * \brief Get pointer to x-nodes array
     *
     * \return Constant pointer to array of x-values (size = m_npts)
     */
    real_type const * x_nodes() const { return m_X; }

    /**
     * \brief Get pointer to y-nodes array for a specific spline
     *
     * \param[in] i Index of the spline
     * \return Constant pointer to array of y-values for spline i (size = m_npts)
     *
     * \throw UTILS_ASSERT if i is out of bounds
     */
    real_type const * y_nodes( integer const i ) const
    {
      UTILS_ASSERT(
        i >= 0 && i < m_nspl,
        "SplineSet[{}]::y_nodes({}) argument out of range [0,{}]\n",
        m_name,
        i,
        m_nspl - 1 );
      return m_Y[i];
    }

    /**
     * \brief Get a specific x-node value
     *
     * \param[in] npt Index of the x-node (0 ≤ npt < m_npts)
     * \return X-value at the specified node
     */
    real_type x_node( integer const npt ) const { return m_X[npt]; }

    /**
     * \brief Get a specific y-node value
     *
     * \param[in] npt Index of the node point
     * \param[in] spl Index of the spline
     * \return Y-value of spline `spl` at node `npt`
     */
    real_type y_node( integer const npt, integer const spl ) const { return m_Y[spl][npt]; }

    /**
     * \brief Get minimum x-value in the domain
     *
     * \return Smallest x-node value (first element of m_X)
     */
    real_type x_min() const { return m_X[0]; }

    /**
     * \brief Get maximum x-value in the domain
     *
     * \return Largest x-node value (last element of m_X)
     */
    real_type x_max() const { return m_X[m_npts - 1]; }

    /**
     * \brief Get minimum y-value of a specific spline
     *
     * \param[in] spl Index of the spline
     * \return Minimum y-value attained by spline `spl` over its domain
     */
    real_type y_min( integer const spl ) const { return m_Ymin[spl]; }

    /**
     * \brief Get maximum y-value of a specific spline
     *
     * \param[in] spl Index of the spline
     * \return Maximum y-value attained by spline `spl` over its domain
     */
    real_type y_max( integer const spl ) const { return m_Ymax[spl]; }

    /**
     * \brief Get minimum y-value of a spline by name
     *
     * \param[in] spl Name of the spline
     * \return Minimum y-value of the specified spline
     *
     * \throw UTILS_ASSERT if spline name not found
     */
    real_type y_min( string_view spl ) const
    {
      integer idx = this->get_position( spl );
      return m_Ymin[idx];
    }

    /**
     * \brief Get maximum y-value of a spline by name
     *
     * \param[in] spl Name of the spline
     * \return Maximum y-value of the specified spline
     *
     * \throw UTILS_ASSERT if spline name not found
     */
    real_type y_max( string_view spl ) const
    {
      integer idx = this->get_position( spl );
      return m_Ymax[idx];
    }

    ///@}

    //! \name Access splines
    ///@{

    /**
     * \brief Get pointer to a spline by index
     *
     * \param[in] i Index of the spline (0 ≤ i < m_nspl)
     * \return Pointer to the i-th spline (non-const)
     *
     * \throw UTILS_ASSERT if i is out of bounds
     *
     * \note The returned pointer remains valid until the SplineSet is modified.
     */
    Spline * get_spline( integer const i ) const
    {
      UTILS_ASSERT(
        i >= 0 && i < m_nspl,
        "SplineSet[{}]::get_spline({}) argument out of range [0,{}]\n",
        m_name,
        i,
        m_nspl - 1 );
      return m_splines[i].get();
    }

    /**
     * \brief Get pointer to a spline by name
     *
     * \param[in] hdr Name of the spline
     * \return Pointer to the requested spline (non-const)
     *
     * \throw UTILS_ASSERT if spline name not found
     */
    Spline * get_spline( string_view hdr ) const
    {
      integer idx = this->get_position( hdr );
      return m_splines[idx].get();
    }

    ///@}

    //! \name Evaluate
    ///@{

#include "SplineSet_eval.hxx"
#include "SplineSet_eval_xs.hxx"

    ///@}

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    //!
    //! \name Evaluate to std vector
    //!
    ///@{

#include "SplineSet_eval_xsv.hxx"

    ///@}

    //! \name Evaluate to vector
    ///@{

#include "SplineSet_eval_xv.hxx"

    ///@}

    //! \name Evaluate to vector using another spline as independent
    ///@{

#include "SplineSet_eval2_sv.hxx"

    ///@}

    //! \name Evaluate to C vector using another spline as independent
    ///@{

#include "SplineSet_eval2_v.hxx"

    ///@}

    // \name Evaluate using another spline as independent using string
    ///@{

#include "SplineSet_eval2_ss.hxx"

    ///@}

    //! \name Evaluate using another spline as independent unsing integer
    ///@{

#include "SplineSet_eval2_ii.hxx"

    ///@}

    //! \name Evaluate onto a Generic Container
    ///@{

#include "SplineSet_eval_gc.hxx"

    ///@}

    //! \name Evaluate to GenericContainer using another spline as independent
    ///@{

#include "SplineSet_eval2_gc.hxx"

    ///@}


    ///@}

    //! \name Build Spline
    ///@{

    /**
     * \brief Create a deep copy of this spline set into another
     *
     * Copies all spline data, names, and types to the target SplineSet.
     * The target is completely overwritten and will be identical to this one.
     *
     * \param[out] S Target SplineSet to receive the copy
     *
     * \par Note:
     * Also copies flags (extrapolation, closed curve, etc.) from each spline.
     */
    void deep_copy_to( SplineSet & S ) const
    {
      std::vector<char const *> headers( m_nspl );
      std::vector<SplineType1D> stype( m_nspl );

      for ( integer i = 0; i < m_nspl; ++i )
      {
        headers[i] = m_splines[i]->name().data();
        stype[i]   = m_splines[i]->type();
      }

      S.build( m_nspl, m_npts, headers.data(), stype.data(), m_X, m_Y, m_Yp );

      // propagate flags
      for ( integer i = 0; i < m_nspl; ++i ) S.m_splines[i]->copy_flags( *m_splines[i] );
    }

#include "SplineSet_setup.hxx"
#include "SplineSet_build.hxx"

    /**
     * \brief Build spline set from GenericContainer
     *
     * \copydetails setup()
     * This is an alias for setup().
     */
    void build( GenericContainer const & gc ) { this->setup( gc ); }

    ///@}

    /**
     * \brief Get the type identifier of this spline set
     *
     * \return SplineType1D::SPLINE_SET enum value
     */
    SplineType1D type() const { return SplineType1D::SPLINE_SET; }

    /**
     * \brief Generate detailed informational string about the spline set
     *
     * \return Multi-line string containing:
     *         - Name of spline set
     *         - Number of points and splines
     *         - For each spline: index, monotonicity, and spline-specific info
     *
     * \par Example Output:
     * \code
     * SplineSet[my_set] n.points=10 n.splines=3
     * Spline n.0 is strictly monotone
     * ConstantSpline[name0] n.points = 10
     * Spline n.1 is NOT monotone
     * CubicSpline[name1] n.points = 10
     * ...
     * \endcode
     */
    string info() const
    {
      string res = fmt::format( "SplineSet[{}] n.points={} n.splines={}", name(), m_npts, m_nspl );

      for ( integer i = 0; i < m_nspl; ++i )
      {
        res += fmt::format( "\nSpline n.{} ", i );
        switch ( m_is_monotone[i] )
        {
          case -2: res += " with NON monotone data\n"; break;
          case -1: res += " is NOT monotone\n"; break;
          case 0: res += " is monotone\n"; break;
          case 1: res += " is strictly monotone\n"; break;
          default: UTILS_ERROR( "SplineSet::info classification: {} not in range {{-2,-1,0,1}}\n", m_is_monotone[i] );
        }
        res += m_splines[i]->info();
      }
      return res;
    }

    /**
     * \brief Print informational string to output stream
     *
     * \param[out] stream Output stream (e.g., std::cout, std::ofstream)
     *
     * \see info() for output format
     */
    void info( ostream_type & stream ) const { stream << this->info() << '\n'; }

    /**
     * \brief Generate a table of spline values for plotting
     *
     * Samples all splines uniformly over the x-domain and writes a
     * tab-separated table to the output stream.
     *
     * \param[out] stream Output stream for the table
     * \param[in] num_points Number of sample points (including endpoints)
     *
     * \par Output Format:
     * ```
     * s    spline0    spline1    ...    splineN
     * 0.0  val00      val01      ...    val0N
     * 0.1  val10      val11      ...    val1N
     * ...  ...        ...        ...    ...
     * 1.0  valM0      valM1      ...    valMN
     * ```
     * where M = num_points-1, N = m_nspl-1
     */
    void dump_table( ostream_type & stream, integer const num_points ) const
    {
      Malloc_real mem( "SplineSet::dump_table" );
      mem.allocate( m_nspl );
      real_type * vals = mem( m_nspl );
      stream << 's';
      for ( integer i = 0; i < m_nspl; ++i ) stream << '\t' << header( i );
      stream << '\n';

      for ( integer j = 0; j < num_points; ++j )
      {
        real_type const s = x_min() + ( ( x_max() - x_min() ) * j ) / ( num_points - 1 );
        this->eval( s, vals, 1 );
        stream << s;
        for ( integer i = 0; i < m_nspl; ++i ) stream << '\t' << vals[i];
        stream << '\n';
      }
    }

#ifdef SPLINES_BACK_COMPATIBILITY
    /**
     * \brief Backward compatibility: get monotonicity information
     * \param[in] i Spline index
     * \return Monotonicity code
     */
    int isMonotone( integer i ) const { return m_is_monotone[i]; }

    /**
     * \brief Backward compatibility: get number of points
     * \return Number of points per spline
     */
    integer numPoints() const { return m_npts; }

    /**
     * \brief Backward compatibility: get number of splines
     * \return Number of splines in set
     */
    integer numSplines() const { return m_nspl; }

    /**
     * \brief Backward compatibility: find spline index by name
     * \param[in] hdr Spline name
     * \return Spline index
     */
    integer getPosition( string_view hdr ) const { return get_position( hdr ); }

    /**
     * \brief Backward compatibility: get x-nodes array
     * \return Pointer to x-values
     */
    real_type const * xNodes() const { return m_X; }

    /**
     * \brief Backward compatibility: get y-nodes array for spline
     * \param[in] i Spline index
     * \return Pointer to y-values for spline i
     */
    real_type const * yNodes( integer i ) const;

    /**
     * \brief Backward compatibility: get specific x-node
     * \param[in] npt Node index
     * \return X-value at node npt
     */
    real_type xNode( integer npt ) const { return this->x_node( npt ); }

    /**
     * \brief Backward compatibility: get specific y-node
     * \param[in] npt Node index
     * \param[in] spl Spline index
     * \return Y-value of spline spl at node npt
     */
    real_type yNode( integer npt, integer spl ) const { return this->y_node( npt, spl ); }

    /**
     * \brief Backward compatibility: get minimum x
     * \return Minimum x-value
     */
    real_type xMin() const { return this->x_min(); }

    /**
     * \brief Backward compatibility: get maximum x
     * \return Maximum x-value
     */
    real_type xMax() const { return this->x_max(); }

    /**
     * \brief Backward compatibility: get minimum y for spline
     * \param[in] spl Spline index
     * \return Minimum y-value of spline spl
     */
    real_type yMin( integer spl ) const { return this->y_min( spl ); }

    /**
     * \brief Backward compatibility: get maximum y for spline
     * \param[in] spl Spline index
     * \return Maximum y-value of spline spl
     */
    real_type yMax( integer spl ) const { return this->y_max( spl ); }

    /**
     * \brief Backward compatibility: get minimum y by name
     * \param[in] spl Spline name
     * \return Minimum y-value of named spline
     */
    real_type yMin( string_view spl ) const { return this->y_min( spl ); }

    /**
     * \brief Backward compatibility: get maximum y by name
     * \param[in] spl Spline name
     * \return Maximum y-value of named spline
     */
    real_type yMax( string_view spl ) const { return this->y_max( spl ); }

    /**
     * \brief Backward compatibility: get spline by index
     * \param[in] i Spline index
     * \return Pointer to spline i
     */
    Spline * getSpline( integer i ) const { return this->get_spline( i ); }

    /**
     * \brief Backward compatibility: get spline by name
     * \param[in] hdr Spline name
     * \return Pointer to named spline
     */
    Spline * getSpline( string_view hdr ) const { return this->get_spline( hdr ); }
#endif
  };


}  // namespace Splines

// EOF: SplineSet.hxx
#endif
