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
 |   ____        _ _            ____       _
 |  / ___| _ __ | (_)_ __   ___/ ___|  ___| |_
 |  \___ \| '_ \| | | '_ \ / _ \___ \ / _ \ __|
 |   ___) | |_) | | | | | |  __/___) |  __/ |_
 |  |____/| .__/|_|_|_| |_|\___|____/ \___|\__|
 |        |_|
\*/

namespace Splines {

  //! Splines Management Class
  class SplineSet {

    SplineSet( SplineSet const & ) = delete;
    SplineSet const & operator = ( SplineSet const & ) = delete;

    class BinarySearch {
    public:
      typedef std::pair<std::string,integer> DATA_TYPE;
    private:
      mutable std::vector<DATA_TYPE> data;
    public:
      BinarySearch() { data.clear(); data.reserve(256); }
      ~BinarySearch() { data.clear(); }

      void clear() { data.clear(); data.reserve(256); }

      integer n_elem() const { return integer(data.size()); }

      DATA_TYPE const &
      get_elem( integer i ) const { return data[size_t(i)]; }

      integer search( std::string const & id ) const;
      void    insert( std::string const & id, integer position );
    };

  protected:

    string const m_name;

    Utils::Malloc<real_type>  m_baseValue;
    Utils::Malloc<real_type*> m_basePointer;
    Utils::Malloc<void*>      m_baseSplines;
    Utils::Malloc<int>        m_baseInt;

    integer m_npts;
    integer m_nspl;

    real_type *  m_X;
    real_type ** m_Y;
    real_type ** m_Yp;
    real_type ** m_Ypp;
    real_type *  m_Ymin;
    real_type *  m_Ymax;

    Spline ** m_splines;
    int     * m_is_monotone;

    BinarySearch m_header_to_position;

  private:

    //! 
    //! find `x` value such that the monotone spline
    //! `(spline[spl])(x)` intersect the value `zeta`
    //! 
    Spline const *
    intersect( integer spl, real_type zeta, real_type & x ) const;

  public:

    //! \name Constructors
    ///@{

    //! spline constructor
    SplineSet( string const & name = "SplineSet" );

    //! spline destructor
    virtual
    ~SplineSet();
    ///@}

    //! \name Info
    ///@{
    string const & name() const { return m_name; }

    string const &
    header( integer i ) const
    { return m_splines[i]->name(); }

    // vectorial values
    //! fill a vector of strings with the names of the splines
    void get_headers( std::vector<std::string> & names ) const;

    string name_list() const;

    // return +1 = strict monotone, 0 weak monotone, -1 non monotone
    int
    isMonotone( integer i ) const
    { return m_is_monotone[i]; }

    //! return the number of support points of the splines
    integer numPoints() const { return m_npts; }

    //! return the number splines in the spline set
    integer numSplines() const { return m_nspl; }

    //! return the column with header(i) == hdr, return -1 if not found
    integer getPosition( char const * hdr ) const;

    //! return the vector of values of x-nodes
    real_type const * xNodes() const { return m_X; }

    //! return the vector of values of x-nodes
    real_type const * yNodes( integer i ) const;

    //! return the i-th node of the spline (x component).
    real_type
    xNode( integer npt ) const
    { return m_X[npt]; }

    //! return the i-th node of the spline (y component).
    real_type
    yNode( integer npt, integer spl ) const
    { return m_Y[spl][npt]; }

    //! return x-minumum spline value
    real_type xMin() const { return m_X[0]; }

    //! return x-maximum spline value
    real_type xMax() const { return m_X[m_npts-1]; }

    //! return y-minumum spline value
    real_type yMin( integer spl ) const { return m_Ymin[size_t(spl)]; }

    //! return y-maximum spline value
    real_type yMax( integer spl ) const { return m_Ymax[size_t(spl)]; }

    //! return y-minumum spline value
    real_type
    yMin( char const * spl ) const {
      integer idx = this->getPosition(spl);
      return m_Ymin[idx];
    }

    //! return y-maximum spline value
    real_type
    yMax( char const * spl ) const {
      integer idx = this->getPosition(spl);
      return m_Ymax[idx];
    }

    ///@}

    //! \name Access splines
    ///@{

    //! Return pointer to the `i`-th spline
    Spline * getSpline( integer i ) const;

    //! Return pointer to the `i`-th spline
    Spline *
    getSpline( char const * hdr ) const {
      integer idx = this->getPosition(hdr);
      return m_splines[idx];
    }

    ///@}

    //! \name Evaluate
    ///@{

    //! Evaluate spline value
    real_type
    operator () ( real_type x, integer spl ) const {
      Spline const * S = this->getSpline(spl);
      return (*S)(x);
    }

    real_type
    eval( real_type x, integer spl ) const {
      Spline const * S = this->getSpline(spl);
      return (*S)(x);
    }

    //! First derivative
    real_type
    D( real_type x, integer spl ) const {
      Spline const * S = this->getSpline(spl);
      return S->D(x);
    }

    real_type
    eval_D( real_type x, integer spl ) const {
      Spline const * S = this->getSpline(spl);
      return S->D(x);
    }

    //! Second derivative
    real_type
    DD( real_type x, integer spl ) const {
      Spline const * S = this->getSpline(spl);
      return S->DD(x);
    }

    real_type
    eval_DD( real_type x, integer spl ) const {
      Spline const * S = this->getSpline(spl);
      return S->DD(x);
    }

    //! Third derivative
    real_type
    DDD( real_type x, integer spl ) const {
      Spline const * S = this->getSpline(spl);
      return S->DDD(x);
    }

    real_type
    eval_DDD( real_type x, integer spl ) const {
      Spline const * S = this->getSpline(spl);
      return S->DDD(x);
    }

    //! 4th derivative
    real_type
    DDDD( real_type x, integer spl ) const {
      Spline const * S = this->getSpline(spl);
      return S->DDDD(x);
    }

    real_type
    eval_DDDD( real_type x, integer spl ) const {
      Spline const * S = this->getSpline(spl);
      return S->DDDD(x);
    }

    //! 5th derivative
    real_type
    DDDDD( real_type x, integer spl ) const {
      Spline const * S = this->getSpline(spl);
      return S->DDDDD(x);
    }

    real_type
    eval_DDDDD( real_type x, integer spl ) const {
      Spline const * S = this->getSpline(spl);
      return S->DDDDD(x);
    }

    //! Evaluate spline value
    real_type
    eval( real_type x, char const * name ) const {
      Spline const * S = this->getSpline(name);
      return (*S)(x);
    }

    //! First derivative
    real_type
    eval_D( real_type x, char const * name ) const {
      Spline const * S = this->getSpline(name);
      return S->D(x);
    }

    //! Second derivative
    real_type
    eval_DD( real_type x, char const * name ) const {
      Spline const * S = this->getSpline(name);
      return S->DD(x);
    }

    //! Third derivative
    real_type
    eval_DDD( real_type x, char const * name ) const {
      Spline const * S = this->getSpline(name);
      return S->DDD(x);
    }

    //! 4th derivative
    real_type
    eval_DDDD( real_type x, char const * name ) const {
      Spline const * S = this->getSpline(name);
      return S->DDDD(x);
    }

    //! 5th derivative
    real_type
    eval_DDDDD( real_type x, char const * name ) const {
      Spline const * S = this->getSpline(name);
      return S->DDDDD(x);
    }

    ///@}

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    //! \name Evaluate to std vector
    ///@{

    //!
    //! Evaluate all the splines at `x`
    //!
    void eval( real_type x, vector<real_type> & vals ) const;

    //!
    //! Evaluate the fist derivative of all the splines at `x`
    //!
    void eval_D( real_type x, vector<real_type> & vals ) const;

    //!
    //! Evaluate the second derivative of all the splines at `x`
    //!
    void eval_DD( real_type x, vector<real_type> & vals ) const;

    //!
    //! Evaluate the third derivative of all the splines at `x`
    //!
    void eval_DDD( real_type x, vector<real_type> & vals ) const;

    //!
    //! Evaluate the 4th derivative of all the splines at `x`
    //!
    void eval_DDDD( real_type x, vector<real_type> & vals ) const;

    //!
    //! Evaluate the 5th derivative of all the splines at `x`
    //!
    void eval_DDDDD( real_type x, vector<real_type> & vals ) const;

    ///@}

    //! \name Evaluate to vector
    ///@{

    //!
    //! Evaluate all the splines at `x`
    //!
    void
    eval(
      real_type         x,
      real_type * const vals,
      integer           incy = 1
    ) const;

    //!
    //! Evaluate the fist derivative of all the splines at `x`
    //!
    void
    eval_D(
      real_type         x,
      real_type * const vals,
      integer           incy = 1
    ) const;

    //!
    //! Evaluate the second derivative of all the splines at `x`
    //!
    void
    eval_DD(
      real_type         x,
      real_type * const vals,
      integer           incy = 1
    ) const;

    //!
    //! Evaluate the third derivative of all the splines at `x`
    //!
    void
    eval_DDD(
      real_type         x,
      real_type * const vals,
      integer           incy = 1
    ) const;

    //!
    //! Evaluate the 4th derivative of all the splines at `x`
    //!
    void
    eval_DDDD(
      real_type         x,
      real_type * const vals,
      integer           incy = 1
    ) const;

    //!
    //! Evaluate the 5th derivative of all the splines at `x`
    //!
    void
    eval_DDDDD(
      real_type         x,
      real_type * const vals,
      integer           incy = 1
    ) const;

    ///@}

    //! \name Evaluate using another spline as independent
    ///@{

    //!
    //! Evaluate all the splines at `zeta` using spline[spl] as independent
    //!
    void
    eval2(
      integer             spl,
      real_type           zeta,
      vector<real_type> & vals
    ) const;

    //!
    //! Evaluate the fist derivative of all the splines
    //! at `zeta` using spline[spl] as independent
    //!
    void
    eval2_D(
      integer             spl,
      real_type           zeta,
      vector<real_type> & vals
    ) const;

    //!
    //! Evaluate the second derivative of all the splines
    //! at `zeta` using spline[spl] as independent
    //!
    void
    eval2_DD(
      integer             spl,
      real_type           zeta,
      vector<real_type> & vals
    ) const;

    //!
    //! Evaluate the third derivative of all the splines
    //! at `zeta` using spline[spl] as independent
    //!
    void
    eval2_DDD(
      integer             spl,
      real_type           zeta,
      vector<real_type> & vals
    ) const;

    ///@}

    //! \name Evaluate using another spline as independent
    ///@{

    //!
    //! Evaluate all the splines at `zeta` using spline[spl] as independent
    //!
    void
    eval2(
      integer           spl,
      real_type         zeta,
      real_type * const vals,
      integer           incy = 1
    ) const;

    //!
    //! Evaluate the fist derivative of all the splines
    //! at `zeta` using spline[spl] as independent
    //!
    void
    eval2_D(
      integer           spl,
      real_type         zeta,
      real_type * const vals,
      integer           incy = 1
    ) const;

    //!
    //! Evaluate the second derivative of all the splines
    //! at `zeta` using spline[spl] as independent
    //!
    void
    eval2_DD(
      integer           spl,
      real_type         zeta,
      real_type * const vals,
      integer           incy = 1
    ) const;

    //!
    //! Evaluate the 3rd derivative of all the splines
    //! at `zeta` using spline[spl] as independent
    //!
    void
    eval2_DDD(
      integer           spl,
      real_type         zeta,
      real_type * const vals,
      integer           incy = 1
    ) const;

    ///@}

    // \name Evaluate using another spline as independent
    ///@{

    //!
    //! Evaluate the spline `name` at `zeta` using
    //! spline `indep` as independent
    //!
    real_type
    eval2(
      real_type    zeta,
      char const * indep,
      char const * name
    ) const;

    //!
    //! Evaluate first derivative of the spline `name`
    //! at `zeta` using spline `indep` as independent
    //!
    real_type
    eval2_D(
      real_type    zeta,
      char const * indep,
      char const * name
    ) const;

    //!
    //! Evaluate second derivative of the spline `name`
    //! at `zeta` using spline `indep` as independent
    //!
    real_type
    eval2_DD(
      real_type    zeta,
      char const * indep,
      char const * name
    ) const;

    //!
    //! Evaluate third derivative of the spline `name`
    //! at `zeta` using spline `indep` as independent
    //!
    real_type
    eval2_DDD(
      real_type    zeta,
      char const * indep,
      char const * name
    ) const;

    ///@}

    //! \name Evaluate using another spline as independent
    ///@{

    //!
    //! Evaluate the spline `spl` at `zeta` using
    //! spline `indep` as independent
    //!
    real_type
    eval2(
      real_type zeta,
      integer   indep,
      integer   spl
    ) const;

    //!
    //! Evaluate first derivative of the spline `spl` 
    //! at `zeta` using spline `indep` as independent
    //!
    real_type
    eval2_D(
      real_type zeta,
      integer   indep,
      integer   spl
    ) const;

    //!
    //! Evaluate second derivative of the spline `spl`
    //! at `zeta` using spline `indep` as independent
    //!
    real_type
    eval2_DD(
      real_type zeta,
      integer   indep,
      integer   spl
    ) const;

    //!
    //! Evaluate third derivative of the spline `spl`
    //! at `zeta` using spline `indep` as independent
    //!
    real_type
    eval2_DDD(
      real_type zeta,
      integer   indep,
      integer   spl
    ) const;

    ///@}

    //! \name Evaluate onto a vector
    ///@{

    //! 
    //! Evaluate all the splines at `x`
    //! and fill a map of values in a GenericContainer
    //! 
    void
    eval( real_type x, GenericContainer & vals ) const;

    //! 
    //! Evaluate all the splines at `x` values contained
    //! in vec and fill a map of vector in a GenericContainer
    //! 
    void
    eval( vec_real_type const & vec, GenericContainer & vals ) const;

    //! 
    //! Evaluate all the splines at `x` and fill a map of
    //! values in a GenericContainer with keys in `columns`
    //! 
    void
    eval(
      real_type               x,
      vec_string_type const & columns,
      GenericContainer      & vals
    ) const;

    //! 
    //! Evaluate all the splines at `x` values contained
    //! in vec and fill a map of vector in a GenericContainer
    //! with keys in `columns`
    //! 
    void
    eval(
      vec_real_type   const & vec,
      vec_string_type const & columns,
      GenericContainer      & vals
    ) const;

    ///@}
  
    //! \name Evaluate to GenericContainer using another spline as independent
    ///@{

    //! 
    //! Evaluate all the splines at `zeta` and fill
    //! a map of values in a GenericContainer and `indep`
    //! as independent spline
    //! 
    void
    eval2(
      real_type          zeta,
      integer            indep,
      GenericContainer & vals
    ) const;

    //! 
    //! Evaluate all the splines at `zeta` values contained in
    //! vec and fill a map of vector in a GenericContainer and
    //! `indep` as independent spline
    //! 
    void
    eval2(
      vec_real_type const & zetas,
      integer               indep,
      GenericContainer    & vals
    ) const;

    //! 
    //! Evaluate all the splines at `zeta` and fill a map
    //! of values in a GenericContainer with keys in `columns`
    //! and `indep` as independent spline
    //! 
    void
    eval2(
      real_type               zeta,
      integer                 indep,
      vec_string_type const & columns,
      GenericContainer      & vals
    ) const;

    //! 
    //! Evaluate all the splines at `zeta` values contained
    //! in vec and fill a map of vector in a GenericContainer
    //! with keys in `columns` and `indep` as independent spline
    //! 
    void
    eval2(
      vec_real_type   const & zetas,
      integer                 indep,
      vec_string_type const & columns,
      GenericContainer      & vals
    ) const;

    ///@}
  
    //! \name Evaluate to GenericContainer using another spline as independent
    ///@{

    //! 
    //! Evaluate all the splines at `zeta` and fill a map
    //! of values in a GenericContainer and `indep` as independent spline
    //! 
    void
    eval2(
      real_type          zeta,
      char const       * indep,
      GenericContainer & vals
    ) const {
      this->eval2( zeta, this->getPosition(indep), vals );
    }

    //! 
    //! Evaluate all the splines at `zeta` values contained
    //! in vec and fill a map of vector in a GenericContainer
    //! and `indep` as independent spline
    //! 
    void
    eval2(
      vec_real_type const & zetas,
      char const          * indep,
      GenericContainer    & vals
    ) const {
      this->eval2( zetas, this->getPosition(indep), vals );
    }

    //! 
    //! Evaluate all the splines at `zeta` and fill a map of
    //! values in a GenericContainer with keys in `columns`
    //! and `indep` as independent spline
    //! 
    void
    eval2(
      real_type               zeta,
      char const            * indep,
      vec_string_type const & columns,
      GenericContainer      & vals
    ) const {
      this->eval2( zeta, this->getPosition(indep), columns, vals );
    }

    //! 
    //! Evaluate all the splines at `zeta` values contained
    //! in vec and fill a map of vector in a GenericContainer
    //! with keys in `columns` and `indep` as independent spline
    //! 
    void
    eval2(
      vec_real_type const   & zetas,
      char const            * indep,
      vec_string_type const & columns,
      GenericContainer      & vals
    ) const {
      this->eval2( zetas, this->getPosition(indep), columns, vals );
    }
    ///@}

    //! \name Evaluate derivative into a Generic Container
    ///@{
    
    //! 
    //! Evaluate all the splines at `x` and fill a map of
    //! values in a GenericContainer
    //! 
    void
    eval_D(
      real_type          x,
      GenericContainer & vals
    ) const;

    //! 
    //! Evaluate all the splines at `x` values contained
    //! in vec and fill a map of vector in a GenericContainer
    //! 
    void
    eval_D(
      vec_real_type const & vec,
      GenericContainer    & vals
    ) const;

    //! 
    //! Evaluate all the splines at `x` and fill a map of
    //! values in a GenericContainer with keys in `columns`
    //! 
    void
    eval_D(
      real_type               x,
      vec_string_type const & columns,
      GenericContainer      & vals
    ) const;

    //! 
    //! Evaluate all the splines at `x` values contained
    //! in vec and fill a map of vector in a GenericContainer
    //! with keys in `columns`
    //! 
    void
    eval_D(
      vec_real_type const   & vec,
      vec_string_type const & columns,
      GenericContainer      & vals
    ) const;

    ///@}

  
    //! \name Evaluate derivative to GenericContainer using another spline as independent
    ///@{

    //! 
    //! Evaluate all the splines at `zeta` and fill a map
    //! of values in a GenericContainer and `indep` as independent spline
    //! 
    void
    eval2_D(
      real_type          zeta,
      integer            indep,
      GenericContainer & vals
    ) const;

    //! 
    //! Evaluate all the splines at `zeta` values contained
    //! in vec and fill a map of vector in a GenericContainer
    //! and `indep` as independent spline
    //! 
    void
    eval2_D(
      vec_real_type const & zetas,
      integer               indep,
      GenericContainer    & vals
    ) const;

    //! 
    //! Evaluate all the splines at `zeta` and fill a map
    //! of values in a GenericContainer with keys in `columns`
    //! and `indep` as independent spline
    //! 
    void
    eval2_D(
      real_type               zeta,
      integer                 indep,
      vec_string_type const & columns,
      GenericContainer      & vals
    ) const;

    //! 
    //! Evaluate all the splines at `zeta` values contained
    //! in vec and fill a map of vector in a GenericContainer
    //! with keys in `columns` and `indep` as independent spline
    //! 
    void
    eval2_D(
      vec_real_type   const & zetas,
      integer                 indep,
      vec_string_type const & columns,
      GenericContainer      & vals
    ) const;

    ///@}

    //! \name Evaluate derivative to GenericContainer using another spline as independent
    ///@{

    //! 
    //! Evaluate all the splines at `zeta` and fill a map
    //! of values in a GenericContainer and `indep` as independent spline
    //! 
    void
    eval2_D(
      real_type          zeta,
      char const       * indep,
      GenericContainer & vals
    ) const {
      this->eval2_D( zeta, this->getPosition(indep), vals );
    }

    //!
    //! Evaluate all the splines at `zeta` values contained in 
    //! vec and fill a map of vector in a GenericContainer and 
    //! `indep` as independent spline
    //!
    void
    eval2_D(
      vec_real_type const & zetas,
      char const          * indep,
      GenericContainer    & vals
    ) const {
      this->eval2_D( zetas, this->getPosition(indep), vals );
    }

    //! 
    //! Evaluate all the splines at `zeta` and fill a map of
    //! values in a GenericContainer with keys in `columns`
    //! and `indep` as independent spline
    //! 
    void
    eval2_D(
      real_type               zeta,
      char const            * indep,
      vec_string_type const & columns,
      GenericContainer      & vals
    ) const {
      this->eval2_D( zeta, this->getPosition(indep), columns, vals );
    }

    //! 
    //! Evaluate all the splines at `zeta` values contained
    //! in vec and fill a map of vector in a GenericContainer
    //! with keys in `columns` and `indep` as independent spline
    //! 
    void
    eval2_D(
      vec_real_type const   & zetas,
      char const            * indep,
      vec_string_type const & columns,
      GenericContainer      & vals
    ) const {
      this->eval2_D( zetas, this->getPosition(indep), columns, vals );
    }

    ///@}

    //! \name Evaluate second derivative to GenericContainer using another spline as independent
    ///@{

    //! 
    //! Evaluate all the splines at `x` and fill a map of
    //! values in a GenericContainer
    //! 
    void
    eval_DD(
      real_type          x,
      GenericContainer & vals
    ) const;

    //! 
    //! Evaluate all the splines at `x` values contained in vec
    //! and fill a map of vector in a GenericContainer
    //! 
    void
    eval_DD(
      vec_real_type const & vec,
      GenericContainer    & vals
    ) const;

    //! 
    //! Evaluate all the splines at `x` and fill a map of
    //! values in a GenericContainer with keys in `columns`
    //! 
    void
    eval_DD(
      real_type               x,
      vec_string_type const & columns,
      GenericContainer      & vals
    ) const;

    //! 
    //! Evaluate all the splines at `x` values contained in vec and
    //! fill a map of vector in a GenericContainer with keys in `columns`
    //! 
    void
    eval_DD(
      vec_real_type   const & vec,
      vec_string_type const & columns,
      GenericContainer      & vals
    ) const;

    //! 
    //! Evaluate all the splines at `zeta` and fill a map of
    //! values in a GenericContainer and `indep` as independent spline
    //! 
    void
    eval2_DD(
      real_type          zeta,
      integer            indep,
      GenericContainer & vals
    ) const;

    //! 
    //! Evaluate all the splines at `zeta` values contained in vec
    //! and fill a map of vector in a GenericContainer
    //! and `indep` as independent spline
    //! 
    void
    eval2_DD(
      vec_real_type const & zetas,
      integer               indep,
      GenericContainer    & vals
    ) const;

    //! 
    //! Evaluate all the splines at `zeta` and fill a map
    //! of values in a GenericContainer with keys in `columns`
    //! and `indep` as independent spline
    //! 
    void
    eval2_DD(
      real_type               zeta,
      integer                 indep,
      vec_string_type const & columns,
      GenericContainer      & vals
    ) const;

    //! 
    //! Evaluate all the splines at `zeta` values contained
    //! in vec and fill a map of vector in a GenericContainer
    //! with keys in `columns` and `indep` as independent spline
    //! 
    void
    eval2_DD(
      vec_real_type   const & zetas,
      integer                 indep,
      vec_string_type const & columns,
      GenericContainer      & vals
    ) const;

    ///@}

    //! \name Evaluate second derivative to GenericContainer using another spline as independent
    ///@{

    //! 
    //! Evaluate all the splines at `zeta` and fill a map of
    //! values in a GenericContainer and `indep` as independent spline
    //! 
    void
    eval2_DD(
      real_type          zeta,
      char const       * indep,
      GenericContainer & vals
    ) const {
      this->eval2_DD( zeta, this->getPosition(indep), vals );
    }

    //! 
    //! Evaluate all the splines at `zeta` values contained in vec
    //! and fill a map of vector in a GenericContainer
    //! and `indep` as independent spline
    //! 
    void
    eval2_DD(
      vec_real_type const & zetas,
      char const          * indep,
      GenericContainer    & vals
    ) const {
      this->eval2_DD( zetas, this->getPosition(indep), vals );
    }

    //! 
    //! Evaluate all the splines at `zeta` and fill a map
    //! of values in a GenericContainer with keys in `columns`
    //! and `indep` as independent spline
    //! 
    void
    eval2_DD(
      real_type               zeta,
      char const            * indep,
      vec_string_type const & columns,
      GenericContainer      & vals
    ) const {
      this->eval2_DD( zeta, this->getPosition(indep), columns, vals );
    }

    //! 
    //! Evaluate all the splines at `zeta` values contained
    //! in vec and fill a map of vector in a GenericContainer
    //! with keys in `columns` and `indep` as independent spline
    //! 
    void
    eval2_DD(
      vec_real_type   const & zetas,
      char            const * indep,
      vec_string_type const & columns,
      GenericContainer      & vals
    ) const {
      this->eval2_DD( zetas, this->getPosition(indep), columns, vals );
    }
    ///@}

    //! \name Evaluate third derivative to GenericContainer
    ///@{

    //! 
    //! Evaluate all the splines at `x` and fill a map
    //! of values in a GenericContainer
    //! 
    void
    eval_DDD(
      real_type          x,
      GenericContainer & vals
    ) const;

    //! 
    //! Evaluate all the splines at `x` values contained
    //! in vec and fill a map of vector in a GenericContainer
    //! 
    void
    eval_DDD(
      vec_real_type const & vec,
      GenericContainer    & vals
    ) const;

    //! 
    //! Evaluate all the splines at `x` and fill a map of
    //! values in a GenericContainer with keys in `columns`
    //! 
    void
    eval_DDD(
      real_type               x,
      vec_string_type const & columns,
      GenericContainer      & vals
    ) const;

    //! 
    //! Evaluate all the splines at `x` values contained in vec
    //! and fill a map of vector in a GenericContainer with keys in `columns`
    //! 
    void
    eval_DDD(
      vec_real_type   const & vec,
      vec_string_type const & columns,
      GenericContainer      & vals
    ) const;
    ///@}

    //! \name Evaluate third derivative to GenericContainer using another spline as independent
    ///@{

    //! 
    //! Evaluate all the splines at `zeta` and fill a map of values
    //! in a GenericContainer and `indep` as independent spline
    //! 
    void
    eval2_DDD(
      real_type          zeta,
      integer            indep,
      GenericContainer & vals
    ) const;

    //! 
    //! Evaluate all the splines at `zeta` values contained in vec
    //! and fill a map of vector in a GenericContainer
    //! and `indep` as independent spline
    //! 
    void
    eval2_DDD(
      vec_real_type const & zetas,
      integer               indep,
      GenericContainer    & vals
    ) const;

    //! 
    //! Evaluate all the splines at `zeta` and fill a map of values
    //! in a GenericContainer with keys in `columns`
    //! and `indep` as independent spline
    //! 
    void
    eval2_DDD(
      real_type               zeta,
      integer                 indep,
      vec_string_type const & columns,
      GenericContainer      & vals
    ) const;

    //! 
    //! Evaluate all the splines at `zeta` values contained
    //! in vec and fill a map of vector in a GenericContainer
    //! with keys in `columns` and `indep` as independent spline
    //! 
    void
    eval2_DDD(
      vec_real_type   const & zetas,
      integer                 indep,
      vec_string_type const & columns,
      GenericContainer      & vals
    ) const;
    ///@}

    //! \name Evaluate third derivative to GenericContainer using another spline as independent
    ///@{

    //! 
    //! Evaluate all the splines at `zeta` and fill a map
    //! of values in a GenericContainer and `indep` as independent spline
    //! 
    void
    eval2_DDD(
      real_type          zeta,
      char const       * indep,
      GenericContainer & vals
    ) const {
      this->eval2_DDD( zeta, this->getPosition(indep), vals );
    }

    //! 
    //! Evaluate all the splines at `zeta` values contained
    //! in vec and fill a map of vector in a GenericContainer
    //! and `indep` as independent spline
    //! 
    void
    eval2_DDD(
      vec_real_type const & zetas,
      char          const * indep,
      GenericContainer    & vals
    ) const {
      this->eval2_DDD( zetas, this->getPosition(indep), vals );
    }

    //! 
    //! Evaluate all the splines at `zeta` and fill a map of
    //! values in a GenericContainer with keys in `columns`
    //! and `indep` as independent spline
    //! 
    void
    eval2_DDD(
      real_type               zeta,
      char            const * indep,
      vec_string_type const & columns,
      GenericContainer      & vals
    ) const {
      this->eval2_DDD( zeta, this->getPosition(indep), columns, vals );
    }

    //! 
    //! Evaluate all the splines at `zeta` values contained
    //! in vec and fill a map of vector in a GenericContainer
    //! with keys in `columns` and `indep` as independent spline
    //! 
    void
    eval2_DDD(
      vec_real_type   const & zetas,
      char            const * indep,
      vec_string_type const & columns,
      GenericContainer      & vals
    ) const {
      this->eval2_DDD( zetas, this->getPosition(indep), columns, vals );
    }

    ///@}

    //! \name Build Spline
    ///@{

    ///////////////////////////////////////////////////////////////////////////
    //! 
    //! Build a set of splines
    //! 
    //! \param nspl    the number of splines
    //! \param npts    the number of points of each splines
    //! \param headers the names of the splines
    //! \param stype   the type of each spline
    //! \param X       pointer to X independent values
    //! \param Y       vector of `nspl` pointers to Y depentendent values.
    //! \param Yp      vector of `nspl` pointers to Y derivative depentendent values.
    //! 
    void
    build(
      integer               nspl,
      integer               npts,
      char         const ** headers,
      SplineType1D const *  stype,
      real_type    const *  X,
      real_type    const ** Y,
      real_type    const ** Yp = nullptr
    );

    void
    setup( GenericContainer const & gc );

    void
    build( GenericContainer const & gc )
    { this->setup(gc); }

    ///@}

    //! Return spline type (as number)
    unsigned
    type() const
    { return SPLINE_SET_TYPE; }

    string
    info() const;

    void
    info( ostream_type & stream ) const
    { stream << this->info() << '\n'; }

    void
    dump_table( ostream_type & s, integer num_points ) const;

  };

}

// EOF: SplineSet.hxx
