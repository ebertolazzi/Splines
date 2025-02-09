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
 |   _     _                       ____        _ _
 |  | |   (_)_ __   ___  __ _ _ __/ ___| _ __ | (_)_ __   ___
 |  | |   | | '_ \ / _ \/ _` | '__\___ \| '_ \| | | '_ \ / _ \
 |  | |___| | | | |  __/ (_| | |   ___) | |_) | | | | | |  __/
 |  |_____|_|_| |_|\___|\__,_|_|  |____/| .__/|_|_|_| |_|\___|
 |                                      |_|
\*/

namespace Splines {

  //! Linear spline class
  class LinearSpline : public Spline {
    Malloc_real m_mem_linear;
    bool        m_external_alloc{false};

  public:

    #ifndef DOXYGEN_SHOULD_SKIP_THIS
    using Spline::build;
    #endif

    //!
    //! Build an empty spline of `LinearSpline` type
    //!
    //! \param name the name of the spline
    //!
    LinearSpline( string_view name = "LinearSpline" );

    //!
    //! Spline destructor.
    //!
    ~LinearSpline() override {}

    //! Use externally allocated memory for `npts` points
    void
    reserve_external(
      integer      n,
      real_type *& p_x,
      real_type *& p_y
    );

    // --------------------------- VIRTUALS -----------------------------------

    real_type eval( real_type x ) const override;
    real_type D( real_type x ) const override;
    real_type DD( real_type ) const override { return 0; }
    real_type DDD( real_type ) const override { return 0; }

    real_type id_eval( integer ni, real_type x ) const override;
    real_type id_D( integer, real_type ) const override;
    real_type id_DD( integer, real_type ) const override { return 0; }
    real_type id_DDD( integer, real_type ) const override { return 0; }

    void write_to_stream( ostream_type & s ) const override;
    SplineType1D type() const override { return SplineType1D::LINEAR; }

    // --------------------------- VIRTUALS -----------------------------------

    void reserve( integer npts ) override;
    void build() override {}
    void clear() override;

    integer // order
    coeffs(
      real_type cfs[],
      real_type nodes[],
      bool      transpose = false
    ) const override;

    integer order() const override;
    void setup( GenericContainer const & gc ) override;

    void
    y_min_max(
      integer   & i_min_pos,
      real_type & x_min_pos,
      real_type & y_min,
      integer   & i_max_pos,
      real_type & x_max_pos,
      real_type & y_max
    ) const override;

    void
    y_min_max(
      vector<integer>   & i_min_pos,
      vector<real_type> & x_min_pos,
      vector<real_type> & y_min,
      vector<integer>   & i_max_pos,
      vector<real_type> & x_max_pos,
      vector<real_type> & y_max
    ) const override;

  };

}

// EOF: SplineLinbear.hxx
