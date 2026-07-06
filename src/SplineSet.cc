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

#include "Splines/Splines.hh"

namespace Splines
{

  void SplineSet::setup( GenericContainer const & gc )
  {
    string const where = fmt::format( "SplineSet[{}]::setup( gc ): ", m_name );

    // ===================================================================
    // FASE 1: Raccolta e tracciamento delle chiavi disponibili
    // ===================================================================
    std::set<std::string> keywords;
    for ( auto const & pair : gc.get_map( where ) ) { keywords.insert( pair.first ); }

    // Strutture dati temporanee per la costruzione
    vec_string_type       spline_type_vec;  // Tipi di spline come stringhe
    vec_real_type         X;                // Nodi x condivisi
    vector<SplineType1D>  stype;            // Tipi di spline convertiti (enum)
    vec_string_type       headers;          // Nomi delle spline
    vector<vec_real_type> Y, Yp;            // Valori y e derivate

    // ===================================================================
    // FASE 2: Lettura e validazione dei campi obbligatori
    // ===================================================================

    // --- Lettura "spline_type" (obbligatorio) ---
    GenericContainer const & gc_stype = gc( "spline_type", where );
    keywords.erase( "spline_type" );

    // --- Lettura "xdata" (obbligatorio) ---
    GenericContainer const & gc_xdata = gc( "xdata", where );
    keywords.erase( "xdata" );

    // --- Lettura "ydata" (obbligatorio) ---
    GenericContainer const & gc_ydata = gc( "ydata", where );
    keywords.erase( "ydata" );

    // ===================================================================
    // FASE 3: Parsing e conversione dei tipi di spline
    // ===================================================================
    gc_stype.copyto_vec_string( spline_type_vec, ( where + "in reading 'spline_type'\n" ) );
    m_nspl = static_cast<integer>( spline_type_vec.size() );

    // Validazione: almeno una spline richiesta
    UTILS_ASSERT( m_nspl > 0, "{} expected at least one spline, got nspl = {}\n", where, m_nspl );

    // Conversione da stringhe a enum SplineType1D
    stype.resize( m_nspl );
    for ( integer spl = 0; spl < m_nspl; ++spl )
    {
      // NOTA: string_to_splineType1D lancia eccezione se tipo non riconosciuto
      stype[spl] = string_to_splineType1D( spline_type_vec[spl] );
    }

    // ===================================================================
    // FASE 4: Lettura headers
    // ===================================================================
    // Per ydata non-MAP gli headers sono obbligatori.
    // Per ydata come MAP gli headers sono opzionali ma, se presenti,
    // definiscono l'ordine stabile con cui associare tipi, dati e boundary.
    bool const has_explicit_headers = gc.exists( "headers" );
    if ( has_explicit_headers )
    {
      GenericContainer const & gc_headers = gc( "headers", where );
      keywords.erase( "headers" );
      gc_headers.copyto_vec_string( headers, ( where + ", reading 'headers'" ) );

      // Validazione: numero di headers deve corrispondere al numero di spline
      UTILS_ASSERT(
        headers.size() == static_cast<size_t>( m_nspl ),
        "{} field 'headers' expected to be of size {} found of size {}\n",
        where,
        m_nspl,
        headers.size() );

      // Validazione: nessun header può essere vuoto o null
      for ( integer spl = 0; spl < m_nspl; ++spl )
      {
        UTILS_ASSERT( !headers[spl].empty(), "{} header[{}] cannot be empty\n", where, spl );
      }
    }
    else if ( GC_type::MAP != gc_ydata.get_type() )
    {
      UTILS_ERROR( "{} missing required field 'headers'\n", where );
    }

    // ===================================================================
    // FASE 5: Lettura e validazione dei dati X (condivisi)
    // ===================================================================
    gc_xdata.copyto_vec_real( X, ( where + "reading 'xdata'" ) );
    m_npts = static_cast<integer>( X.size() );

    // Validazione: almeno 2 punti necessari per una spline
    UTILS_ASSERT( m_npts > 1, "{} expected at least 2 points, got npts = {}\n", where, m_npts );

    // Validazione: X deve essere strettamente crescente
    for ( integer i = 1; i < m_npts; ++i )
    {
      UTILS_ASSERT(
        X[i] > X[i - 1],
        "{} xdata must be strictly increasing: X[{}] = {} <= X[{}] = {}\n",
        where,
        i - 1,
        X[i - 1],
        i,
        X[i] );
    }

    // Pre-allocazione per le Y e Yp di tutte le spline
    Y.resize( m_nspl );
    Yp.resize( m_nspl );

    // ===================================================================
    // FASE 6: Lettura dei dati Y (supporta formati multipli)
    // ===================================================================
    switch ( gc_ydata.get_type() )
    {
      // ---------------------------------------------------------------
      // CASO 1: Vettore singolo → Una sola spline
      // ---------------------------------------------------------------
      case GC_type::VEC_BOOL:
      case GC_type::VEC_INTEGER:
      case GC_type::VEC_LONG:
      case GC_type::VEC_REAL:
      case GC_type::VEC_COMPLEX:
      {
        vec_real_type data;
        gc_ydata.copyto_vec_real( data, ( where + "reading 'ydata'" ) );

        // Deve esserci esattamente una spline
        UTILS_ASSERT(
          m_nspl == 1,
          "{} number of splines [{}] is incompatible with ydata as vector (expected 1 spline)\n",
          where,
          m_nspl );

        // Numero di punti deve corrispondere (o n-1 per CONSTANT)
        integer expected_npts = ( stype[0] == SplineType1D::CONSTANT ) ? m_npts - 1 : m_npts;
        UTILS_ASSERT(
          static_cast<size_t>( expected_npts ) == data.size(),
          "{} number of points [{}] differs from ydata size [{}] for spline type '{}'\n",
          where,
          expected_npts,
          data.size(),
          spline_type_vec[0] );

        // Usa move semantics invece di copy
        Y[0] = std::move( data );
      }
      break;

      // ---------------------------------------------------------------
      // CASO 2: Matrice → Ogni colonna è una spline
      // ---------------------------------------------------------------
      case GC_type::MAT_INTEGER:
      case GC_type::MAT_LONG:
      case GC_type::MAT_REAL:
      {
        mat_real_type data;
        gc_ydata.copyto_mat_real( data, ( where + "reading 'ydata'" ) );

        // Numero colonne deve corrispondere al numero di spline
        UTILS_ASSERT(
          static_cast<size_t>( m_nspl ) == data.num_cols(),
          "{} number of splines [{}] differs from ydata columns [{}]\n",
          where,
          m_nspl,
          data.num_cols() );

        // Numero righe deve corrispondere al numero di punti
        UTILS_ASSERT(
          static_cast<size_t>( m_npts ) == data.num_rows(),
          "{} number of points [{}] differs from ydata rows [{}]\n",
          where,
          m_npts,
          data.num_rows() );

        // Estrazione delle colonne
        for ( integer i = 0; i < m_nspl; ++i ) { data.get_column( static_cast<unsigned>( i ), Y[i] ); }
      }
      break;

      // ---------------------------------------------------------------
      // CASO 3: Vettore di vettori → Array di spline
      // ---------------------------------------------------------------
      case GC_type::VECTOR:
      {
        vector_type const & data = gc_ydata.get_vector();

        // Dimensione del vettore deve corrispondere al numero di spline
        UTILS_ASSERT(
          static_cast<size_t>( m_nspl ) == data.size(),
          "{} field 'ydata' expected of size {} found of size {}\n",
          where,
          m_nspl,
          data.size() );

        string msg1 = where + " reading 'ydata' columns";
        for ( integer spl = 0; spl < m_nspl; ++spl )
        {
          GenericContainer const & datai = data[spl];

          // Calcola numero atteso di righe (n-1 per CONSTANT, n altrimenti)
          integer expected_npts = m_npts;
          if ( stype[spl] == SplineType1D::CONSTANT ) { --expected_npts; }

          datai.copyto_vec_real( Y[spl], msg1 );

          // Validazione dimensione per questa spline
          UTILS_ASSERT(
            static_cast<size_t>( expected_npts ) == Y[spl].size(),
            "{} column {} of 'ydata' of type '{}' expected of size {} found of size {}\n",
            where,
            spl,
            spline_type_vec[spl],
            expected_npts,
            Y[spl].size() );
        }
      }
      break;

      // ---------------------------------------------------------------
      // CASO 4: Mappa → Chiavi = nomi spline, valori = dati
      // ---------------------------------------------------------------
      case GC_type::MAP:
      {
        map_type const & data = gc_ydata.get_map();

        // Dimensione mappa deve corrispondere al numero di spline
        UTILS_ASSERT(
          data.size() == static_cast<size_t>( m_nspl ),
          "{} field 'ydata' expected of size {} found of size {}\n",
          where,
          m_nspl,
          data.size() );

        string  msg1 = where + " reading 'ydata' columns";
        if ( headers.empty() )
        {
          headers.clear();
          headers.reserve( data.size() );

          integer spl = 0;
          for ( auto const & [name, container] : data )
          {
            headers.emplace_back( name );
            UTILS_ASSERT( !name.empty(), "{} spline name cannot be empty in ydata map\n", where );

            integer expected_npts = m_npts;
            if ( stype[spl] == SplineType1D::CONSTANT ) { --expected_npts; }

            container.copyto_vec_real( Y[spl], msg1 );

            UTILS_ASSERT(
              static_cast<size_t>( expected_npts ) == Y[spl].size(),
              "{} column '{}' of 'ydata' of type '{}' expected of size {} found of size {}\n",
              where,
              name,
              spline_type_vec[spl],
              expected_npts,
              Y[spl].size() );

            ++spl;
          }
        }
        else
        {
          for ( integer spl = 0; spl < m_nspl; ++spl )
          {
            string const & name = headers[spl];
            auto const     it   = data.find( name );

            UTILS_ASSERT( !name.empty(), "{} header[{}] cannot be empty\n", where, spl );
            UTILS_ASSERT( it != data.end(), "{} column '{}' listed in 'headers' was not found in 'ydata'\n", where, name );

            integer expected_npts = m_npts;
            if ( stype[spl] == SplineType1D::CONSTANT ) { --expected_npts; }

            it->second.copyto_vec_real( Y[spl], msg1 );

            UTILS_ASSERT(
              static_cast<size_t>( expected_npts ) == Y[spl].size(),
              "{} column '{}' of 'ydata' of type '{}' expected of size {} found of size {}\n",
              where,
              name,
              spline_type_vec[spl],
              expected_npts,
              Y[spl].size() );
          }
        }
      }
      break;

      // ---------------------------------------------------------------
      // CASO DEFAULT: Tipo non supportato
      // ---------------------------------------------------------------
      default:
        UTILS_ERROR(
          "{} field 'ydata' expected to be of type:\n"
          "  - vec_[int/long/real]_type (single spline)\n"
          "  - mat_[int/long/real]_type (matrix of splines)\n"
          "  - vector_type (vector of spline data)\n"
          "  - map_type (named splines)\n"
          "Found: '{}'\n",
          where,
          gc_ydata.get_type_name() );
        break;
    }

    // ===================================================================
    // FASE 7: Lettura derivate (opzionale, per spline Hermite)
    // ===================================================================
    if ( gc.exists( "ypdata" ) )
    {
      GenericContainer const & gc_ypdata = gc( "ypdata", where );
      keywords.erase( "ypdata" );

      // Le derivate devono essere fornite come MAP (nome → valori)
      UTILS_ASSERT(
        GC_type::MAP == gc_ypdata.get_type(),
        "{} field 'ypdata' expected to be of type 'map_type' found: '{}'\n",
        where,
        gc_ypdata.get_type_name() );

      // ---------------------------------------------------------------
      // Costruzione della mappa nome → indice (usa unordered_map)
      // ---------------------------------------------------------------
      std::unordered_map<string, integer> h_to_pos;
      h_to_pos.reserve( headers.size() );

      for ( integer idx = 0; idx < static_cast<integer>( headers.size() ); ++idx ) { h_to_pos[headers[idx]] = idx; }

      string           msg1 = where + " reading 'ypdata' columns";
      map_type const & data = gc_ypdata.get_map();

      // Per ogni derivata fornita nella mappa
      for ( auto const & [name, container] : data )
      {
        // Cerca la posizione corrispondente negli headers
        auto it_pos = h_to_pos.find( name );
        UTILS_ASSERT(
          it_pos != h_to_pos.end(),
          "{} column '{}' of 'ypdata' does not match any spline name in headers\n",
          where,
          name );

        integer spl = it_pos->second;

        // Calcola numero atteso di punti
        integer expected_npts = m_npts;
        if ( stype[spl] == SplineType1D::CONSTANT ) { --expected_npts; }

        container.copyto_vec_real( Yp[spl], msg1 );

        // FIXED: Confronta con Yp[spl].size() invece di Y[spl].size()
        // FIXED: Typo "or type" → "of type"
        UTILS_ASSERT(
          static_cast<size_t>( expected_npts ) == Yp[spl].size(),
          "{} column '{}' of 'ypdata' of type '{}' expected of size {} found of size {}\n",
          where,
          name,
          spline_type_vec[spl],
          expected_npts,
          Yp[spl].size() );
      }
    }

    // ===================================================================
    // FASE 8: Preparazione puntatori per chiamata a build()
    // ===================================================================
    // NOTA: build() richiede array di puntatori C-style (char**, real_type**)
    // Questa sezione prepara questi array dai vettori C++ moderni

    Utils::Malloc<void *> mem( where );
    mem.allocate( 3 * m_nspl );

    void ** pp__headers = mem( m_nspl );
    void ** pp__Y       = mem( m_nspl );
    void ** pp__Yp      = mem( m_nspl );

    // Conversione da vector<string>/vector<vector<real_type>> a puntatori C
    // ATTENZIONE: headers deve rimanere valido fino alla fine di build()
    for ( integer spl = 0; spl < m_nspl; ++spl )
    {
      // Puntatore al nome della spline (const char*)
      // SICUREZZA: headers[spl] deve rimanere valido durante build()
      pp__headers[spl] = const_cast<void *>( reinterpret_cast<void const *>( headers[spl].c_str() ) );

      // Puntatore ai dati Y
      pp__Y[spl] = reinterpret_cast<void *>( Y[spl].data() );

      // Puntatore ai dati Yp (nullptr se non forniti)
      pp__Yp[spl] = reinterpret_cast<void *>( Yp[spl].empty() ? nullptr : Yp[spl].data() );
    }

    // ===================================================================
    // FASE 9: Chiamata a build() con i dati preparati
    // ===================================================================
    // NOTA: I cast multipli sono necessari per l'interfaccia C-style di build()
    // ma sono sicuri perché i dati sorgente rimangono validi per tutta la chiamata
    this->build(
      m_nspl,
      m_npts,
      reinterpret_cast<char const **>( const_cast<void const **>( pp__headers ) ),
      stype.data(),
      X.data(),
      reinterpret_cast<real_type const **>( const_cast<void const **>( pp__Y ) ),
      reinterpret_cast<real_type const **>( const_cast<void const **>( pp__Yp ) ) );

    // mem viene deallocato automaticamente qui (RAII)

    // ===================================================================
    // FASE 10: Configurazione condizioni al contorno (opzionale)
    // ===================================================================
    if ( gc.exists( "boundary" ) )
    {
      GenericContainer const & gc_boundary = gc( "boundary" );
      keywords.erase( "boundary" );

      integer ne = static_cast<integer>( gc_boundary.get_num_elements() );

      // Validazione: una configurazione boundary per ogni spline
      // FIXED: Ordine corretto dei parametri nel messaggio
      UTILS_ASSERT(
        ne == m_nspl,
        "{} field 'boundary' expected a generic vector of size: {} but is of size: {}\n",
        where,
        m_nspl,
        ne );

      // Configurazione boundary per ogni spline
      for ( integer ispl = 0; ispl < ne; ++ispl )
      {
        Spline *                 S    = m_splines[ispl].get();
        GenericContainer const & item = gc_boundary( ispl, "SplineSet boundary data" );

        // Raccolta keywords per questo elemento boundary
        std::set<std::string> keywords2;
        for ( auto const & pair : item.get_map( where ) ) { keywords2.insert( pair.first ); }

        // ---------------------------------------------------------------
        // Configurazione closed/open
        // ---------------------------------------------------------------
        bool is_closed = false;
        item.get_if_exists( "closed", is_closed );
        keywords2.erase( "closed" );

        if ( is_closed )
        {
          // Spline periodica (endpoint = startpoint)
          S->make_closed();
        }
        else
        {
          // Spline aperta (gestione extrapolation)
          keywords2.erase( "extend" );
          keywords2.erase( "can_extend" );
          S->make_opened();

          bool can_extend = false;
          // Supporta sia "extend" che "can_extend" come alias
          if ( !item.get_if_exists( "extend", can_extend ) ) { item.get_if_exists( "can_extend", can_extend ); }

          if ( can_extend )
          {
            // Permette valutazione oltre il range dei dati
            S->make_unbounded();

            bool extend_constant = false;
            keywords2.erase( "extend_constant" );
            item.get_if_exists( "extend_constant", extend_constant );

            if ( extend_constant )
            {
              // Extrapolazione costante (usa valore endpoint)
              S->make_extended_constant();
            }
            else
            {
              // Extrapolazione con continuità della derivata
              S->make_extended_not_constant();
            }
          }
          else
          {
            // Non permette valutazione oltre il range
            S->make_bounded();
          }
        }

        // ---------------------------------------------------------------
        // Warning per chiavi non utilizzate in questo boundary
        // ---------------------------------------------------------------
        if ( !keywords2.empty() )
        {
          // Usa std::accumulate per costruire la stringa in modo efficiente
          string unused_keys{ std::accumulate(
            keywords2.begin(),
            keywords2.end(),
            string{},
            []( string acc, string const & key ) { return acc.empty() ? key : acc + " " + key; } ) };

          UTILS_WARNING( false, "{} spline N.{} of {} has unused boundary keys: {}\n", where, ispl, ne, unused_keys );
        }
      }
    }

    // ===================================================================
    // FASE 11: Warning per chiavi non utilizzate nel GC principale
    // ===================================================================
    if ( !keywords.empty() )
    {
      // Usa std::accumulate per costruire la stringa in modo efficiente
      string unused_keys{ std::accumulate(
        keywords.begin(),
        keywords.end(),
        string{},
        []( string acc, string const & key ) { return acc.empty() ? key : acc + " " + key; } ) };

      UTILS_WARNING( false, "{} unused configuration keys: {}\n", where, unused_keys );
    }
  }

  void SplineSet::build(
    integer const           nspl,
    integer const           npts,
    char const * const      headers[],
    SplineType1D const      stype[],
    real_type const         X[],
    real_type const * const Y[],
    real_type const * const Yp[] )
  {
    string const msg = fmt::format( "SplineSet[{}]::build(...):", m_name );

    // Validazione input
    UTILS_ASSERT( nspl > 0, "{} expected positive nspl = {}\n", msg, nspl );
    UTILS_ASSERT( npts > 1, "{} expected npts = {} greater than 1", msg, npts );
    UTILS_ASSERT( headers != nullptr, "{} headers array is null\n", msg );
    UTILS_ASSERT( stype != nullptr, "{} stype array is null\n", msg );
    UTILS_ASSERT( X != nullptr, "{} X array is null\n", msg );
    UTILS_ASSERT( Y != nullptr, "{} Y array is null\n", msg );

    m_nspl = nspl;
    m_npts = npts;

    // Allocazione memoria per le strutture principali
    m_splines.resize( m_nspl );
    m_is_monotone = m_mem_int.realloc( m_nspl );
    m_header_to_position.clear();

    // ===================================================================
    // FASE 1: Calcolo della memoria totale necessaria
    // ===================================================================
    integer mem = npts;  // Spazio base per X

    for ( integer spl = 0; spl < nspl; ++spl )
    {
      UTILS_ASSERT( headers[spl] != nullptr, "{} headers[{}] array is null\n", msg, spl );

      switch ( stype[spl] )
      {
        // Spline quintiche: necessitano Y, Yp, Ypp (3 * npts)
        case SplineType1D::QUINTIC_CUBIC:
        case SplineType1D::QUINTIC_AKIMA:
        case SplineType1D::QUINTIC_VANLEER:
        case SplineType1D::QUINTIC_PCHIP:
          mem += npts;  // Spazio per Ypp
          [[fallthrough]];

        // Spline cubiche/hermitiane: necessitano Y, Yp (2 * npts)
        case SplineType1D::CUBIC:
        case SplineType1D::AKIMA:
        case SplineType1D::VANLEER:
        case SplineType1D::PCHIP:
        case SplineType1D::HERMITE:
          mem += npts;  // Spazio per Yp
          [[fallthrough]];

        // Spline lineari/costanti: necessitano solo Y (npts)
        case SplineType1D::CONSTANT:
        case SplineType1D::LINEAR:
          mem += npts;  // Spazio per Y
          break;

        // Tipi non supportati
        case SplineType1D::SPLINE_SET:
        case SplineType1D::SPLINE_VEC:
          UTILS_ERROR(
            "{} At spline n.{} named {} cannot be done for type = {}\n",
            msg,
            spl,
            headers[spl],
            to_string( stype[spl] ) );
      }
    }

    // ===================================================================
    // FASE 2: Allocazione memoria
    // ===================================================================
    m_mem.reallocate( mem + 2 * nspl );
    m_mem_p.reallocate( 3 * nspl );  // Per Y, Yp, Ypp pointers

    // Assegnazione dei puntatori alle sezioni di memoria
    m_Y    = m_mem_p( m_nspl );
    m_Yp   = m_mem_p( m_nspl );
    m_Ypp  = m_mem_p( m_nspl );
    m_X    = m_mem( m_npts );
    m_Ymin = m_mem( m_nspl );
    m_Ymax = m_mem( m_nspl );

    // Copia dei nodi X (condivisi da tutte le spline)
    if ( npts > 0 ) std::memcpy( m_X, X, npts * sizeof( *m_X ) );

    // ===================================================================
    // FASE 3: Costruzione delle singole spline
    // ===================================================================
    for ( integer spl = 0; spl < nspl; ++spl )
    {
      UTILS_ASSERT( Y[spl] != nullptr, "{} Y[{}] array is null\n", msg, spl );

      // Riferimenti ai puntatori per questa spline
      real_type *& pY   = m_Y[spl];
      real_type *& pYp  = m_Yp[spl];
      real_type *& pYpp = m_Ypp[spl];

      // Allocazione e copia dei valori Y
      pY = m_mem( m_npts );
      integer y_npts = npts;
      if ( stype[spl] == SplineType1D::CONSTANT ) --y_npts;
      if ( y_npts > 0 ) std::memcpy( pY, Y[spl], y_npts * sizeof( *pY ) );
      if ( stype[spl] == SplineType1D::CONSTANT ) pY[npts - 1] = pY[y_npts - 1];

      // ---------------------------------------------------------------
      // Calcolo min/max per questa spline
      // ---------------------------------------------------------------
      // NOTA: Per CONSTANT si esclude l'ultimo punto (comportamento corretto?)
      if ( stype[spl] == SplineType1D::CONSTANT )
      {
        m_Ymin[spl] = *std::min_element( pY, pY + npts - 1 );
        m_Ymax[spl] = *std::max_element( pY, pY + npts - 1 );
      }
      else
      {
        m_Ymin[spl] = *std::min_element( pY, pY + npts );
        m_Ymax[spl] = *std::max_element( pY, pY + npts );
      }

      // Inizializzazione puntatori derivate
      pYpp = pYp = nullptr;

      // ---------------------------------------------------------------
      // Allocazione memoria per derivate se necessario
      // ---------------------------------------------------------------
      switch ( stype[spl] )
      {
        case SplineType1D::QUINTIC_CUBIC:
        case SplineType1D::QUINTIC_AKIMA:
        case SplineType1D::QUINTIC_VANLEER:
        case SplineType1D::QUINTIC_PCHIP:
          pYpp = m_mem( m_npts );  // Alloca seconda derivata
          [[fallthrough]];

        case SplineType1D::CUBIC:
        case SplineType1D::AKIMA:
        case SplineType1D::VANLEER:
        case SplineType1D::PCHIP:
        case SplineType1D::HERMITE:
          pYp = m_mem( m_npts );  // Alloca prima derivata

          // Per HERMITE le derivate devono essere fornite dall'utente
          if ( stype[spl] == SplineType1D::HERMITE )
          {
            UTILS_ASSERT(
              Yp != nullptr && Yp[spl] != nullptr,
              "{} At spline n.{} named {}\n"
              "expect to find derivative values",
              msg,
              spl,
              headers[spl] );
            if ( npts > 0 ) std::memcpy( pYp, Yp[spl], npts * sizeof( *pYp ) );
          }
          [[fallthrough]];

        case SplineType1D::CONSTANT:
        case SplineType1D::LINEAR:
        case SplineType1D::SPLINE_SET:
        case SplineType1D::SPLINE_VEC: break;
      }

      // ---------------------------------------------------------------
      // Creazione e costruzione della spline specifica
      // ---------------------------------------------------------------
      string_view               h = headers[spl];
      std::unique_ptr<Spline> & s = m_splines[spl];

      m_is_monotone[spl] = -1;  // Valore di default (sconosciuto)

      // MIGLIORAMENTO: Questo switch è molto ripetitivo
      // Potrebbe essere refactorizzato con template o funzioni helper
      switch ( stype[spl] )
      {
        case SplineType1D::CONSTANT:
        {
          auto S = std::make_unique<ConstantSpline>( h );
          S->reserve_external( m_npts, m_X, pY );
          S->m_npts = m_npts;
          S->build();
          s = std::move( S );
        }
        break;

        case SplineType1D::LINEAR:
        {
          auto S = std::make_unique<LinearSpline>( h );
          S->reserve_external( m_npts, m_X, pY );
          S->m_npts = m_npts;
          S->build();

          // Check manuale monotonia per spline lineari
          // NOTA: Potrebbe essere spostato in LinearSpline::is_monotone()
          integer flag = 1;  // 1 = strettamente monotona crescente
          for ( integer j = 1; j < m_npts; ++j )
          {
            if ( pY[j - 1] > pY[j] )
            {
              flag = -1;  // Non monotona
              break;
            }
            // MIGLIORAMENTO: La condizione è complicata, meglio separare
            if ( Utils::is_zero( pY[j - 1] - pY[j] ) && m_X[j - 1] < m_X[j] )
              flag = 0;  // Monotona non stretta (plateau)
          }
          m_is_monotone[spl] = flag;
          s                  = std::move( S );
        }
        break;

        // NOTA: Tutti i casi seguenti sono quasi identici
        // Solo il tipo di spline cambia
        case SplineType1D::CUBIC:
        {
          auto S = std::make_unique<CubicSpline>( h );
          S->reserve_external( m_npts, m_X, pY, pYp );
          S->m_npts = m_npts;
          S->build();
          m_is_monotone[spl] = S->is_monotone();
          s                  = std::move( S );
        }
        break;

        case SplineType1D::AKIMA:
        {
          auto S = std::make_unique<AkimaSpline>( h );
          S->reserve_external( m_npts, m_X, pY, pYp );
          S->m_npts = m_npts;
          S->build();
          m_is_monotone[spl] = S->is_monotone();
          s                  = std::move( S );
        }
        break;

        case SplineType1D::VANLEER:
        {
          auto S = std::make_unique<VanLeerSpline>( h );
          S->reserve_external( m_npts, m_X, pY, pYp );
          S->m_npts = m_npts;
          S->build();
          m_is_monotone[spl] = S->is_monotone();
          s                  = std::move( S );
        }
        break;

        case SplineType1D::PCHIP:
        {
          auto S = std::make_unique<PchipSpline>( h );
          S->reserve_external( m_npts, m_X, pY, pYp );
          S->m_npts = m_npts;
          S->build();
          m_is_monotone[spl] = S->is_monotone();
          s                  = std::move( S );
        }
        break;

        case SplineType1D::HERMITE:
        {
          auto S = std::make_unique<HermiteSpline>( h );
          S->reserve_external( m_npts, m_X, pY, pYp );
          S->m_npts = m_npts;
          S->build();
          m_is_monotone[spl] = S->is_monotone();
          s                  = std::move( S );
        }
        break;

        // Casi quintic con seconda derivata
        case SplineType1D::QUINTIC_CUBIC:
        {
          auto S = std::make_unique<QuinticSpline>( Spline_sub_type::CUBIC, h );
          S->reserve_external( m_npts, m_X, pY, pYp, pYpp );
          S->m_npts = m_npts;
          S->build();
          m_is_monotone[spl] = S->is_monotone();
          s                  = std::move( S );
        }
        break;

        case SplineType1D::QUINTIC_AKIMA:
        {
          auto S = std::make_unique<QuinticSpline>( Spline_sub_type::AKIMA, h );
          S->reserve_external( m_npts, m_X, pY, pYp, pYpp );
          S->m_npts = m_npts;
          S->build();
          m_is_monotone[spl] = S->is_monotone();
          s                  = std::move( S );
        }
        break;

        case SplineType1D::QUINTIC_VANLEER:
        {
          auto S = std::make_unique<QuinticSpline>( Spline_sub_type::VANLEER, h );
          S->reserve_external( m_npts, m_X, pY, pYp, pYpp );
          S->m_npts = m_npts;
          S->build();
          m_is_monotone[spl] = S->is_monotone();
          s                  = std::move( S );
        }
        break;

        case SplineType1D::QUINTIC_PCHIP:
        {
          auto S = std::make_unique<QuinticSpline>( Spline_sub_type::PCHIP, h );
          S->reserve_external( m_npts, m_X, pY, pYp, pYpp );
          S->m_npts = m_npts;
          S->build();
          m_is_monotone[spl] = S->is_monotone();
          s                  = std::move( S );
        }
        break;

        case SplineType1D::SPLINE_SET:
        case SplineType1D::SPLINE_VEC:
          UTILS_ERROR(
            "{} At spline n.{} named {}\n"
            "{} not allowed as spline type\n"
            "in SplineSet::build for {}-th spline\n",
            msg,
            spl,
            headers[spl],
            to_string( stype[spl] ),
            spl );
      }

      // ---------------------------------------------------------------
      // Registrazione del mapping nome -> indice
      // ---------------------------------------------------------------
      // ERRORE POTENZIALE: .data() potrebbe non essere sicuro
      // MIGLIORAMENTO: Usare s->name() direttamente se possibile
      // ATTENZIONE: Possibile duplicazione di nomi non gestita
      m_header_to_position.insert( { s->name().data(), static_cast<integer>( spl ) } );
    }

    // ===================================================================
    // FASE 4: Verifica finale allocazione memoria
    // ===================================================================
    // Controllo che tutta la memoria allocata sia stata effettivamente usata
    m_mem.must_be_empty( "SplineSet::build, baseValue" );
    m_mem_p.must_be_empty( "SplineSet::build, basePointer" );
  }

  void SplineSet::eval( real_type const x, GenericContainer & gc ) const
  {
    map_type & vals = gc.set_map();
    for ( auto const & [fst, snd] : m_header_to_position ) vals[fst] = m_splines[snd]->eval( x );
  }

  void SplineSet::eval_D( real_type const x, GenericContainer & gc ) const
  {
    map_type & vals = gc.set_map();
    for ( auto const & [fst, snd] : m_header_to_position ) vals[fst] = m_splines[snd]->eval_D( x );
  }

  void SplineSet::eval_DD( real_type const x, GenericContainer & gc ) const
  {
    map_type & vals = gc.set_map();
    for ( auto const & [fst, snd] : m_header_to_position ) vals[fst] = m_splines[snd]->eval_DD( x );
  }

  void SplineSet::eval_DDD( real_type const x, GenericContainer & gc ) const
  {
    map_type & vals = gc.set_map();
    for ( auto const & [fst, snd] : m_header_to_position ) vals[fst] = m_splines[snd]->eval_DDD( x );
  }

  void SplineSet::eval( vec_real_type const & vec, GenericContainer & gc ) const
  {
    integer const npts = static_cast<integer>( vec.size() );
    map_type &    vals = gc.set_map();

    for ( auto const & [fst, snd] : m_header_to_position )
    {
      vec_real_type & v     = vals[fst].set_vec_real( npts );
      Spline const *  p_spl = m_splines[snd].get();
      for ( integer i = 0; i < npts; ++i ) v[i] = p_spl->eval( vec[i] );
    }
  }

  void SplineSet::eval_D( vec_real_type const & vec, GenericContainer & gc ) const
  {
    integer const npts = static_cast<integer>( vec.size() );
    map_type &    vals = gc.set_map();

    for ( auto const & [fst, snd] : m_header_to_position )
    {
      vec_real_type & v     = vals[fst].set_vec_real( npts );
      Spline const *  p_spl = m_splines[snd].get();
      for ( integer i = 0; i < npts; ++i ) v[i] = p_spl->eval_D( vec[i] );
    }
  }

  void SplineSet::eval_DD( vec_real_type const & vec, GenericContainer & gc ) const
  {
    integer const npts = static_cast<integer>( vec.size() );
    map_type &    vals = gc.set_map();

    for ( auto const & [fst, snd] : m_header_to_position )
    {
      vec_real_type & v     = vals[fst].set_vec_real( npts );
      Spline const *  p_spl = m_splines[snd].get();
      for ( integer i = 0; i < npts; ++i ) v[i] = p_spl->eval_DD( vec[i] );
    }
  }

  void SplineSet::eval_DDD( vec_real_type const & vec, GenericContainer & gc ) const
  {
    integer const npts = static_cast<integer>( vec.size() );
    map_type &    vals = gc.set_map();

    for ( auto const & [fst, snd] : m_header_to_position )
    {
      vec_real_type & v     = vals[fst].set_vec_real( npts );
      Spline const *  p_spl = m_splines[snd].get();
      for ( integer i = 0; i < npts; ++i ) v[i] = p_spl->eval_DDD( vec[i] );
    }
  }

  void SplineSet::eval( real_type const x, vec_string_type const & columns, GenericContainer & gc ) const
  {
    map_type & vals = gc.set_map();
    for ( auto const & S : columns )
    {
      Spline const * p_spl = get_spline( S );
      vals[S]              = p_spl->eval( x );
    }
  }

  void SplineSet::eval_D( real_type const x, vec_string_type const & columns, GenericContainer & gc ) const
  {
    map_type & vals = gc.set_map();
    for ( auto const & S : columns )
    {
      Spline const * p_spl = get_spline( S );
      vals[S]              = p_spl->eval_D( x );
    }
  }

  void SplineSet::eval_DD( real_type const x, vec_string_type const & columns, GenericContainer & gc ) const
  {
    map_type & vals = gc.set_map();
    for ( auto const & S : columns )
    {
      Spline const * p_spl = get_spline( S );
      vals[S]              = p_spl->eval_DD( x );
    }
  }

  void SplineSet::eval_DDD( real_type const x, vec_string_type const & columns, GenericContainer & gc ) const
  {
    map_type & vals = gc.set_map();
    for ( auto const & S : columns )
    {
      Spline const * p_spl = get_spline( S );
      vals[S]              = p_spl->eval_DDD( x );
    }
  }

  void SplineSet::eval( vec_real_type const & vec, vec_string_type const & columns, GenericContainer & gc ) const
  {
    integer const npts = static_cast<integer>( vec.size() );
    map_type &    vals = gc.set_map();

    // Pre-allocation and cache spline pointers
    std::vector<Spline const *> spline_ptrs;
    spline_ptrs.reserve( columns.size() );

    for ( auto const & S : columns )
    {
      spline_ptrs.push_back( get_spline( S ) );
      vals[S].set_vec_real( npts );
    }

    // Evaluate at all points
    for ( integer i = 0; i < npts; ++i )
    {
      for ( size_t j = 0; j < columns.size(); ++j )
      {
        vec_real_type & v = vals[columns[j]].get_vec_real();
        v[i]              = spline_ptrs[j]->eval( vec[i] );
      }
    }
  }

  void SplineSet::eval_D( vec_real_type const & vec, vec_string_type const & columns, GenericContainer & gc ) const
  {
    integer const npts = static_cast<integer>( vec.size() );
    map_type &    vals = gc.set_map();

    // Pre-allocation and cache spline pointers
    std::vector<Spline const *> spline_ptrs;
    spline_ptrs.reserve( columns.size() );

    for ( auto const & S : columns )
    {
      spline_ptrs.push_back( get_spline( S ) );
      vals[S].set_vec_real( npts );
    }

    // Evaluate at all points
    for ( integer i = 0; i < npts; ++i )
    {
      for ( size_t j = 0; j < columns.size(); ++j )
      {
        vec_real_type & v = vals[columns[j]].get_vec_real();
        v[i]              = spline_ptrs[j]->eval_D( vec[i] );
      }
    }
  }

  void SplineSet::eval_DD( vec_real_type const & vec, vec_string_type const & columns, GenericContainer & gc ) const
  {
    integer const npts = static_cast<integer>( vec.size() );
    map_type &    vals = gc.set_map();

    // Pre-allocation and cache spline pointers
    std::vector<Spline const *> spline_ptrs;
    spline_ptrs.reserve( columns.size() );

    for ( auto const & S : columns )
    {
      spline_ptrs.push_back( get_spline( S ) );
      vals[S].set_vec_real( npts );
    }

    // Evaluate at all points
    for ( integer i = 0; i < npts; ++i )
    {
      for ( size_t j = 0; j < columns.size(); ++j )
      {
        vec_real_type & v = vals[columns[j]].get_vec_real();
        v[i]              = spline_ptrs[j]->eval_DD( vec[i] );
      }
    }
  }

  void SplineSet::eval_DDD( vec_real_type const & vec, vec_string_type const & columns, GenericContainer & gc ) const
  {
    integer const npts = static_cast<integer>( vec.size() );
    map_type &    vals = gc.set_map();

    // Pre-allocation and cache spline pointers
    std::vector<Spline const *> spline_ptrs;
    spline_ptrs.reserve( columns.size() );

    for ( auto const & S : columns )
    {
      spline_ptrs.push_back( get_spline( S ) );
      vals[S].set_vec_real( npts );
    }

    // Evaluate at all points
    for ( integer i = 0; i < npts; ++i )
    {
      for ( size_t j = 0; j < columns.size(); ++j )
      {
        vec_real_type & v = vals[columns[j]].get_vec_real();
        v[i]              = spline_ptrs[j]->eval_DDD( vec[i] );
      }
    }
  }

  void SplineSet::eval2( integer const spl, real_type const zeta, real_type & x, real_type vals[], integer const incy )
    const
  {
    intersect( spl, zeta, x );
    size_t ii = 0;
    for ( integer i = 0; i < m_nspl; ++i, ii += incy ) vals[ii] = m_splines[i]->eval( x );
  }

  void SplineSet::eval2_D(
    integer const   spl,
    real_type const zeta,
    real_type &     x,
    real_type       vals[],
    integer const   incy ) const
  {
    Spline const *  S  = intersect( spl, zeta, x );
    real_type const ds = S->D( x );
    size_t          ii = 0;
    for ( integer i = 0; i < m_nspl; ++i, ii += incy ) vals[ii] = m_splines[i]->D( x ) / ds;
  }

  void SplineSet::eval2_DD(
    integer const   spl,
    real_type const zeta,
    real_type &     x,
    real_type       vals[],
    integer const   incy ) const
  {
    Spline const *  S   = intersect( spl, zeta, x );
    real_type const dt  = 1 / S->D( x );
    real_type const dt2 = dt * dt;
    real_type const ddt = -S->DD( x ) * ( dt * dt2 );
    size_t          ii  = 0;
    for ( integer i = 0; i < m_nspl; ++i, ii += incy )
    {
      auto & Si = m_splines[i];
      vals[ii]  = Si->DD( x ) * dt2 + Si->D( x ) * ddt;
    }
  }

  void SplineSet::eval2_DDD(
    integer const   spl,
    real_type const zeta,
    real_type &     x,
    real_type       vals[],
    integer const   incy ) const
  {
    Spline const *  S    = intersect( spl, zeta, x );
    real_type const dt   = 1 / S->D( x );
    real_type const dt3  = dt * dt * dt;
    real_type const ddt  = -S->DD( x ) * dt3;
    real_type const dddt = 3 * ( ddt * ddt ) / dt - S->DDD( x ) * ( dt * dt3 );
    size_t          ii   = 0;
    for ( integer i = 0; i < m_nspl; ++i, ii += incy )
    {
      auto & Si = m_splines[i];
      vals[ii]  = Si->DDD( x ) * dt3 + 3 * Si->DD( x ) * dt * ddt + Si->D( x ) * dddt;
    }
  }

}  // namespace Splines
