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
 |      Versione Ottimizzata con Eigen - 2026                               |
 |                                                                          |
\*--------------------------------------------------------------------------*/

/**
 * @file SplineCubic.hxx
 * @brief Implementazione di spline cubiche con boundary conditions multiple
 *
 * Questo file contiene l'implementazione di spline cubiche interpolanti
 * utilizzando la libreria Eigen per operazioni vettoriali efficienti.
 * Supporta diverse condizioni al contorno (boundary conditions).
 */

namespace Splines
{

  /**
   * @enum CubicSpline_BC
   * @brief Tipi di condizioni al contorno per spline cubiche
   *
   * Definisce le possibili condizioni al contorno applicabili
   * agli estremi della spline cubica.
   */
  using CubicSpline_BC = enum class CubicSpline_BC : integer {
    EXTRAPOLATE      = 0,  ///< Estrapolazione usando derivate di ordine superiore
    NATURAL          = 1,  ///< Derivata seconda nulla agli estremi
    PARABOLIC_RUNOUT = 2,  ///< Derivata terza nulla agli estremi (continuità cubica)
    NOT_A_KNOT       = 3   ///< Condizione "not-a-knot" di De Boor
  };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  /**
   * @brief Costruisce una spline cubica calcolando le derivate prime
   *
   * Questa funzione risolve il sistema tridiagonale per le derivate seconde (Z)
   * e calcola le derivate prime (Yp) necessarie per l'interpolazione cubica.
   * Utilizza Eigen per operazioni vettoriali efficienti e gestisce automaticamente
   * la memoria.
   *
   * @param[in]  X    Array delle coordinate x (strettamente crescenti)
   * @param[in]  Y    Array delle coordinate y corrispondenti
   * @param[out] Yp   Array delle derivate prime calcolate
   * @param[in]  npts Numero di punti dati
   * @param[in]  bc0  Condizione al contorno iniziale (x = X[0])
   * @param[in]  bcn  Condizione al contorno finale (x = X[npts-1])
   *
   * @note Il sistema tridiagonale viene risolto usando l'algoritmo di Thomas
   * @note Per condizioni NOT_A_KNOT viene usata una correzione iterativa (Sherman-Morrison)
   *
   * @throw UTILS_ASSERT Se npts < 2
   *
   * @par Dettagli matematici
   * La spline cubica su [X[i], X[i+1]] è definita come:
   * \f[
   * S_i(x) = Y[i] + Yp[i](x-X[i]) + \frac{Z[i]}{2}(x-X[i])^2 +
   *          \frac{Z[i+1]-Z[i]}{6h_i}(x-X[i])^3
   * \f]
   * dove h_i = X[i+1] - X[i] e Z[i] è la derivata seconda in X[i].
   */
  void CubicSpline_build(
    real_type const      X[],
    real_type const      Y[],
    real_type            Yp[],
    real_type            Ypp[],
    integer const        npts,
    CubicSpline_BC const bc0,
    CubicSpline_BC const bcn );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  /**
   * @class CubicSpline
   * @brief Classe per la gestione di spline cubiche interpolanti
   *
   * Questa classe eredita da CubicSplineBase e implementa spline cubiche
   * con supporto per diverse condizioni al contorno. La spline risultante
   * garantisce continuità C² (derivate seconde continue) all'interno del
   * dominio, mentre le condizioni al contorno determinano il comportamento
   * agli estremi.
   *
   * @par Caratteristiche
   * - Interpolazione cubica a tratti con continuità C²
   * - 4 tipi di boundary conditions supportate
   * - Gestione automatica di segmenti monotoni multipli
   * - Ottimizzazione con Eigen per efficienza computazionale
   *
   * @par Uso tipico
   * @code
   * CubicSpline spline("my_spline");
   * spline.set_initial_BC(CubicSpline_BC::NATURAL);
   * spline.set_final_BC(CubicSpline_BC::NATURAL);
   * spline.build(x_data, y_data);
   * real_type y = spline.eval(x_query);
   * @endcode
   */
  class CubicSpline final : public CubicSplineBase
  {
  private:
    CubicSpline_BC m_bc0 = CubicSpline_BC::EXTRAPOLATE;  ///< Boundary condition iniziale
    CubicSpline_BC m_bcn = CubicSpline_BC::EXTRAPOLATE;  ///< Boundary condition finale

  public:
    //!
    //! @name Constructors
    //!
    ///@{

    using CubicSplineBase::build;
    using CubicSplineBase::reserve;

    /**
     * @brief Costruisce una spline vuota di tipo CubicSpline
     *
     * @param[in] name Nome identificativo della spline (opzionale)
     */
    explicit CubicSpline( string_view name = "CubicSpline" ) : CubicSplineBase( name ) {}

    /**
     * @brief Distruttore della spline
     */
    ~CubicSpline() override = default;

    ///@}

    //!
    //! @name Setup
    //!
    ///@{

    /**
     * @brief Imposta la condizione al contorno iniziale
     *
     * Definisce come la spline si comporta all'estremo sinistro (x = X[0]).
     * La scelta influenza la curvatura e la smoothness della spline.
     *
     * @param[in] bc0 Tipo di boundary condition iniziale
     *
     * @see CubicSpline_BC per i tipi disponibili
     */
    void set_initial_BC( CubicSpline_BC bc0 ) { m_bc0 = bc0; }

    /**
     * @brief Imposta la condizione al contorno finale
     *
     * Definisce come la spline si comporta all'estremo destro (x = X[npts-1]).
     *
     * @param[in] bcn Tipo di boundary condition finale
     *
     * @see CubicSpline_BC per i tipi disponibili
     */
    void set_final_BC( CubicSpline_BC bcn ) { m_bcn = bcn; }

    // --------------------------- VIRTUALS -----------------------------------

    /**
     * @brief Costruisce la spline a partire dai dati impostati
     *
     * Questo metodo:
     * 1. Verifica la validità dei dati (NaN check)
     * 2. Suddivide i dati in segmenti monotoni crescenti
     * 3. Costruisce una spline cubica per ogni segmento
     * 4. Applica le boundary conditions appropriate
     *
     * @throw UTILS_ASSERT Se npts < 2
     * @throw UTILS_ASSERT Se i dati contengono NaN
     *
     * @note La spline viene costruita anche per dati non monotoni,
     *       dividendoli automaticamente in segmenti monotoni
     */
    void build() override;

    /**
     * @brief Configura la spline usando un GenericContainer
     *
     * Carica i dati e le impostazioni da un container generico.
     *
     * @param[in] gc GenericContainer con i dati e le opzioni
     *
     * @par Campi richiesti nel GenericContainer
     * - `"xdata"`: vettore con le coordinate x
     * - `"ydata"`: vettore con le coordinate y
     *
     * @par Campi opzionali
     * - `"bc_begin"`: boundary condition iniziale ("extrapolate", "natural",
     *                 "parabolic", "not_a_knot")
     * - `"bc_end"`: boundary condition finale (stessi valori di bc_begin)
     *
     * @throw UTILS_ERROR Se i campi richiesti mancano o hanno formato errato
     * @warning Emette warning se campi opzionali mancano (usa default)
     * @warning Emette warning se ci sono campi non riconosciuti
     *
     * @code
     * GenericContainer gc;
     * gc["xdata"] << vector<double>{0, 1, 2, 3};
     * gc["ydata"] << vector<double>{0, 1, 4, 9};
     * gc["bc_begin"] = "natural";
     * gc["bc_end"] = "natural";
     * spline.setup(gc);
     * @endcode
     */
    void setup( GenericContainer const & gc ) override;

    ///@}

    /**
     * @brief Ritorna il tipo di spline
     *
     * @return SplineType1D::CUBIC
     */
    [[nodiscard]] SplineType1D type() const override { return SplineType1D::CUBIC; }

#ifdef SPLINES_BACK_COMPATIBILITY
    /// @deprecated Usa set_initial_BC invece
    void setInitialBC( CubicSpline_BC bc0 ) { m_bc0 = bc0; }

    /// @deprecated Usa set_final_BC invece
    void setFinalBC( CubicSpline_BC bcn ) { m_bcn = bcn; }
#endif
  };

}  // namespace Splines

// EOF: SplineCubic.hxx
