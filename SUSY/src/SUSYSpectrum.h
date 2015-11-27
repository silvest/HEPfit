/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef SUSYSPECTRUM_H
#define	SUSYSPECTRUM_H

#include <gslpp.h>

class SUSY;

/**
 * @class SUSYSpectrum
 * @ingroup SUSY
 * @brief A class for calculating the Higgs and sparticle spectra at tree level. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details The SUSYSpectrum class calculates the all the sparticle masses and their mixing matrices at tree-level.
 */
class SUSYSpectrum {
public:

    /**
     * @brief A SUSYSpectrum constructor.
     * @param[in] SUSY_in An object of SUSY class.
     */
    SUSYSpectrum(const SUSY& SUSY_in);

    /**
     * @brief Computes the Higgs spectrum at tree level.
     */
    bool CalcHiggs(double mh[4], gslpp::complex& saeff_i);

    /**
     * @brief Computes the chargino spectrum at tree level.
     */
    bool CalcChargino(gslpp::matrix<gslpp::complex>& U_i, gslpp::matrix<gslpp::complex>& V_i, gslpp::vector<double>& mch_i);

    /**
     * @brief Computes the neutralino spectrum at tree level.
     */
    bool CalcNeutralino(gslpp::matrix<gslpp::complex>& N_i, gslpp::vector<double>& mneu_i);

    /**
     * @brief Computes the up-type squark spectrum at tree level.
     */
    bool CalcSup(gslpp::matrix<gslpp::complex>& Ru_i, gslpp::vector<double>& m_su2_i);

    /**
     * @brief Computes the down-type squark spectrum at tree level.
     */
    bool CalcSdown(gslpp::matrix<gslpp::complex>& Rd_i, gslpp::vector<double>& m_sd2_i);

    /**
     * @brief Computes the sneutrino spectrum at tree level.
     */
    bool CalcSneutrino(gslpp::matrix<gslpp::complex>& Rn_i, gslpp::vector<double>& m_sn2_i);

    /**
     * @brief Computes the charged-slepton spectrum at tree level.
     */
    bool CalcSelectron(gslpp::matrix<gslpp::complex>& Rl_i, gslpp::vector<double>& m_se2_i);

    /**
     * @brief Computes the SUSY spectrum without the Higgs part at tree level.
     */
    bool CalcSpectrum();

    ///////////////////////////////////////////////////////////////////////////

    /**
     * @brief Gets the light Higgs mass.
     * @return The light Higgs mass.
     */
    double getMHl() const
    {
        return mh[0];
    }

    /**
     * @brief Gets the heavy Higgs mass.
     * @return The heavy Higgs mass.
     */
    double getMHh() const
    {
        return mh[1];
    }

    /**
     * @brief Gets the pseudo-scalar Higgs mass.
     * @return The pseudo-scalar Higgs mass.
     */
    double getMHa() const
    {
        return mh[2];
    }

    /**
     * @brief Gets the charged Higgs mass.
     * @return The charged Higgs mass.
     */
    double getMHp() const
    {
        return mh[3];
    }

    ///////////////////////////////////////////////////////////////////////////

    /**
     * @brief Gets the Chargino spectrum at tree-level.
     * @return the Chargino mass eigenvalues
     */
    gslpp::matrix<gslpp::complex> getMchargino() const
    {
        return Mchargino;
    }

    /**
     * @brief Gets the Chargino spectrum at tree-level.
     * @return the Chargino mass eigenvalues
     */
    gslpp::vector<double> getMch() const
    {
        return mch;
    }

    /**
     * @brief Gets the Chargino mixing matrix U.
     * @return the Chargino mixing matrix U
     */
    gslpp::matrix<gslpp::complex> getU() const
    {
        return U;
    }

    /**
     * @brief Gets the Chargino mixing matrix V.
     * @return the Chargino mixing matrix V
     */
    gslpp::matrix<gslpp::complex> getV() const
    {
        return V;
    }

    ///////////////////////////////////////////////////////////////////////////

    /**
     * @brief Gets the Neutralino spectrum at tree-level.
     * @return the Neutralino mass eigenvalues
     */
    gslpp::matrix<gslpp::complex> getMneutralino() const
    {
        return Mneutralino;
    }

    /**
     * @brief Gets the Neutralino spectrum at tree-level.
     * @return the Neutralino mass eigenvalues
     */
    gslpp::vector<double> getMneu() const
    {
        return mneu;
    }

    /**
     * @brief Gets the Neutralino mixing matrix.
     * @return the Neutralino mixing matrix
     */
    gslpp::matrix<gslpp::complex> getN() const
    {
        return N;
    }

    ///////////////////////////////////////////////////////////////////////////

    /**
     * @brief Gets the Up-squark mass matrix at tree-level.
     * @return the Up-squark mass matrix
     */
    gslpp::matrix<gslpp::complex> getMsup2() const
    {
        return Msup2;
    }

    /**
     * @brief Gets the Down-squark mass matrix at tree-level.
     * @return the Down-squark mass matrix
     */
    gslpp::matrix<gslpp::complex> getMsdown2() const
    {
        return Msdown2;
    }
    
    /**
     * @brief Gets the Up-squark spectrum at tree-level.
     * @return the Up-squark mass-squared eigenvalues
     */
    gslpp::vector<double> getMsu2() const
    {
        return m_su2;
    }

    /**
     * @brief Gets the Down-squark spectrum at tree-level.
     * @return the Down-squark mass-squared eigenvalues
     */
    gslpp::vector<double> getMsd2() const
    {
        return m_sd2;
    }

    /**
     * @brief Gets the Up-squark mixing matrix.
     * @return the Up-squark mixing matrix
     */
    gslpp::matrix<gslpp::complex> getRu() const
    {
        return Ru;
    }

    /**
     * @brief Gets the Down-squark mixing matrix.
     * @return the Down-squark mixing matrix
     */
    gslpp::matrix<gslpp::complex> getRd() const
    {
        return Rd;
    }

    ///////////////////////////////////////////////////////////////////////////
    

    /**
     * @brief Gets the Sneutrino mass matrix at tree-level.
     * @return the Sneutrino mass matrix
     */
    gslpp::matrix<gslpp::complex> getMsneutrino2() const
    {
        return Msneutrino2;
    }

    /**
     * @brief Gets the Slepton mass matrix at tree-level.
     * @return the Slepton mass matrix
     */
    gslpp::matrix<gslpp::complex> getMselectron2() const
    {
        return Mselectron2;
    }

    /**
     * @brief Gets the Sneutrino spectrum at tree-level.
     * @return the Sneutrino mass-squared eigenvalues
     */
    gslpp::vector<double> getMsn2() const
    {
        return m_sn2;
    }

    /**
     * @brief Gets the Slepton spectrum at tree-level.
     * @return the Slepton mass-squared eigenvalues
     */
    gslpp::vector<double> getMse2() const
    {
        return m_se2;
    }

    /**
     * @brief Gets the Sneutrino mixing matrix.
     * @return the Sneutrino mixing matrix
     */
    gslpp::matrix<gslpp::complex> getRn() const
    {
        return Rn;
    }

    /**
     * @brief Gets the Slepton mixing matrix.
     * @return the Slepton mixing matrix
     */
    gslpp::matrix<gslpp::complex> getRl() const
    {
        return Rl;
    }
    
    /**
     * Sort sfermion masses in increasing order.
     * @param[in,out] m_sf2 A vector of sfermion mass squared.
     * @param[in,out] Rf The corresponding rotation matrix.
     */
    void SortSfermionMasses(gslpp::vector<double>& m_sf2, gslpp::matrix<gslpp::complex>& Rf) const;

    ///////////////////////////////////////////////////////////////////////////

private:
    const SUSY& mySUSY;

    /**
     * @brief Stores the tree-level Higgs spectrum. 
     */
    double mh[4];

    /**
     * @brief Stores the Sine of tree-level  CP-even mixing angle.
     */
    gslpp::complex saeff;

    /**
     * @brief Stores the tree-level Chargino and Neutralino mass matrix.
     */
    gslpp::matrix<gslpp::complex> Mchargino, Mneutralino;
    
    /**
     * @brief Stores the tree-level Chargino and Neutralino mixing matrices.
     */
    gslpp::matrix<gslpp::complex> U, V, N;

    /**
     * @brief Stores the tree-level Up-squark, Down-squark, Sneutrino, and Slepton mass matrix.
     */
    gslpp::matrix<gslpp::complex> Msup2, Msdown2, Msneutrino2, Mselectron2;

    /**
     * @brief Stores the tree-level Up-squark, Down-squark, Sneutrino, and Slepton mass-squared eigenvalues.
     */
    gslpp::vector<double> mch, mneu, m_su2, m_sd2, m_sn2, m_se2;

    /**
     * @brief Stores the tree-level Up-squark, Down-squark, Sneutrino, and Slepton mixing matrices.
     */
    gslpp::matrix<gslpp::complex> Ru, Rd, Rn, Rl;

};

#endif	/* SUSYSPECTRUM_H */

