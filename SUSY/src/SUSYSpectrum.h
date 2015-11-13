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
 * @details
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

    gslpp::matrix<gslpp::complex> getMchargino() const
    {
        return Mchargino;
    }

    gslpp::vector<double> getMch() const
    {
        return mch;
    }

    gslpp::matrix<gslpp::complex> getU() const
    {
        return U;
    }

    gslpp::matrix<gslpp::complex> getV() const
    {
        return V;
    }

    ///////////////////////////////////////////////////////////////////////////

    gslpp::matrix<gslpp::complex> getMneutralino() const
    {
        return Mneutralino;
    }

    gslpp::vector<double> getMneu() const
    {
        return mneu;
    }

    gslpp::matrix<gslpp::complex> getN() const
    {
        return N;
    }

    ///////////////////////////////////////////////////////////////////////////

    gslpp::matrix<gslpp::complex> getMsup2() const
    {
        return Msup2;
    }

    gslpp::matrix<gslpp::complex> getMsdown2() const
    {
        return Msdown2;
    }
    
    gslpp::vector<double> getMsu2() const
    {
        return m_su2;
    }

    gslpp::vector<double> getMsd2() const
    {
        return m_sd2;
    }

    gslpp::matrix<gslpp::complex> getRu() const
    {
        return Ru;
    }

    gslpp::matrix<gslpp::complex> getRd() const
    {
        return Rd;
    }

    ///////////////////////////////////////////////////////////////////////////

    gslpp::matrix<gslpp::complex> getMsneutrino2() const
    {
        return Msneutrino2;
    }

    gslpp::matrix<gslpp::complex> getMselectron2() const
    {
        return Mselectron2;
    }

    gslpp::vector<double> getMsn2() const
    {
        return m_sn2;
    }

    gslpp::vector<double> getMse2() const
    {
        return m_se2;
    }

    gslpp::matrix<gslpp::complex> getRn() const
    {
        return Rn;
    }

    gslpp::matrix<gslpp::complex> getRl() const
    {
        return Rl;
    }

    ///////////////////////////////////////////////////////////////////////////

private:
    const SUSY& mySUSY;

    double mh[4];
    gslpp::complex saeff;
    gslpp::matrix<gslpp::complex> Mchargino, Mneutralino;
    gslpp::matrix<gslpp::complex> Msup2, Msdown2, Msneutrino2, Mselectron2;

    gslpp::vector<double> mch, mneu, m_su2, m_sd2, m_sn2, m_se2;
    gslpp::matrix<gslpp::complex> U, V, N, Ru, Rd, Rn, Rl;

};

#endif	/* SUSYSPECTRUM_H */

