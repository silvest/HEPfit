/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef SUSYSPECTRUM_H
#define	SUSYSPECTRUM_H

#include <gslpp.h>
#include "SUSY.h"

using namespace gslpp;

/**
 * @class SUSYSpectrum
 * @ingroup SUSY
 * @brief A class for calculating the Higgs and sparticle spectra at tree level. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details
 */
class SUSYSpectrum {
public:

    /**
     * @brief A SUSYSpectrum constructor.
     * @param[in] SUSY_in An object of SUSY class.
     */
    SUSYSpectrum(SUSY& SUSY_in);

    /**
     * @brief Computes the Higgs spectrum at tree level.
     */
    void CalcHiggs();

    /**
     * @brief Computes the chargino spectrum at tree level.
     */
    void CalcChargino();

    /**
     * @brief Computes the neutralino spectrum at tree level.
     */
    void CalcNeutralino();

    /**
     * @brief Computes the up-type squark spectrum at tree level.
     */
    void CalcSup();

    /**
     * @brief Computes the down-type squark spectrum at tree level.
     */
    void CalcSdown();

    /**
     * @brief Computes the sneutrino spectrum at tree level.
     */
    void CalcSneutrino();

    /**
     * @brief Computes the charged-slepton spectrum at tree level.
     */
    void CalcSelectron();

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

    matrix<complex> getMchargino() const
    {
        return Mchargino;
    }

    vector<double> getMch() const
    {
        return mch;
    }

    matrix<complex> getU() const
    {
        return U;
    }

    matrix<complex> getV() const
    {
        return V;
    }

    ///////////////////////////////////////////////////////////////////////////

    matrix<complex> getMneutralino() const
    {
        return Mneutralino;
    }

    vector<double> getMneu() const
    {
        return mneu;
    }

    matrix<complex> getN() const
    {
        return N;
    }

    ///////////////////////////////////////////////////////////////////////////

    matrix<complex> getMsup2() const
    {
        return Msup2;
    }

    matrix<complex> getMsdown2() const
    {
        return Msdown2;
    }
    
    vector<double> getMsu2() const
    {
        return m_su2;
    }

    vector<double> getMsd2() const
    {
        return m_sd2;
    }

    matrix<complex> getRu() const
    {
        return Ru;
    }

    matrix<complex> getRd() const
    {
        return Rd;
    }

    ///////////////////////////////////////////////////////////////////////////

    matrix<complex> getMsneutrino2() const
    {
        return Msneutrino2;
    }

    matrix<complex> getMselectron2() const
    {
        return Mselectron2;
    }

    vector<double> getMsn2() const
    {
        return m_sn2;
    }

    vector<double> getMse2() const
    {
        return m_se2;
    }

    matrix<complex> getRn() const
    {
        return Rn;
    }

    matrix<complex> getRl() const
    {
        return Rl;
    }

    ///////////////////////////////////////////////////////////////////////////

private:
    SUSY& mySUSY;

    double mh[4];
    matrix<complex> Mchargino, Mneutralino;
    matrix<complex> Msup2, Msdown2, Msneutrino2, Mselectron2;

    vector<double> mch, mneu, m_su2, m_sd2, m_sn2, m_se2;
    matrix<complex> U, V, N, Ru, Rd, Rn, Rl;

};

#endif	/* SUSYSPECTRUM_H */

