/* 
 * Copyright (C) 2017 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GENERALTHDMMATCHING_H
#define	GENERALTHDMMATCHING_H

#include "StandardModelMatching.h"

class GeneralTHDM;

/**
 * @class GeneralTHDMMatching
 * @ingroup GeneralTHDM
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class GeneralTHDMMatching : public StandardModelMatching {
public:
    GeneralTHDMMatching(const GeneralTHDM & GeneralTHDM_i);

    /**
     * @brief Wilson coefficient for \f$ (g-2)_{\mu} \f$.
     * @return
     */
    virtual std::vector<WilsonCoefficient>& CMgminus2mu();

    /** Calculates the muon g-2 at LO**/
    /**
     * @brief Calculates amplitudes for \f$ (g-2)_{\mu} \f$ at one loop from \cite Broggio:2014mna.
     * @return 
     */    
    virtual double gminus2muLO();

    /** Calculates the NLO contribution to the muon g-2**/
    /**
     * @brief Calculates amplitudes for \f$ (g-2)_{\mu} \f$ at approximate two-loop from \cite Broggio:2014mna.
     * @return 
     */    
    virtual double gminus2muNLO();

    /**
     *
     * @brief Updates to new GTHDM parameter set.
     * @return
     */
    void updateGTHDMParameters();

private:
    const GeneralTHDM & myGTHDM;

    WilsonCoefficient mcgminus2mu;
    
    double GF, mMU;

};

#endif	/* GENERALTHDMMATCHING_H */
