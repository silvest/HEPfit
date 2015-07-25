/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LEP2SIGMAHADRON_H
#define	LEP2SIGMAHADRON_H

#include "LEP2ThObservable.h"

/**
 * @class LEP2sigmaHadron
 * @ingroup EW
 * @brief A class for the hadronic cross section of @f$e^+e^-@f$ annihilations above the @f$Z@f$ pole.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class LEP2sigmaHadron : public LEP2ThObservable {
public:

    /**
     * @brief LEP2sigmaHadron constructor
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     * @param[in] bSigmaForAFB_i true for the denominator of A_FB
     * @param[in] bSigmaForR_i true for the denominator of R_b or R_c
     */
    LEP2sigmaHadron(const StandardModel& SM_i, const double sqrt_s_i,
                    const bool bSigmaForAFB_i=false,
                    const bool bSigmaForR_i=false) 
    : LEP2ThObservable(SM_i, sqrt_s_i, bSigmaForAFB_i, bSigmaForR_i) 
    {
    }

    /**
     * @return the cross section for e^+ e^- -> hadrons at sqrt_s in pb
     */
    double computeThValue();

private:

};

#endif	/* LEP2SIGMAHADRON_H */

