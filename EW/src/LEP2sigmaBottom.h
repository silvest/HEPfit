/* 
 * File:   LEP2sigmaBottom.h
 * Author: mishima
 */

#ifndef LEP2SIGMABOTTOM_H
#define	LEP2SIGMABOTTOM_H

#include "LEP2ThObservable.h"


class LEP2sigmaBottom : public LEP2ThObservable {
public:

    /**
     * @brief LEP2sigmaBottom constructor
     * @param[in] EW_i an object of EW class
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     * @param[in] bSigmaForAFB_i true for the denominator of A_FB
     */
    LEP2sigmaBottom(const EW& EW_i, const double sqrt_s_i, 
                    const bool bSigmaForAFB_i=false) : LEP2ThObservable(EW_i, sqrt_s_i, bSigmaForAFB_i) {
        q_flavor = StandardModel::BOTTOM;
    }

    /**
     * @return the cross section for e^+ e^- -> b bbar at sqrt_s in pb
     */
    double getThValue();

private:
    
};

#endif	/* LEP2SIGMABOTTOM_H */

