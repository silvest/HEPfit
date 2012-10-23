/* 
 * File:   LEP2sigmaTau.h
 * Author: mishima
 */

#ifndef LEP2SIGMATAU_H
#define	LEP2SIGMATAU_H

#include "LEP2ThObservable.h"


class LEP2sigmaTau : public LEP2ThObservable {
public:

    /**
     * @brief LEP2sigmaTau constructor
     * @param[in] EW_i an object of EW class
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     * @param[in] bSigmaForAFB_i true for the denominator of A_FB
     */
    LEP2sigmaTau(const EW& EW_i, const double sqrt_s_i, 
                 const bool bSigmaForAFB_i=false) : LEP2ThObservable(EW_i, sqrt_s_i, bSigmaForAFB_i) {
        l_flavor = StandardModel::TAU;
    }

    /**
     * @return the cross section for e^+ e^- -> tau^+ tau^- at sqrt_s in pb
     */
    double getThValue();

private:

};

#endif	/* LEP2SIGMATAU_H */

