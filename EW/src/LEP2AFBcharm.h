/* 
 * File:   LEP2AFBcharm.h
 * Author: mishima
 */

#ifndef LEP2AFBCHARM_H
#define	LEP2AFBCHARM_H

#include "LEP2ThObservable.h"


class LEP2AFBcharm : public LEP2ThObservable {
public:

    /**
     * @brief LEP2AFBcharm constructor
     * @param[in] EW_i an object of EW class
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     */
    LEP2AFBcharm(const EW& EW_i, const double sqrt_s_i) : LEP2ThObservable(EW_i, sqrt_s_i) {
        q_flavor = StandardModel::CHARM;
    }

    /**
     * @return the forward-backward asymmetry for e^+ e^- -> c cbar at sqrt_s
     */
    double getThValue();

private:

};

#endif	/* LEP2AFBCHARM_H */

