/* 
 * File:   LEP2AFBtau.h
 * Author: mishima
 */

#ifndef LEP2AFBTAU_H
#define	LEP2AFBTAU_H

#include "LEP2ThObservable.h"


class LEP2AFBtau : public LEP2ThObservable {
public:

    /**
     * @brief LEP2AFBtau constructor
     * @param[in] EW_i an object of EW class
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     */
    LEP2AFBtau(const EW& EW_i, const double sqrt_s_i) : LEP2ThObservable(EW_i, sqrt_s_i) {
    }

    /**
     * @return the forward-backward asymmetry for e^+ e^- -> tau^+ tau^- at sqrt_s
     */
    double getThValue();

private:

};

#endif	/* LEP2AFBTAU_H */

