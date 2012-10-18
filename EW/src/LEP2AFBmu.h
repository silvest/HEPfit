/* 
 * File:   LEP2AFBmu.h
 * Author: mishima
 */

#ifndef LEP2AFBMU_H
#define	LEP2AFBMU_H

#include "LEP2ThObservable.h"


class LEP2AFBmu : public LEP2ThObservable {
public:

    /**
     * @brief LEP2AFBmu constructor
     * @param[in] EW_i an object of EW class
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     */
    LEP2AFBmu(const EW& EW_i, const double sqrt_s_i) : LEP2ThObservable(EW_i, sqrt_s_i) {
        l_flavor = StandardModel::MU;
    }

    /**
     * @return the forward-backward asymmetry for e^+ e^- -> mu^+ mu^- at sqrt_s
     */
    double getThValue();

private:
    
};

#endif	/* LEP2AFBMU_H */

