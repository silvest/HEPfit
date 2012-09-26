/* 
 * File:   LEP2AFBmu.h
 * Author: mishima
 */

#ifndef LEP2AFBMU_H
#define	LEP2AFBMU_H

#include <ThObservable.h>
#include "EW.h"


class LEP2AFBmu : public ThObservable {
public:

    /**
     * @brief LEP2AFBmu constructor
     * @param[in] EW_i an object of EW class
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     */
    LEP2AFBmu(const EW& EW_i, const double sqrt_s_i) : ThObservable(EW_i), 
            myEW(EW_i), sqrt_s(sqrt_s_i) {};

    /**
     * @return the forward-backward asymmetry for e^+ e^- -> mu^+ mu^- at sqrt_s
     */
    double getThValue();


private:
    const EW& myEW;
    const double sqrt_s;
};

#endif	/* LEP2AFBMU_H */

