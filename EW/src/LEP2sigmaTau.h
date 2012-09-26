/* 
 * File:   LEP2sigmaTau.h
 * Author: mishima
 */

#ifndef LEP2SIGMATAU_H
#define	LEP2SIGMATAU_H

#include <ThObservable.h>
#include "EW.h"


class LEP2sigmaTau : public ThObservable {
public:

    /**
     * @brief LEP2sigmaTau constructor
     * @param[in] EW_i an object of EW class
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     */
    LEP2sigmaTau(const EW& EW_i, const double sqrt_s_i) : ThObservable(EW_i), 
            myEW(EW_i), sqrt_s(sqrt_s_i) {};

    /**
     * @return the cross section for e^+ e^- -> tau^+ tau^- at sqrt_s in pb
     */
    double getThValue();


private:
    const EW& myEW;
    const double sqrt_s;
};

#endif	/* LEP2SIGMATAU_H */

