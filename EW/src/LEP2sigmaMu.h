/* 
 * File:   LEP2sigmaMu.h
 * Author: mishima
 */

#ifndef LEP2SIGMAMU_H
#define	LEP2SIGMAMU_H

#include <ThObservable.h>
#include "EW.h"


class LEP2sigmaMu : public ThObservable {
public:

    /**
     * @brief LEP2sigmaMu constructor
     * @param[in] EW_i an object of EW class
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     */
    LEP2sigmaMu(const EW& EW_i, const double sqrt_s_i) : ThObservable(EW_i), 
            myEW(EW_i), sqrt_s(sqrt_s_i) {};

    /**
     * @return the cross section for e^+ e^- -> mu^+ mu^- at sqrt_s in pb
     */
    double getThValue();

    
private:
    const EW& myEW;
    double sqrt_s;
};

#endif	/* LEP2SIGMAMU_H */

