/* 
 * File:   LEP2sigmaMu.h
 * Author: mishima
 */

#ifndef LEP2SIGMAMU_H
#define	LEP2SIGMAMU_H

#include "LEP2ThObservable.h"


class LEP2sigmaMu : public LEP2ThObservable {
public:

    /**
     * @brief LEP2sigmaMu constructor
     * @param[in] EW_i an object of EW class
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     */
    LEP2sigmaMu(const EW& EW_i, const double sqrt_s_i) : LEP2ThObservable(EW_i, sqrt_s_i) {
        l_flavor = StandardModel::MU;
    }

    /**
     * @return the cross section for e^+ e^- -> mu^+ mu^- at sqrt_s in pb
     */
    double getThValue();

private:
    
};

#endif	/* LEP2SIGMAMU_H */

