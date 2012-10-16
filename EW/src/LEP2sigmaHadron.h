/* 
 * File:   LEP2sigmaHadron.h
 * Author: mishima
 */

#ifndef LEP2SIGMAHADRON_H
#define	LEP2SIGMAHADRON_H

#include "LEP2ThObservable.h"


class LEP2sigmaHadron : public LEP2ThObservable {
public:

    /**
     * @brief LEP2sigmaHadron constructor
     * @param[in] EW_i an object of EW class
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     */
    LEP2sigmaHadron(const EW& EW_i, const double sqrt_s_i) : LEP2ThObservable(EW_i, sqrt_s_i) {
    }

    /**
     * @return the cross section for e^+ e^- -> hadrons at sqrt_s in pb
     */
    double getThValue();

private:

};

#endif	/* LEP2SIGMAHADRON_H */

