/* 
 * File:   sin2thetaEff.h
 * Author: mishima
 *
 * Created on June 9, 2011, 3:45 PM
 */

#ifndef SIN2THETAEFF_H
#define	SIN2THETAEFF_H

#include <ThObservable.h>
#include "EW.h"

class sin2thetaEff : public ThObservable {
public:

    /**
     * @brief sin2thetaEff constructor
     * @param[in] myEW an object of EW class
     */
    sin2thetaEff(const EW& myEW);

    /**
     * @return the effective weak mixing angle for a leptonic channel
     */
    double getThValue();

private:

};

#endif	/* SIN2THETAEFF_H */

