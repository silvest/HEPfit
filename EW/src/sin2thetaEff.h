/* 
 * File:   sin2thetaEff.h
 * Author: mishima
 */

#ifndef SIN2THETAEFF_H
#define	SIN2THETAEFF_H

#include <ThObservable.h>
#include "EW.h"


class sin2thetaEff : public ThObservable {
public:

    /**
     * @brief sin2thetaEff constructor
     * @param[in] EW_i an object of EW class
     */
    sin2thetaEff(const EW& EW_i);

    /**
     * @return the effective weak mixing angle for a leptonic channel
     */
    double getThValue();

    
private:
    double sin2_theta_eff;
    
};

#endif	/* SIN2THETAEFF_H */

