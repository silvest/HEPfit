/* 
 * File:   sin2thetaEff.h
 * Author: mishima
 */

#ifndef SIN2THETAEFF_H
#define	SIN2THETAEFF_H

#include <ThObservable.h>
#include "EW.h"
#include "EW_CHMN.h"


class sin2thetaEff : public ThObservable {
public:

    /**
     * @brief sin2thetaEff constructor
     * @param[in] EW_i an object of EW class
     * @param[in] bCHMN_i true if using EW_CHMN class 
     */
    sin2thetaEff(const EW& EW_i, const bool bCHMN_i=false) : ThObservable(EW_i), 
            myEW(EW_i), myEW_CHMN(EW_i.getSM()), bCHMN(bCHMN_i) {};

    /**
     * @return the effective weak mixing angle for a leptonic channel
     */
    double getThValue();

    
private:
    const EW& myEW;
    const EW_CHMN myEW_CHMN;
    const bool bCHMN;
};

#endif	/* SIN2THETAEFF_H */

