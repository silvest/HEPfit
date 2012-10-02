/* 
 * File:   AFBcharm.h
 * Author: mishima
 */

#ifndef AFBCHARM_H
#define	AFBCHARM_H

#include <ThObservable.h>
#include "EW.h"
#include "EW_CHMN.h"


class AFBcharm : public ThObservable {
public:

    /**
     * @brief AFBcharm constructor
     * @param[in] EW_i an object of EW class
     * @param[in] bCHMN_i true if using EW_CHMN class 
     */
    AFBcharm(const EW& EW_i, const bool bCHMN_i=false) : ThObservable(EW_i), 
            myEW(EW_i), myEW_CHMN(EW_i.getSM()), bCHMN(bCHMN_i) {};

    /**
     * @return the forward-backward asymmetry of the c-cbar channel
     */
    double getThValue();

    
private:
    const EW& myEW;
    const EW_CHMN myEW_CHMN;
    const bool bCHMN;
};

#endif	/* AFBCHARM_H */

