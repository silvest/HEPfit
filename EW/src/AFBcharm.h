/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef AFBCHARM_H
#define	AFBCHARM_H

#include <stdexcept>
#include <ThObservable.h>
#include "EW.h"
#include "EW_CHMN.h"


class AFBcharm : public ThObservable {
public:

    /**
     * @brief AFBcharm constructor
     * @param[in] EW_i an object of EW class
     * @param[in] bCHMN_i true if using EW_CHMN class 
     * @param[in] bBURGESS_i true if using the formula in hep-ph/9411257 by C.P. Burgess
     */
    AFBcharm(const EW& EW_i, const bool bCHMN_i=false, const bool bBURGESS_i=false) : ThObservable(EW_i), 
            myEW(EW_i), myEW_CHMN(EW_i.getSM()), bCHMN(bCHMN_i), bBURGESS(bBURGESS_i) {
        if (bCHMN && bBURGESS)
            throw std::runtime_error("bCHMN and bBURGESS cannot be set to true simultaneously in AFBcharm()");
    };
    /**
     * @return the forward-backward asymmetry of the c-cbar channel
     */
    double getThValue();

    
private:
    const EW& myEW;
    const EW_CHMN myEW_CHMN;
    const bool bCHMN, bBURGESS;
};

#endif	/* AFBCHARM_H */

