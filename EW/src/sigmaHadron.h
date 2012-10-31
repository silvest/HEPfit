/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef SIGMAHADRON_H
#define	SIGMAHADRON_H

#include <ThObservable.h>
#include "EW.h"
#include "EW_CHMN.h"


class sigmaHadron : public ThObservable {
public:

    /**
     * @brief sigmaHadron constructor
     * @param[in] EW_i an object of EW class
     * @param[in] bCHMN_i true if using EW_CHMN class 
     */
    sigmaHadron(const EW& EW_i, const bool bCHMN_i=false) : ThObservable(EW_i), 
            myEW(EW_i), myEW_CHMN(EW_i.getSM()), bCHMN(bCHMN_i) {};

    /**
     * @return the hadronic cross section in nb
     */
    double getThValue();


private:
    const EW& myEW;
    const EW_CHMN myEW_CHMN;
    const bool bCHMN;
};

#endif	/* SIGMAHADRON_H */

