/*
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GMINUS2_H
#define	GMINUS2_H

#include "gslpp.h"
#include "ThObservable.h"
#include "LeptonFlavour.h"

class gminus2 : public ThObservable {
public:
    /**
     * constructor
     * @param LeptonFlavour
     */
    gminus2(const StandardModel& SM_i);

    /**
     *
     * @brief 
     * @return
     */
    double computeThValue();

protected:

private:
    
};

class gminus2_mu : public gminus2 {
public:
    
    /**
     * @brief Constructor.
     */
    gminus2_mu(const StandardModel& SM_i);
    
    /**
     * @return gminus2_mu
     */
    double computeThValue();
    
private:
    const StandardModel& mySM;

};

#endif	/* GMINUS2_H */
