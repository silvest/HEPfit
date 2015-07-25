/*
 * Copyright (C) 2013 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LI_LJ_GAMMA_H
#define	LI_LJ_GAMMA_H

#include <ThObservable.h>
#include "LeptonFlavour.h"

class li_lj_gamma : public ThObservable {
public:
    /**
     * constructor
     * @param LeptonFlavour
     */
    li_lj_gamma(const StandardModel& SM_i);
    
    /**
     *
     * @brief 
     * @return
     */
    double computeThValue();
    
    
protected:
    
private:
    const StandardModel& mySM;
    
};

#endif	/* LI_LJ_GAMMA_H */

