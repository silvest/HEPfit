/*
 * Copyright (C) 2012 SusyFit Collaboration
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
    li_lj_gamma(LeptonFlavour& LeptonFlavour_i);
    
    /**
     *
     * @brief 
     * @return
     */
    double computeThValue();
    
    
protected:
    
private:
    LeptonFlavour& myLeptonFlavour;
    
};

#endif	/* LI_LJ_GAMMA_H */

