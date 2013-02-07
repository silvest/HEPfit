/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef HEFFDF1_H
#define	HEFFDF1_H

#include <StandardModel.h>
#include <StandardModelMatching.h>
#include <WilsonCoefficient.h>
#include "EvolDF1bsg.h"
#include <sstream>

using namespace gslpp;

class HeffDF1bsg {
public:
    /**
     * @brief constructor
     * @param SM
     * @param SM_Matching
     */
    HeffDF1bsg(const StandardModel & SM, StandardModelMatching & SM_Matching);
    
    /**
     * 
     * @brief destructor
     */
    virtual ~HeffDF1bsg();
    
    /**
     * 
     * @param mu is the low energy scale
     * @param scheme indicates the renormalization scheme
     * @return the effective hamiltonian at the scale mu B -> s gamma decay, Misiak basis, Chetyrkin et al hep-ph/9612313
     */
    vector<complex>** ComputeCoeffBsg(double mu, schemes scheme = NDR);
    
    
    EvolDF1bsg getUDF1() const {
        return u;
    }

    const StandardModel& GetModel() const {
        return model;
    }
    
private :
    const StandardModel& model;
    ModelMatching& modelmatching;
    WilsonCoefficient coeffbsg;
    EvolDF1bsg u;
    
};

#endif	/* HEFFDF1_H */

