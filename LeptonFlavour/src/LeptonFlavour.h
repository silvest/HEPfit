/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LEPTONFLAVOUR_H
#define	LEPTONFLAVOUR_H

#include <ThObsType.h>
#include <StandardModel.h>
//#include "HeffDF2.h"

/**
 * @addtogroup LeptonFlavour
 * @brief A project for Lepton Flavour observables.
 * @{
 */

using namespace gslpp;

class LeptonFlavour : public ThObsType {
public:

    LeptonFlavour(const StandardModel& SM_i) : ThObsType(SM_i) 
            //, HDF2(SM_i), HDS1(SM_i), HDB1(SM_i) 
    {   
        if(!SM_i.IsModelInitialized())
            throw std::runtime_error("Model not initialized "); 
    };

    /*const HeffDF2& getHDF2() const {
        return HDF2;
    }*/
    
    /*vector<complex>** ComputeCoeffBd(double mu, schemes scheme = NDR) {
        return HDF2.ComputeCoeffBd(mu, scheme);
    }*/
    
private:
    //HeffDF2 HDF2; 
};

/**
 * @}
 */

#endif	/* LEPTONFLAVOUR_H */

