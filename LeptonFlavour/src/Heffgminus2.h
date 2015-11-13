/*
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef HEFFGMINUS2_H
#define	HEFFGMINUS2_H

#include <StandardModel.h>
#include <StandardModelMatching.h>
#include <WilsonCoefficient.h>

using namespace gslpp;

class Heffgminus2 {
public:
    /**
     * @brief constructor
     * @param SM
     * @param modelmatching
     */
    Heffgminus2(const StandardModel & SM_i);
    
    /**
     * 
     * @brief destructor
     */
    virtual ~Heffgminus2();
    
    /**
     * 
     * @param 
     * @param 
     * @return
     */
    vector<complex>** ComputeCoeffgm2mu();

    const StandardModel& GetModel() const {
        return model;
    }
    
private :
    const StandardModel& model;

    WilsonCoefficient coeffgminus2mu;

};

#endif	/* HEFFGMINUS2_H */
