/*
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef HEFFDL1_H
#define	HEFFDL1_H

#include <StandardModel.h>
#include <StandardModelMatching.h>
#include <WilsonCoefficient.h>

using namespace gslpp;

class HeffDL1 {
public:
    /**
     * @brief constructor
     * @param SM
     * @param modelmatching
     */
    HeffDL1(const StandardModel & SM_i);
    
    /**
     * 
     * @brief destructor
     */
    virtual ~HeffDL1();
    
    /**
     * 
     * @param 
     * @param 
     * @return
     */
    vector<complex>** ComputeCoeffDL1();
    
    
    WilsonCoefficient getCoeffDL1() const {
        return coeffDL1;
    }
    

    const StandardModel& GetModel() const {
        return model;
    }
    
private :
    const StandardModel& model;
    
    WilsonCoefficient coeffDL1;
    
    //gslpp::vector<complex> nlep, nlep2, nlepCC;
};

#endif	/* HEFFDL1_H */

