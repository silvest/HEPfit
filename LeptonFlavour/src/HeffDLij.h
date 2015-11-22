/*
 * Copyright (C) 2013 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef HEFFDLIJ_H
#define	HEFFDLIJ_H

#include "StandardModel.h"
#include "StandardModelMatching.h"
#include "WilsonCoefficient.h"

class HeffDLij {
public:
    /**
     * @brief constructor
     * @param SM
     * @param modelmatching
     */
    HeffDLij(const StandardModel & SM_i);
    
    /**
     * 
     * @brief destructor
     */
    virtual ~HeffDLij();
    
    /**
     * 
     * @param 
     * @param 
     * @return
     */
    gslpp::vector<gslpp::complex>** ComputeCoeffDLij(int li_lj);

    const StandardModel& GetModel() const {
        return model;
    }
    
private :
    const StandardModel& model;

    WilsonCoefficient coeffDLij_1;
    WilsonCoefficient coeffDLij_2;
    WilsonCoefficient coeffDLij_3;

    //gslpp::vector<gslpp::complex> nlep, nlep2, nlepCC;
};

#endif	/* HEFFDLIJ_H */

