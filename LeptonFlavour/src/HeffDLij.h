/*
 * Copyright (C) 2013 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef HEFFDLIJ_H
#define	HEFFDLIJ_H

#include <StandardModel.h>
#include <StandardModelMatching.h>
#include <WilsonCoefficient.h>

using namespace gslpp;

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
    vector<complex>** ComputeCoeffDLij(int li_lj);

    const StandardModel& GetModel() const {
        return model;
    }
    
private :
    const StandardModel& model;

    WilsonCoefficient coeffDLij_1;
    WilsonCoefficient coeffDLij_2;
    WilsonCoefficient coeffDLij_3;

    //gslpp::vector<complex> nlep, nlep2, nlepCC;
};

#endif	/* HEFFDLIJ_H */

